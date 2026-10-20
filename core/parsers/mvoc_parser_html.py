# -*- coding: utf-8 -*-

import sys
import json
import re
from bs4 import BeautifulSoup
from playwright.sync_api import sync_playwright
from core.tools import preprocess, convert_value

# Fetch mVOC result page HTML using Playwright
def _fetch_mvoc_html(pubchem_id):
    with sync_playwright() as p:
        # Launch headless Chromium
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()
        # Open mVOC search page
        page.goto(
            "https://bioinformatics.charite.de/mvoc/index.php?site=mvocsearch"
        )
        # Submit PubChem CID
        page.fill('input[name="compound_id"]', str(pubchem_id))
        page.click('input[name="B3"]')
        # Wait for results page
        page.wait_for_load_state("load")
        page.wait_for_timeout(2000)
        # Get rendered HTML
        html = page.content()
        browser.close()
    return html

# Load mVOC properties into molecule object
def fill_from_mvoc(mol):

    if mol.id_pubchem is None:
        print("[mVOC] CID PubChem unavailable")
        return "[WARN] [mVOC]: CID PubChem unavailable"

    # Retrieve HTML from mVOC
    try:
        html = _fetch_mvoc_html(mol.id_pubchem)

    except Exception as e:
        raise RuntimeError(
            f"[mVOC] Playwright navigation error: {e}")

    # Parse HTML
    soup = BeautifulSoup(html, "html.parser")
    soup_text = soup.get_text()

    result = {}

    # Extract general properties from tables
    rows = soup.find_all("tr")

    for row in rows:

        cells = row.find_all(["td", "th"])

        if len(cells) >= 2:

            key = cells[0].text.strip()
            val = cells[1].text.strip()

            if "LogP" in key:
                result["logP"] = val

            if "IUPAC" in key:
                result["iupac"] = val

            if "SMILES" in key:
                result["smiles"] = val

            if "Molweight" in key:
                result["Molweight"] = val

            if "organism" in key.lower():
                result["in_organism"] = val

    # Extract Henry constant
    henry_idx = html.lower().find("henry")

    if henry_idx > 0:

        henry_text = html[henry_idx : henry_idx + 1000]

        match = re.search(
            r"([\d.]+[xX]10[-\d]+)\s*atm", henry_text)

        if match:
            result["henry_raw"] = match.group(1)

        else:
            match_simple = re.search(
                r"estimated\s*as\s*([\d.]+)\s*atm", henry_text, re.IGNORECASE)
            if match_simple:
                result["henry_raw"] = match_simple.group(1)

    # Extract boiling point
    for table in soup.find_all("table"):
        text = table.get_text()
        if "boiling" in text.lower():
            match = re.search(r"(-?\d+(?:\.\d+)?)\s*\u00b0C", text)

            if match:
                result["boiling_point_C"] = float(match.group(1))
                break
    
    for table in soup.find_all("table"):
        text = table.get_text()
        if "melting" in text.lower():
            match = re.search(r"(-?\d+(?:\.\d+)?)\s*\u00b0C", text)

            if match:
                result["melting_point_C"] = float(match.group(1))
                break

    # Extract vapor pressure from dedicated table
    vp_found = False

    for table in soup.find_all("table"):

        headers = [th.text.strip() for th in table.find_all("th")]

        if "Pressure" in headers and "Reference" in headers:
            rows_vp = table.find_all("tr")[1:]
            for row in rows_vp:
                cols = [c.text.strip() for c in row.find_all("td")]
                if not cols or not cols[0]:
                    continue
                match = re.search(
                    r"([\d.]+)\s*(mm\s*Hg|Pa|atm|kPa)\s*at\s*([\d.]+)\s*(deg\s*C|°C|K)", cols[0], re.IGNORECASE,)

                if match:
                    result["vapor_pressure_value"] = float(match.group(1))
                    result["vapor_pressure_unit"] = match.group(2).strip()

                    unit_T = match.group(4).strip().lower()
                    T_val = float(match.group(3))

                    if "k" in unit_T:
                        result["vapor_pressure_temp_K"] = T_val
                    else:
                        result["vapor_pressure_temp_K"] = T_val + 273.15

                    vp_found = True
                    break

            if vp_found:
                break

    # Fallback vapor pressure extraction
    if not vp_found:
        vp_pattern = (r"vapor\s*pressure[,\s]+"
            r"([\d.]+(?:[xX]10[+-]?\d+|e[+-]?\d+)?)\s*"
            r"(mm\s*Hg|Pa|atm|kPa)")

        vp_match = re.search(vp_pattern, soup_text, re.IGNORECASE)

        if vp_match:
            val_str = vp_match.group(1).lower()

            if "x10" in val_str:
                val_str = val_str.replace("x10", "e")

            result["vapor_pressure_value"] = float(val_str)
            result["vapor_pressure_unit"] = vp_match.group(2).strip()
            result["vapor_pressure_temp_K"] = 298.15

            vp_found = True

    if not result:

        print("[mVOC]: Molecule not found in database")
        return "[WARN] [mVOC]: No data available for this molecule"

    # Populate molecule properties
    mol.molar_mass.mvoc = (
        float(result["Molweight"])
        if "Molweight" in result
        else None)

    mol.logp.mvoc = (
        float(result["logP"])
        if "logP" in result
        else None)

    mol.boilingtemp.mvoc = (
        float(result["boiling_point_C"])
        if "boiling_point_C" in result
        else None)

    mol.meltingtemp.mvoc = (
        float(result["melting_point_C"])
        if "melting_point_C" in result
        else None)

    if result.get("henry_raw"):

        raw_val = preprocess(result["henry_raw"])
        hvpc = float(raw_val)

        mol.henryconstant.mvoc = (1.83089 / (101325 * hvpc))

        mol.henryconstant.mvoc_temperature = (298.15 if mol.henryconstant.mvoc is not None else None)

    if "vapor_pressure_value" in result:

        mol.vaporpressure.mvoc = convert_value(result["vapor_pressure_value"],"mmHg","Pa",)

        mol.vaporpressure.mvoc_temperature = (result["vapor_pressure_temp_K"])

    print("[mVOC]: OK")

    return "[OK] [mVOC] data loaded"