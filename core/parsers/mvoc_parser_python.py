# -*- coding: utf-8 -*-

import subprocess
import sys
import json
from core.tools import preprocess, convert_value

# Autonomous script to simulate browser search via Playwright
MVOC_SCRIPT = """
import asyncio
from playwright.async_api import async_playwright
from bs4 import BeautifulSoup
import re
import json

# Automates a headless chromium browser to fetch mvoc search results for a pubchem id
async def search_mvoc_async(pubchem_id):
    async with async_playwright() as p:
        browser = await p.chromium.launch(headless=True)
        page = await browser.new_page()
        await page.goto("https://bioinformatics.charite.de/mvoc/index.php?site=mvocsearch")
        await page.fill('input[name=\\"compound_id\\"]', str(pubchem_id))
        await page.click('input[name=\\"B3\\"]')
        await page.wait_for_load_state('load')
        await page.wait_for_timeout(2000)
        html = await page.content()
        browser.close()
    return html

html = asyncio.run(search_mvoc_async({pubchem_id}))
soup = BeautifulSoup(html, 'html.parser')
soup_text = soup.get_text()

result = {{}}

rows = soup.find_all('tr')
for row in rows:
    cells = row.find_all(['td', 'th'])
    if len(cells) >= 2:
        key = cells[0].text.strip()
        val = cells[1].text.strip()
        if 'LogP' in key : result['logP'] = val
        if 'IUPAC' in key : result['iupac'] = val
        if 'SMILES' in key : result['smiles'] = val
        if 'Molweight' in key : result['Molweight'] = val
        if 'organism' in key.lower(): result['in_organism'] = val

henry_idx = html.lower().find("henry")
if henry_idx > 0:
    henry_text = html[henry_idx:henry_idx+1000]
    
    match = re.search(r'([\d.]+[xX]10[-\d]+)\s*atm', henry_text)
    if match:
        result['henry_raw'] = match.group(1)
    else:
        match_simple = re.search(r'estimated\s*as\s*([\d.]+)\s*atm', henry_text, re.IGNORECASE)
        if match_simple:
            result['henry_raw'] = match_simple.group(1)

for table in soup.find_all('table'):
    text = table.get_text()
    if 'boiling' in text.lower():
        match = re.search(r'(-?\\d+(?:\\.\\d+)?)\\s*\\u00b0C', text)
        if match:
            result['boiling_point_C'] = float(match.group(1))
            break

vp_found = False
for table in soup.find_all('table'):
    headers = [th.text.strip() for th in table.find_all('th')]
    if 'Pressure' in headers and 'Reference' in headers:
        rows_vp = table.find_all('tr')[1:]  
        for row in rows_vp:
            cols = [c.text.strip() for c in row.find_all('td')]
            if not cols or not cols[0]:
                continue
            match = re.search(
                r'([\d.]+)\s*(mm\s*Hg|Pa|atm|kPa)\s*at\s*([\d.]+)\s*(deg\s*C|°C|K)',
                cols[0],
                re.IGNORECASE)
            if match:
                result['vapor_pressure_value'] = float(match.group(1))
                result['vapor_pressure_unit']  = match.group(2).strip()
                unit_T = match.group(4).strip().lower()
                T_val  = float(match.group(3))
                if 'k' in unit_T:
                    result['vapor_pressure_temp_K'] = T_val 
                else:
                    result['vapor_pressure_temp_K'] = T_val+273.15  
                vp_found = True
                break
        if vp_found:
            break

if not vp_found:
    # Scan text for specific pattern format (e.g. 4.66X10+5 mm Hg)
    # Handles both scientific notation formats: X10+5 and e+05
    vp_pattern = r'vapor\s*pressure[,\s]+([\d.]+(?:[xX]10[+-]?\d+|e[+-]?\d+)?)\s*(mm\s*Hg|Pa|atm|kPa)'
    vp_match = re.search(vp_pattern, soup_text, re.IGNORECASE)
    
    if vp_match:
        val_str = vp_match.group(1).lower()
        # Convert scientific format token "X10+5" into valid python float expression
        if 'x10' in val_str:
            val_str = val_str.replace('x10', 'e')
        
        result['vapor_pressure_value'] = float(val_str)
        result['vapor_pressure_unit'] = vp_match.group(2).strip()
        # Set default value fallback to 25°C for document database entry points
        result['vapor_pressure_temp_K'] = 298.15
        vp_found = True

print(json.dumps(result, ensure_ascii=False))
"""

# Executes the background Playwright scraper and maps retrieved mVOC properties onto the molecule
def fill_from_mvoc(mol):
    if mol.id_pubchem is None:
        print("[mVOC] CID PubChem unavailable")
        return

    script = MVOC_SCRIPT.format(pubchem_id=mol.id_pubchem)

    result = subprocess.run(
        [sys.executable, '-c', script],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"[mVOC] Subprocess error: {result.stderr}")
    
    if not result.stdout.strip():
        return "[WARN] [mVOC]: Molecule not found"
        
    try:
        data = json.loads(result.stdout.strip().split('\n')[-1])
    except json.JSONDecodeError:
        raise ValueError("[mVOC] Invalid JSON format")
    
    if not data:
        print("[mVOC]: Molecule not found in database")
        return "[WARN] [mVOC]: No data available for this molecule"

    # Property mappings
    mol.molar_mass.mvoc = float(data['Molweight']) if 'Molweight' in data else None   
    mol.logp.mvoc = float(data['logP']) if 'logP' in data else None
    mol.boilingtemp.mvoc = float(data['boiling_point_C']) if 'boiling_point_C' in data else None

    if data.get('henry_raw'):
        raw_val = preprocess(data['henry_raw'])
        hvpc = float(raw_val)
        mol.henryconstant.mvoc = 1.83089 / (101325 * hvpc) # Convert to atm-1
        mol.henryconstant.mvoc_temperature = 298.15 if mol.henryconstant.mvoc is not None else None

    if 'vapor_pressure_value' in data:
        mol.vaporpressure.mvoc = convert_value(data['vapor_pressure_value'], 'mmHg', 'Pa')
        mol.vaporpressure.mvoc_temperature = data['vapor_pressure_temp_K']
    
    print("[mVOC]: OK")
    return "[OK] [mVOC] data loaded"