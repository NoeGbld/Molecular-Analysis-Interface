# -*- coding: utf-8 -*-

import requests
from bs4 import BeautifulSoup
from core.tools import henry_at_T, preprocess

# Scrapes Henry's law constants and parameters from the Sander online database using a CAS number
def get_henry_from_sander(cas):
    url = f"https://www.henrys-law.org/henry/casrn/{cas}"
    headers = {'User-Agent': 'Mozilla/5.0'}
    
    r = requests.get(url, headers=headers)
    if r.status_code != 200:
        print(f"[Sander] Error code: {r.status_code}")
        return None
    
    soup = BeautifulSoup(r.text, 'html.parser')
    tables = soup.find_all('table')
    
    if len(tables) < 5:
        print("[Sander] Unexpected HTML structure")
        return None
    
    table = tables[4]
    rows = table.find_all('tr')
    
    results = []
    for row in rows[1:]:
        cols = [c.text.strip() for c in row.find_all('td')]
        cols = [c.replace('\xa0', ' ') for c in cols]
        if cols and cols[0]:   
            try:
                results.append({
                    'Hscp': float(preprocess(cols[0])) if cols[0] else None,
                    'dlnHs_dT': float(preprocess(cols[1])) if len(cols) > 1 and cols[1] else None,
                    'reference': cols[2] if len(cols) > 2 else None,
                    'type': cols[3] if len(cols) > 3 else None,
                    'notes': cols[4] if len(cols) > 4 else None,
                })
            except Exception:
                continue
    return results if results else None

# Main interface function to extract Sander datasets and compute temperature-adjusted Henry constants
def fill_from_sander(mol):
    if mol.CAS is None:
        raise ValueError("CAS number missing: cannot query Sander.")

    results = get_henry_from_sander(mol.CAS)

    if not results:
        return "[WARN] No data found in Sander for this CAS."

    # Priority 1: Search for dataset containing both the constant and its temperature coefficient
    best = next((d for d in results if d.get("dlnHs_dT") is not None), None)
    
    # Priority 2: Fallback to dataset with only the raw constant available
    if best is None:
        best = next((d for d in results if d.get("Hscp") is not None), None)

    # Convert to atm-1 via Sander's tables and adjust temperature parameters
    if best["dlnHs_dT"] is not None:
        mol.henryconstant.sander = henry_at_T(best["Hscp"], best["dlnHs_dT"], mol.temperature) * 1.83089
        mol.henryconstant.sander_temperature = mol.temperature
    else:
        mol.henryconstant.sander = best["Hscp"] * 1.83089
        mol.henryconstant.sander_temperature = 298.15

    mol.henryconstant.sander_ref = best.get("reference")
    mol.henryconstant.sander_all = results

    print("[Sander]: OK")
    return "[OK] Sander data loaded."