# -*- coding: utf-8 -*-

import requests
import pubchempy as pcp
from core.tools import convert_all_multiple
import re

# Extracts the CAS registry number from a list of compound synonyms
def extract_cas(synonyms):
    for s in synonyms:
        if re.match(r'^\d{2,7}-\d{2}-\d$', s):  
            return s
    return None

def fetch_cas_from_xrefs(cid):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/RegistryID/JSON"
    try:
        response = requests.get(url, timeout=5)
        if response.status_code == 200:
            data = response.json()
            registry_ids = data.get("InformationList", {}).get("Information", [{}])[0].get("RegistryID", [])
            
            for rid in registry_ids:
                if re.match(r'^\d{2,7}-\d{2}-\d$', rid):
                    return rid
    except Exception as e:
        print(f"Error retrieving PubChem XRefs : {e}")
    return None

# Recursively searches PubChem JSON data to extract strings matching a specific heading
def extract_pubchem_property(data, property_name):
    results = []

    def walk(obj):
        if isinstance(obj, dict):
            if obj.get("TOCHeading", "").strip().lower() == property_name.lower():
                infos = obj.get("Information", [])
                for item in infos:
                    try:
                        val = item["Value"]["StringWithMarkup"][0]["String"]
                        results.append(val)
                    except:
                        pass

            for v in obj.values():
                walk(v)

        elif isinstance(obj, list):
            for i in obj:
                walk(i)

    walk(data)
    return results

# Retrieves and filters a specific property dataset ensuring mandatory keys are present
def get_property(data, property_name, required_keys, default_keys={}):
    raw = extract_pubchem_property(data, property_name)
    
    converted = []
    for v in raw:
        multiple = convert_all_multiple(v)
        if multiple:
            converted.extend(multiple)
    
    converted = [v for v in converted if v is not None]
    
    result = []
    for v in converted:
        if not all(k in v for k in required_keys):
            continue
        for k, default in default_keys.items():
            if k not in v:
                v[k] = default
        result.append(v)
    
    return result

# Extracts properties that include specific temperature fields via Regex
def get_property_with_temp(data, property_name):
    # Retrieve raw entries for the specified property
    raw_entries = extract_pubchem_property(data, property_name)
    results = []
    
    # Regex to capture a value (e.g. 25.5) and a temperature (e.g. 25 degC or 300 K)
    pattern = r"([\d\.]+)\s*(?:kJ/mol|kcal/mol).*(?:at|@)?\s*([\d\.]+)\s*(degC|K|C)"
    
    for entry in raw_entries:
        match = re.search(pattern, entry, re.IGNORECASE)
        if match:
            val, temp, unit = match.groups()
            results.append({
                'value': float(val),
                'temp': float(temp),
                'unit': unit
            })
    return results

# Parses GHS (security) classification data to extract signal words, h-codes, and p-codes
def get_ghs(data):
    raw = extract_pubchem_property(data, "GHS Classification")
    
    signal_words = set()
    h_codes = {}
    p_codes = set()
    
    for text in raw:
        if text.strip() in ['Danger', 'Warning']:
            signal_words.add(text.strip())
        
        h_matches = re.findall(r'(H\d{3}(?:\+H\d{3})*)[:\s]+([^\[]+)', text)
        for code, description in h_matches:
            h_codes[code] = description.strip()
        
        p_matches = re.findall(r'P\d{3}(?:\+P\d{3})*', text)
        for code in p_matches:
            p_codes.add(code)
    
    return {
        'signal_word' : list(signal_words),
        'H_codes' : h_codes,
        'P_codes' : sorted(list(p_codes))}

# Selects the data entry whose experimental temperature is closest to the target
def get_closest_to_temperature(data_list, target_T_celsius, value_key, temp_key='degC'):
    if not data_list:
        return None
    
    with_temp = [v for v in data_list if temp_key in v and v[temp_key] is not None]
    
    if not with_temp:
        return data_list[0] if data_list else None
    
    closest = min(with_temp, 
                  key=lambda v: abs(round(v[temp_key], 1) - round(target_T_celsius, 1)))
    return closest

# Fixes formatting errors in chemical formulas involving lowercase letters after numbers
def fix_formula(formula):
    formula = re.sub(r'(\d)([a-z])', 
                      lambda m: m.group(1) + m.group(2).upper(), formula)
    return formula

# Main interface function to query PubChem API and populate the molecule object
def fill_from_pubchem(mol) -> None:
    # 1 Perform a single lookup of the compound on PubChem
    # Accept either a molecular name/formula or a bare numeric CID
    if mol.name.strip().isdigit():
        try:
            c = pcp.Compound.from_cid(int(mol.name.strip()))
        except Exception as e:
            raise ValueError(f"PubChem service unreachable: {str(e)}")
        if c is None or c.cid is None:
            raise ValueError(f"[PubChem] CID '{mol.name}' not found")
    else:
        try:
            compounds = pcp.get_compounds(mol.name, "name")
        except Exception as e:
            raise ValueError(f"PubChem service unreachable: {str(e)}")
        # 2 Raise an error immediately if no compound is found
        if not compounds:
            raise ValueError(f"[PubChem] molecule '{mol.name}' not found")
        c = compounds[0]
    
    # Assign molecular identifiers
    mol.id_pubchem = c.cid
    synonyms = c.synonyms
    mol.formula = fix_formula(c.molecular_formula)
    mol.CAS = extract_cas(synonyms)
    mol.SMILES = c.smiles
    mol.iupac_name = c.iupac_name

    # If the user searched by CID, resolve mol.name to the canonical name so
    # downstream parsers (NIST, CDMS, etc.) receive a human-readable identifier
    if mol.name.strip().isdigit():
        mol.name = mol.iupac_name or (synonyms[0] if synonyms else mol.name)

    print(f"[PubChem] OK : CAS={mol.CAS} CID={mol.id_pubchem}")

    mol.molar_mass.pubchem = float(c.molecular_weight) if c.molecular_weight else None
    mol.logp.pubchem = float(c.xlogp) if c.xlogp else None

    # Fetch extended structural properties via the PUG View JSON API
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{mol.id_pubchem}/JSON"
    data = requests.get(url).json()

    # Boiling Point extraction
    bp = get_property(
        data,
        "Boiling Point",
        required_keys=['degC'],
        default_keys={'pascal': 101325.0})
    
    # Melting Point extraction
    mp = get_property(
        data,
        "Melting Point",
        required_keys=['degC'],
        default_keys={'pascal': 101325.0})

    if bp:
        bp_temps = [b['degC'] for b in bp]
        mp = [v for v in mp if not any(abs(v['degC'] - t) < 1.0 for t in bp_temps)]
    
    # Set the target ambient temperature calculation in Celsius
    T_celsius = mol.temperature - 273.15

    # Vapor Pressure extraction
    raw_vp = extract_pubchem_property(data, "Vapor Pressure")
    cleaned_vp = []
    for v in raw_vp:
        v_clean = v.replace('[', '').replace(']', '')
        v_clean = re.sub(r'[xX]10\s*\+?(-?\d+)', r'e\1', v_clean)
        cleaned_vp.append(v_clean)
    converted_vp = []
    for v in cleaned_vp:
        multiple = convert_all_multiple(v)
        if multiple:
            converted_vp.extend(multiple)       
    vp = [v for v in converted_vp if v is not None and all(k in v for k in ['pascal', 'degC'])]
    vp_best = get_closest_to_temperature(vp, T_celsius, 'pascal')
    
    # Vaporization Enthalpy extraction
    raw_hv = extract_pubchem_property(data, "Heat of Vaporization")
    cleaned_hv = []
    for v in raw_hv:
        v_clean = v.replace('[', '').replace(']', '')
        v_clean = re.sub(r'[xX]10\s*\+?(-?\d+)', r'e\1', v_clean)
        cleaned_hv.append(v_clean)
    converted_hv = []
    for v in cleaned_hv:
        multiple = convert_all_multiple(v)
        if multiple:
            converted_hv.extend(multiple)       
    hv = [v for v in converted_hv if v is not None and 'kilojoule / mole' in v]
    if hv:
        hv_best = get_closest_to_temperature(hv, T_celsius, 'kilojoule / mole', temp_key='degC')
        mol.vaporizationenthalpy.pubchem = float(hv_best['kilojoule / mole'])        

        if 'degC' in hv_best:
            mol.vaporizationenthalpy.pubchem_temperature = float(hv_best['degC'] + 273.15) 
        else:
            mol.vaporizationenthalpy.pubchem_temperature = 298.15
    
    # Fusion Enthalpy extraction
    raw_hf = extract_pubchem_property(data, "Heat of Fusion")
    cleaned_hf = []
    for v in raw_hf:
        v_clean = v.replace('[', '').replace(']', '')
        v_clean = re.sub(r'[xX]10\s*\+?(-?\d+)', r'e\1', v_clean)
        cleaned_hf.append(v_clean)
    converted_hf = []
    for v in cleaned_hf:
        multiple = convert_all_multiple(v)
        if multiple:
            converted_hf.extend(multiple)       
    hf = [v for v in converted_hf if v is not None and 'kilojoule / mole' in v]
    if hf:
        hf_best = get_closest_to_temperature(hf, T_celsius, 'kilojoule / mole', temp_key='degC')
        mol.fusionenthalpy.pubchem = float(hf_best['kilojoule / mole'])        

        if 'degC' in hf_best:
            mol.fusionenthalpy.pubchem_temperature = float(hf_best['degC'] + 273.15) 
        else:
            mol.fusionenthalpy.pubchem_temperature = 298.15
    
    # Dipole Moment extraction
    dm = get_property(data, "Dipole Moment", required_keys=['debye'])

    # Map remaining extracted endpoints to molecule object
    mol.boilingtemp.pubchem = float(bp[0]['degC']) if bp else None
    mol.meltingtemp.pubchem = float(mp[0]['degC']) if mp else None
    mol.vaporpressure.pubchem = float(vp_best['pascal']) if vp_best else None
    mol.vaporpressure.pubchem_temperature = float(vp_best['degC'] + 273.15) if vp_best else None
    mol.dipolemoment.pubchem = float(dm[0]['debye']) if dm else None
    
    # GHS Hazard Profiles assignment
    ghs = get_ghs(data)
    mol.security.signal_word = ghs['signal_word']
    mol.security.H_codes = ghs['H_codes']
    mol.security.P_codes = ghs['P_codes']