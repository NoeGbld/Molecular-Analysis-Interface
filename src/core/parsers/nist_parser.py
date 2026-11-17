# -*- coding: utf-8 -*-

import re
import requests
import numpy as np
from bs4 import BeautifulSoup
from core.tools import convert_value, henry_at_T

BASE_URL = "https://webbook.nist.gov/cgi/cbook.cgi"
HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36"
    )
}

session = requests.Session()


def fetch_nist_html(mol, mask):
    # Fetches HTML content from the NIST WebBook database for a specific molecule
    params = {"Units": "SI"}

    if hasattr(mol, 'CAS') and mol.CAS:
        params["ID"] = str(mol.CAS).strip()
    else:
        params["Name"] = mol.name

    if mask and str(mask).strip():
        params["Mask"] = str(mask).strip()

    try:
        response = session.get(
            BASE_URL,
            params=params,
            headers=HEADERS,
            timeout=10,
        )
        response.raise_for_status()

        if "Species List" in response.text:
            soup = BeautifulSoup(response.text, "html.parser")
            link = soup.find("a", href=re.compile(r"cbook\.cgi\?ID="))

            if not link:
                return None

            extra_mask = f"&Mask={mask}" if mask else ""
            url = f"https://webbook.nist.gov{link['href']}&Units=SI{extra_mask}"
            response = session.get(url, headers=HEADERS, timeout=10)
            return response.text

        if "did not match any species" in response.text or "not found" in response.text.lower():
            return None

        return response.text

    except requests.RequestException:
        return None


def parse_nist_tables(html: str, mol_temperature: float, data_dict: dict | None = None) -> dict:
    # Parses NIST WebBook HTML tables to extract chemical and thermodynamic properties
    default_structure = {
        "antoine_params": [], "henry_k0": None, "henry_alpha": None, 
        "dh_vap_list": [], "dh_fus_list": [], "molar_mass": None, 
        "CAS": None, "boiling_point_K": None, "melting_point_K": None, "dipole_moment": None
    }
    
    if not html:
        return data_dict or default_structure

    soup = BeautifulSoup(html, "html.parser")

    if data_dict is None:
        data_dict = default_structure

    for key in ["molar_mass", "CAS", "boiling_point_K", "melting_point_K", "dipole_moment"]:
        if key not in data_dict:
            data_dict[key] = None

    for key in ["dh_vap_A", "dh_vap_alpha", "dh_vap_beta", "dh_vap_Tc", "dh_vap_tmin", "dh_vap_tmax"]:
        if key not in data_dict:
            data_dict[key] = None
            
    if "dh_vap_list" not in data_dict: data_dict["dh_vap_list"] = []
    if "dh_fus_list" not in data_dict: data_dict["dh_fus_list"] = []

    for element in soup.find_all(["li", "p"]):
        text = element.get_text(" ", strip=True).lower()
        
        if ("poids" in text and "moléculaire" in text) or "molecular weight" in text:
            m_mass = re.search(r"([\d.]+)", text)
            if m_mass and data_dict["molar_mass"] is None:
                data_dict["molar_mass"] = float(m_mass.group(1))

    # Step 1: Table identification
    tables = soup.find_all("table")
    is_ln_form = False
    
    for table in tables:
        table_text = table.get_text().lower()
        
        if "critical pressure" in table_text or "triple point pressure" in table_text:
            continue
            
        is_vap_table = any(x in table_text for x in ["vap h", "vaporization enthalpy", "vaph", "enthalpy of vaporization", "sublimation enthalpy", "sub h"]) or \
                       (any(x in table_text for x in ["a (kj/mol)", "α", "β"]) and "temperature" in table_text)
                       
        is_fus_table = any(x in table_text for x in ["fus h", "fusion enthalpy", "fush", "enthalpy of fusion"])

        # Specific header check to reliably identify an Antoine parameters table
        headers = [th.get_text().strip().lower() for th in table.find_all("th")]
        is_antoine_table = ("temperature (k)" in headers and "a" in headers and "b" in headers and "c" in headers)

        # Step 2: Data extraction based on table type
        if is_vap_table or is_fus_table or is_antoine_table:
            if is_antoine_table:
                # Check logarithmic base form (ln vs log10) in the preceding section text
                prev_text = ""
                for sibling in table.previous_siblings:
                    if sibling.name in ["h2", "h3", "p", "div"]:
                        prev_text = sibling.get_text().lower()
                        break
                if "ln(p)" in prev_text or "ln(p)" in table_text:
                    is_ln_form = True

            rows = table.find_all("tr")
            for row in rows:
                cells = [c.get_text(" ", strip=True) for c in row.find_all(["td", "th"])]
                if len(cells) < 2:
                    continue
                
                label = cells[0].lower().strip()
                
                # Case 1: Row contains a temperature range formatted as "T_min to T_max"
                if "to" in label and len(cells) >= 4:
                    try:
                        clean_range_text = re.sub(r'\s+', ' ', label.replace('\xa0', ' ')).strip()
                        t_range = clean_range_text.split(" to ")
                        if len(t_range) == 2:
                            t_min = float(re.search(r"([\d.]+)", t_range[0]).group(1).strip("."))
                            t_max = float(re.search(r"([\d.]+)", t_range[1]).group(1).strip("."))
                            
                            values = []
                            for c in cells[1:4]:
                                m = re.search(r"(-?[\d.]+)", c.replace(" ", ""))
                                if m: 
                                    values.append(float(m.group(1)))
                            
                            if len(values) >= 3:
                                if is_antoine_table:
                                    record = {
                                        "t_min": t_min, "t_max": t_max,
                                        "A": values[0], "B": values[1], "C": values[2],
                                        "is_ln": is_ln_form
                                    }
                                    if record not in data_dict["antoine_params"]:
                                        data_dict["antoine_params"].append(record)
                                elif is_vap_table and data_dict["dh_vap_A"] is None:
                                    # Horizontal enthalpy coefficients format
                                    data_dict["dh_vap_tmin"] = t_min
                                    data_dict["dh_vap_tmax"] = t_max
                                    data_dict["dh_vap_A"] = values[0]
                                    data_dict["dh_vap_beta"] = values[1]
                                    data_dict["dh_vap_Tc"] = values[2]
                    except Exception:
                        pass
                    continue
                
                # Case 2: Vertical format for classic enthalpy variables
                val_text = "".join(cells[1:]).replace(" ", "").replace(",", ".").strip()
                m_val = re.search(r"(-?[\d.]+)", val_text)
                
                if m_val and not is_antoine_table:
                    val = float(m_val.group(1).strip("."))
                    if is_vap_table:
                        if "tc" in label or "t c" in label or "t_c" in label or "critical" in label:
                            data_dict["dh_vap_Tc"] = val
                        elif "temperature" in label or "range" in label:
                            t_match = re.findall(r"([\d.]+)", val_text)
                            if len(t_match) >= 2:
                                if data_dict["dh_vap_tmin"] is None:
                                    data_dict["dh_vap_tmin"] = float(t_match[0].strip("."))
                                    data_dict["dh_vap_tmax"] = float(t_match[1].strip("."))
                        elif label == "a" or "a (" in label or "a =" in label:
                            data_dict["dh_vap_A"] = val
                        elif label in ["α", "alpha"] or "α" in label or "alpha" in label:
                            data_dict["dh_vap_alpha"] = val
                        elif label in ["β", "beta"] or "β" in label or "beta" in label:
                            data_dict["dh_vap_beta"] = val

                # Case 3: Single point experimental data tables (value vs temperature)
                if re.match(r"^\d", cells[0].strip()) and not is_antoine_table:
                    if not any(x in label for x in ["a", "α", "alpha", "β", "beta", "t", "critical"]):
                        try:
                            value_clean = re.search(r"(-?[\d.]+)", cells[0].replace(" ", "")).group(1)
                            temp_clean = re.search(r"(-?[\d.]+)", cells[1].replace(" ", "")).group(1)
                            v_data = float(value_clean.strip("."))
                            t_data = float(temp_clean.strip("."))
                            
                            if is_vap_table and (t_data, v_data) not in data_dict["dh_vap_list"]:
                                data_dict["dh_vap_list"].append((t_data, v_data))
                            elif is_fus_table and (t_data, v_data) not in data_dict["dh_fus_list"]:
                                data_dict["dh_fus_list"].append((t_data, v_data))
                        except Exception:
                            pass

    # Step 3: Other properties (boiling point, dipole, Henry's law)
    for row in soup.find_all("tr"):
        cells = [c.get_text(" ", strip=True) for c in row.find_all(["td", "th"])]
        if len(cells) < 2:
            continue
        row_text = " ".join(cells).lower()
        parent_table = row.find_parent("table")
        table_text = parent_table.get_text().lower() if parent_table else ""

        if "boil" in row_text or "sublimation point" in row_text or "triple point" in row_text:
            if data_dict["boiling_point_K"] is None:
                m = re.search(r"([\d.]+)", cells[1].replace(" ", ""))
                if m: 
                    data_dict["boiling_point_K"] = float(m.group(1))
        
        elif "fus" in row_text or "melting point" in row_text or "fusion temperature" in row_text:
            if data_dict["melting_point_K"] is None:
                m = re.search(r"([\d.]+)", cells[1].replace(" ", ""))
                if m: 
                    data_dict["melting_point_K"] = float(m.group(1))
            
        elif "dipole" in row_text:
            m = re.search(r"([\d.]+)", cells[1])
            if m: 
                data_dict["dipole_moment"] = float(m.group(1))

        elif data_dict["henry_k0"] is None:
            if "henry" in table_text or "mol/(kg" in row_text:
                try:
                    m1 = re.search(r"([\d.]+)", cells[0])
                    if m1:
                        data_dict["henry_k0"] = float(m1.group(1))
                        m2 = re.search(r"([\d.]+)", cells[1]) if len(cells) > 1 else None
                        data_dict["henry_alpha"] = float(m2.group(1)) if m2 else None
                except Exception:
                    pass

    return data_dict


def fill_from_nist(mol) -> str:
    data = {
        "antoine_params": [], "henry_k0": None, "henry_alpha": None,
        "dh_vap_list": [], "dh_fus_list": [], "molar_mass": None,
        "CAS": None, "boiling_point_K": None, "melting_point_K": None, "dipole_moment": None
    }

    html_retrieved = False

    for mask in ("", "4", "10"):
        html = fetch_nist_html(mol, mask=mask)
        if html:
            html_retrieved = True
            data = parse_nist_tables(
                html,
                mol_temperature=mol.temperature,
                data_dict=data,
            )
    
    data_found = any([
        data.get("antoine_params"),
        data.get("henry_k0") is not None,
        data.get("dh_vap_list"),
        data.get("dh_fus_list"),
        data.get("boiling_point_K") is not None,
        data.get("melting_point_K") is not None,
        data.get("dipole_moment") is not None,
        data.get("dh_vap_A") is not None
    ])

    if not html_retrieved or not data_found:
        print(f"[NIST]: NO DATA FOUND FOR {mol.name}")
        return f"[WARN] [NIST] : No data available for this molecule"

    if not hasattr(mol, 'molar_mass') or mol.molar_mass is None:
        mol.molar_mass = {}
    if isinstance(mol.molar_mass, dict):
        mol.molar_mass["nist"] = data.get("molar_mass")
    else:
        mol.molar_mass.nist = data.get("molar_mass")

    if not hasattr(mol, 'boilingtemp') or mol.boilingtemp is None:
        mol.boilingtemp = {}
    converted_bp = convert_value(data["boiling_point_K"], "K", "degC") if data.get("boiling_point_K") else None
    if isinstance(mol.boilingtemp, dict):
        mol.boilingtemp["nist"] = converted_bp
    else:
        mol.boilingtemp.nist = converted_bp
        
    if not hasattr(mol, 'meltingtemp') or mol.meltingtemp is None:
        mol.meltingtemp = {}
    converted_mp = convert_value(data["melting_point_K"], "K", "degC") if data.get("melting_point_K") else None
    if isinstance(mol.meltingtemp, dict):
        mol.meltingtemp["nist"] = converted_mp
    else:
        mol.meltingtemp.nist = converted_mp

    if not hasattr(mol, 'dipolemoment') or mol.dipolemoment is None:
        mol.dipolemoment = {}
    if isinstance(mol.dipolemoment, dict):
        mol.dipolemoment["nist"] = data.get("dipole_moment")
    else:
        mol.dipolemoment.nist = data.get("dipole_moment")

    if not hasattr(mol, 'fusionenthalpy') or mol.fusionenthalpy is None:
        mol.fusionenthalpy = {}
        
    fus_val, fus_temp = None, None
    if data.get("dh_fus_list"):
        closest_fus = min(data["dh_fus_list"], key=lambda x: abs(x[0] - mol.temperature))
        fus_val = closest_fus[1]
        fus_temp = closest_fus[0]

    if isinstance(mol.fusionenthalpy, dict):
        mol.fusionenthalpy["nist"] = fus_val
        mol.fusionenthalpy["nist_temperature"] = fus_temp
    else:
        mol.fusionenthalpy.nist = fus_val
        mol.fusionenthalpy.nist_temperature = fus_temp

    if not hasattr(mol, 'vaporizationenthalpy') or mol.vaporizationenthalpy is None:
        mol.vaporizationenthalpy = {}

    has_coeffs = (data.get("dh_vap_A") is not None and 
                  data.get("dh_vap_beta") is not None and 
                  data.get("dh_vap_Tc") is not None)

    in_range = True
    if data.get("dh_vap_tmin") is not None and data.get("dh_vap_tmax") is not None:
        in_range = (data["dh_vap_tmin"] <= mol.temperature <= data["dh_vap_tmax"])

    h_vap_final = None
    t_vap_final = None
    
    if has_coeffs and in_range:
        A = data["dh_vap_A"]
        beta = data["dh_vap_beta"]
        Tc = data["dh_vap_Tc"]
        alpha = data.get("dh_vap_alpha") if data.get("dh_vap_alpha") is not None else 0.0
        
        Tr = mol.temperature / Tc
        if 0 < Tr < 1:
            h_vap_final = float(A * np.exp(-alpha * Tr) * ((1.0 - Tr) ** beta))
            t_vap_final = mol.temperature

    if h_vap_final is None and data.get("dh_vap_list"):
        closest_vap = min(data["dh_vap_list"], key=lambda x: abs(x[0] - mol.temperature))
        h_vap_final = closest_vap[1]
        t_vap_final = closest_vap[0]

    if isinstance(mol.vaporizationenthalpy, dict):
        mol.vaporizationenthalpy["nist"] = h_vap_final
        mol.vaporizationenthalpy["nist_temperature"] = t_vap_final
    else:
        mol.vaporizationenthalpy.nist = h_vap_final
        mol.vaporizationenthalpy.nist_temperature = t_vap_final

    if not hasattr(mol, 'vaporpressure') or mol.vaporpressure is None:
        mol.vaporpressure = {}

    v_pres = None
    t_pres = None
    if data["antoine_params"]:
        T = mol.temperature
        perfect_match = [p for p in data["antoine_params"] if p["t_min"] <= T <= p["t_max"]]
        
        if perfect_match:
            p = perfect_match[0]
            t_pres = mol.temperature
        else:
            p = min(data["antoine_params"], key=lambda x: abs(((x["t_min"] + x["t_max"]) / 2.0) - T))
            if mol.temperature > p["t_max"]:
                T = p["t_max"]
            elif mol.temperature < p["t_min"]:
                T = p["t_min"]
            t_pres = T

        try:
            denominator = T + p["C"]
            if abs(denominator) > 1e-5:
                exponent = p["A"] - (p["B"] / denominator)
                if p.get("is_ln", False):
                    p_bar = np.exp(exponent)
                else:
                    p_bar = 10 ** exponent
                v_pres = p_bar * 1e5
        except (ZeroDivisionError, OverflowError):
            v_pres = None

    if isinstance(mol.vaporpressure, dict):
        mol.vaporpressure["nist"] = v_pres
        mol.vaporpressure["nist_temperature"] = t_pres
    else:
        mol.vaporpressure.nist = v_pres
        mol.vaporpressure.nist_temperature = t_pres

    if not hasattr(mol, 'henryconstant') or mol.henryconstant is None:
        mol.henryconstant = {}

    h_val = None
    h_temp = None
    if data["henry_k0"] is not None:
        H0 = data["henry_k0"]
        alpha = data["henry_alpha"]
        T = mol.temperature
        h_temp = T
        try:
            H_T = henry_at_T(H0, alpha, T) if alpha else H0
            h_val = H_T * 1825.4 * 1e-5
        except Exception:
            h_val = H0 * 1825.4 * 1e-5

    if isinstance(mol.henryconstant, dict):
        mol.henryconstant["nist"] = h_val
        mol.henryconstant["nist_temperature"] = h_temp
    else:
        mol.henryconstant.nist = h_val
        mol.henryconstant.nist_temperature = h_temp
        
    print("[NIST]: OK")
    return "[OK] [NIST] : Data successfully retrieved"