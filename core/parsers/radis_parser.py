# -*- coding: utf-8 -*-

import re
import os
import sys
import threading
import time
from collections import Counter
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import numpy as np
from scipy import constants
import pubchempy as pcp

# ── Constante de conversion HITRAN → nm²·MHz ─────────────────────────────────
C_HITRAN_TO_NM2MHZ = 2.99792458e18   # S[nm²·MHz] = S[cm/molec] × c[cm·MHz]

# ── Numéros HITRAN (isotopologue principal) ───────────────────────────────────
HITRAN_MOLECULES = {
    "H2O":1,"CO2":2,"O3":3,"N2O":4,"CO":5,"CH4":6,"O2":7,"NO":8,
    "SO2":9,"NO2":10,"NH3":11,"HNO3":12,"OH":13,"HF":14,"HCL":15,
    "HBR":16,"HI":17,"CLO":18,"OCS":19,"H2CO":20,"HOCL":21,"N2":22,
    "HCN":23,"CH3CL":24,"H2O2":25,"C2H2":26,"C2H6":27,"PH3":28,
    "COF2":29,"SF6":30,"H2S":31,"HCOOH":32,"HO2":33,"O":34,
    "CLONO2":35,"NO+":36,"HOBR":37,"C2H4":38,"CH3OH":39,"CH3BR":40,
    "CH3CN":41,"CF4":42,"C4H2":43,"HC3N":44,"H2":45,"CS":46,
    "SO3":47,"C2N2":48,"COCL2":49,"SO":50,"C3H4":51,"CH3":52,"CS2":53,
}

@contextmanager
def suppress_output():
    """ Contexte pour rendre muets les affichages intempestifs """
    with open(os.devnull, 'w') as devnull:
        with redirect_stdout(devnull), redirect_stderr(devnull):
            try:
                yield
            except Exception:
                raise

def run_with_timeout(func, args, timeout_duration=45.0):
    """ Exécute une fonction réseau dans un thread démon de manière stable """
    local_results = {}
    def worker():
        try:
            local_results['data'] = func(*args)
            local_results['success'] = True
        except Exception as e:
            local_results['success'] = False
            local_results['error'] = e

    thread = threading.Thread(target=worker)
    thread.daemon = True
    thread.start()
    thread.join(timeout=timeout_duration)
    
    if thread.is_alive():
        # On abandonne proprement le thread réseau en tâche de fond (safe car daemon=True)
        raise TimeoutError(f'Timeout after {timeout_duration}s.')
    if not local_results.get('success', False):
        raise local_results.get('error', RuntimeError('Unknown failure.'))
    return local_results['data']

def normalize_formula(formula):
    """ Compte précisément les atomes d'une formule (gère la casse et l'ordre) """
    if not formula or not isinstance(formula, str):
        return Counter()
    
    # Nettoie et standardise la casse de la chaîne reçue (ex: ch3cl -> CH3Cl)
    clean_formula = re.sub(r'([a-zA-Z])', lambda m: m.group(1).upper(), formula)
    clean_formula = clean_formula.replace("CL", "Cl").replace("BR", "Br").replace("SO", "So")
    
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', clean_formula)
    counter = Counter()
    for atom, n in tokens:
        counter[atom] += int(n) if n else 1
    return counter

def get_formula_from_name(name):
    """ Interroge PubChem pour traduire un nom commun en formule brute """
    try:
        results = pcp.get_compounds(name, 'name')
        if results:
            return results[0].molecular_formula
    except Exception:
        pass
    return None

def calculate_doppler_width(mol, nu_GHz):
    """ Extrait la masse molaire et calcule la largeur Doppler associée """
    Mw = None
    if hasattr(mol, 'molar_mass'):
        if isinstance(mol.molar_mass, dict):
            Mw = mol.molar_mass.get('chemical') or mol.molar_mass.get('pubchem') or mol.molar_mass.get('nist')
        else:
            for attr in ('chemical', 'pubchem', 'nist'):
                if hasattr(mol.molar_mass, attr):
                    Mw = getattr(mol.molar_mass, attr)
                    if Mw: break
    
    if Mw is None:
        return None
        
    T = getattr(mol, 'temperature', 298.15)
    m = (Mw / 1000.0) / constants.Avogadro
    nu_hz = nu_GHz * 1e9
    delta = (nu_hz / constants.c) * np.sqrt(8 * constants.k * T * np.log(2) / m)
    return delta * 1e-6

def _init_spectre(mol):
    """ Initialise la structure mol.spectre si elle est absente """
    class SpectreData:
        def __init__(self):
            self.nu_cm = np.array([])
            self.nu_GHz = np.array([])
            self.log_int = np.array([])
            self.coef = np.array([])
            self.hitran = None
            self.databank = 'none'
            self.doppler_width = None
            
    if not hasattr(mol, 'spectre') or mol.spectre is None:
        mol.spectre = SpectreData()

def _store_lines(mol, nu_GHz, log_int_nm2MHz, databank):
    """ Stocke et organise les raies collectées dans l'objet molécule """
    nu_GHz = np.array(nu_GHz, dtype=float)
    log_int = np.array(log_int_nm2MHz, dtype=float)
    nu_cm = nu_GHz * 1e9 / (constants.c * 100)

    imax = np.argmax(log_int) if len(log_int) > 0 else 0

    mol.spectre.nu_cm = nu_cm
    mol.spectre.nu_GHz = nu_GHz
    mol.spectre.log_int = log_int
    mol.spectre.coef = log_int 
    mol.spectre.hitran = nu_GHz[imax] if len(nu_GHz) > 0 else None
    mol.spectre.databank = databank
    mol.spectre.doppler_width = (
        calculate_doppler_width(mol, nu_GHz[imax]) if len(nu_GHz) > 0 else None
    )

# ── Helpers de sélection d'état fondamental (CDMS / JPL) ────────────────────
def _base_formula(s):
    return str(s).split(",")[0].strip().upper()

EXCITED_MARKERS = ("v2=", "v3=", "v1=", "nu2", "nu3", "nu1", "hot", " l=")

def _is_excited(mol_str):
    return any(m in mol_str.strip().lower() for m in EXCITED_MARKERS)

def _is_v0(mol_str):
    return "v=0" in mol_str.lower()

def _select_ground_state(matches):
    v0 = [r for r in matches if _is_v0(str(r["molecule"]))]
    if v0: return v0[0]
    ground = [r for r in matches if not _is_excited(str(r["molecule"]))]
    pool = ground if ground else list(matches)
    return min(pool, key=lambda r: len(str(r["molecule"])))


# ── Priorité 1 : CDMS ────────────────────────────────────────────────────────
def _fetch_cdms(mol, fmin_GHz, fmax_GHz):
    import astropy.units as u
    from astroquery.linelists.cdms import CDMS

    working = mol.formula or mol.name
    query = working.upper().strip()
    species = CDMS.get_species_table()

    hits = [r for r in species if _base_formula(r["molecule"]) == query or _base_formula(r["Name"]) == query]
    if not hits:
        raise LookupError(f'CDMS: molecule {query} not found.')

    selected = _select_ground_state(hits)
    tag = str(int(selected["tag"])).zfill(6)
    name = selected["Name"]

    table = CDMS.query_lines(
        min_frequency = fmin_GHz * u.GHz,
        max_frequency = fmax_GHz * u.GHz,
        min_strength = -6.0,
        molecule = f"{tag} {name}",
    )
    if len(table) == 0:
        raise LookupError('CDMS: no lines returned.')

    freqs_GHz = np.array(table["FREQ"].data, dtype=float) / 1e3
    log_ints = np.array(table["LGINT"].data, dtype=float)

    T = getattr(mol, 'temperature', 300.0)
    if T != 300.0:
        elo = np.array(table["ELO"].data, dtype=float)
        log_ints += elo * (1.0/(0.695 * 300.0) - 1.0/(0.695 * T)) / np.log(10)

    return freqs_GHz, log_ints, f'cdms ({name})'


# ── Priorité 2 : JPL ─────────────────────────────────────────────────────────
def _fetch_jpl(mol, fmin_GHz, fmax_GHz):
    import astropy.units as u
    from astroquery.jplspec import JPLSpec

    working = mol.formula or mol.name
    query = working.upper().strip()
    species = JPLSpec.get_species_table()

    hits = []
    for row in species:
        for col in ("NAME", "name", "MOLECULE", "molecule"):
            if col in row.colnames and _base_formula(row[col]) == query:
                hits.append(row)
                break

    if not hits:
        raise LookupError(f'JPL: molecule {query} not found.')

    selected = _select_ground_state(hits)
    jpl_tag = None
    for col in ("TAG", "tag"):
        if col in selected.colnames:
            jpl_tag = int(selected[col])
            break

    table = JPLSpec.query_lines(
        min_frequency = fmin_GHz * u.GHz,
        max_frequency = fmax_GHz * u.GHz,
        min_strength = -6.0,
        molecule = jpl_tag,
    )
    if len(table) == 0:
        raise LookupError('JPL: no lines returned.')

    freqs_GHz = np.array(table["FREQ"].data, dtype=float) / 1e3
    log_ints = np.array(table["LGINT"].data, dtype=float)

    T = getattr(mol, 'temperature', 300.0)
    if T != 300.0:
        elo = np.array(table["ELO"].data, dtype=float)
        log_ints += elo * (1.0/(0.695 * 300.0) - 1.0/(0.695 * T)) / np.log(10)

    name = selected.get("NAME") or selected.get("name") or working
    return freqs_GHz, log_ints, f'jpl ({name})'


# ── Priorité 3 : HITRAN ───────────────────────────────────────────────────────
def _fetch_hitran(mol, fmin_GHz, fmax_GHz):
    import astropy.units as u
    from astroquery.hitran import Hitran

    working = (mol.formula or mol.name).upper().strip()
    mol_num = HITRAN_MOLECULES.get(working)
    if mol_num is None:
        raise LookupError(f'HITRAN: molecule {working} not found.')

    # Conversion propre et sécurisée vers l'unité attendue par HITRAN (cm^-1)
    fmin_cm = (fmin_GHz * u.GHz).to(u.cm**-1, equivalencies=u.spectral())
    fmax_cm = (fmax_GHz * u.GHz).to(u.cm**-1, equivalencies=u.spectral())

    wmin = min(fmin_cm.value, fmax_cm.value)
    wmax = max(fmin_cm.value, fmax_cm.value)

    table = Hitran.query_lines(
        molecule_number = mol_num,
        isotopologue_number = 1,
        min_frequency = wmin * u.cm**-1,
        max_frequency = wmax * u.cm**-1,
    )
    if len(table) == 0:
        raise LookupError(f'HITRAN: no lines returned for {working} between {wmin:.2f} and {wmax:.2f} cm-1')

    # Reconversion des résultats vers l'unité GHz
    freqs_GHz = (np.array(table["nu"].data, dtype=float) * u.cm**-1).to(u.GHz, equivalencies=u.spectral()).value
    sw = np.array(table["sw"].data, dtype=float)
    elo = np.array(table["elower"].data, dtype=float)

    sw_nm2MHz = sw * C_HITRAN_TO_NM2MHZ
    log_ints = np.log10(np.where(sw_nm2MHz > 0, sw_nm2MHz, 1e-99))
    
    T = getattr(mol, 'temperature', 296.0)
    if T != 296.0:
        log_ints += elo * (1.0/(0.695 * 296.0) - 1.0/(0.695 * T)) / np.log(10)

    return freqs_GHz, log_ints, f'hitran ({working})'


# ── Point d'entrée principal ──────────────────────────────────────────────────
def fill_from_radis(mol):
    """
    Tente de charger le spectre d'une molécule via une cascade CDMS -> JPL -> HITRAN.
    Gère intelligemment l'isomérie de formule brute et l'alignement des unités.
    """
    _init_spectre(mol)

    # 1. Récupération de l'identifiant initial
    working = getattr(mol, 'formula', None) or getattr(mol, 'name', None)
    if not working:
        print("[LINES] Erreur : La molécule n'a ni nom ni formule définie.")
        return "[WARN] [LINES]: No identifier found."
        
    working = str(working).strip()

    # 2. Si c'est un nom vernaculaire (ex: "Acétonitrile"), on obtient la formule par PubChem
    if not re.match(r'^[A-Z]', working):
        pubchem = get_formula_from_name(working)
        if pubchem:
            working = pubchem

    # 3. RECONNAISSANCE ET ALIGNEMENT AUTOMATIQUE PAR FORMULE BRUTE
    # Permet de traduire automatiquement 'C2H3N' en 'CH3CN' ou 'CH4O' en 'CH3OH'
    working_counts = normalize_formula(working)
    if working_counts:
        for known_mol in HITRAN_MOLECULES.keys():
            if normalize_formula(known_mol) == working_counts:
                if working != known_mol:
                    print(f"[DEBUG RADIS] Molécule détectée par isomérie : {working} -> alignée sur {known_mol}")
                working = known_mol
                break

    # Sauvegarde de la formule standardisée sur l'objet molécule
    mol.formula = working

    # 4. SÉCURITÉ ET INTERPRÉTATION DES UNITÉS DE FRÉQUENCE
    fmin_raw = mol.fmin
    fmax_raw = mol.fmax

    if fmin_raw > 1e6:
        # Si la valeur est géante, elle était stockée en Hz -> Passage en GHz
        fmin_GHz = fmin_raw * 1e-9
        fmax_GHz = fmax_raw * 1e-9
    else:
        # Sinon, elle était déjà exprimée en GHz
        fmin_GHz = fmin_raw
        fmax_GHz = fmax_raw

    print(f"[DEBUG RADIS] Requête finale pour {working} : {fmin_GHz:.4f} GHz à {fmax_GHz:.4f} GHz")

    # Cascade des bases de données
    sources = [
        ('CDMS',   _fetch_cdms),
        ('JPL',    _fetch_jpl),
        ('HITRAN', _fetch_hitran),
    ]

    for source_name, fetch_func in sources:
        try:
            nu_GHz, log_int, label = run_with_timeout(
                fetch_func, [mol, fmin_GHz, fmax_GHz], timeout_duration=30.0
            )
            _store_lines(mol, nu_GHz, log_int, label)
            print(f'[LINES] OK : {working} — {len(nu_GHz)} raies via {source_name}')
            return f'[OK] [LINES]: {working} via {source_name} ({len(nu_GHz)} raies)'
        except Exception as e:
            print(f'[LINES] {source_name} failed: {e}')
            continue

    # Échec complet de la cascade
    print(f'[LINES] All sources failed for {working}. Spectre vide.')
    return f'[WARN] [LINES]: No lines found for {working}.'