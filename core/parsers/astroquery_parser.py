# -*- coding: utf-8 -*-

import re
import os
import sys
import ctypes
import threading
import time
from collections import Counter
from contextlib import contextmanager, redirect_stdout, redirect_stderr
import numpy as np
from scipy import constants
import pubchempy as pcp

# Molecule databases
radis_molecules = [
    'OCS', 'AlCl', 'PF3', 'H2CS', 'C2H6', 'CH', 'CH3CN', 'C2H4', 'BeH',
    'AlH', 'HI', 'CS2', 'CS', 'TiO', 'C2', 'SF6', 'ScH', 'NaCl', 'HOBr',
    'SiH4', 'COCl2', 'HBr', 'KF', 'H2CO', 'HO2', 'NH', 'CaCl', 'AlF',
    'YO', 'H2S', 'ScO', 'MgF', 'CH3', 'VO', 'HBO', 'LiH_p', 'ClONO2',
    'SiN', 'PO', 'NO+', 'CN', 'C2H2', 'H3+', 'CaH', 'NaF', 'MgH',
    'PH', 'NaOH', 'SH', 'NH3', 'ZrO', 'LiH', 'S2', 'N2O', 'CH3I',
    'NO2', 'H2', 'ClO', 'HCOOH', 'CaF', 'CaOH', 'CrH', 'CH3Cl',
    'FeH', 'SO3', 'COF2', 'O3', 'ClNO2', 'AsH3', 'HC3N', 'MgO',
    'NF3', 'SO', 'HF', 'CO', 'COFCl', 'LiOH', 'H2O2', 'HONO', 'PS',
    'C4H2', 'LiF', 'C3', 'SiS', 'SiO2', 'NiH', 'HCN', 'CH3Br', 'CF4',
    'SiH2', 'LiCl', 'PN', 'CH3F', 'CO2', 'CP', 'GeH4', 'PH3', 'CaO',
    'TiH', 'O2', 'SO2', 'N2', 'HOCl', 'NS', 'NaH', 'H2O', 'KOH',
    'NaO', 'OH', 'SiO', 'C2N2', 'NO', 'HNO3', 'KCl', 'CH3OH', 'SiH',
    'LaO', 'AlO', 'HCl', 'CH4'
]

HITRAN_MOLECULES = {
    "H2O":1,   "CO2":2,    "O3":3,    "N2O":4,   "CO":5,    "CH4":6,
    "O2":7,    "NO":8,     "SO2":9,   "NO2":10,  "NH3":11,  "HNO3":12,
    "OH":13,   "HF":14,    "HCl":15,  "HBr":16,  "HI":17,   "ClO":18,
    "OCS":19,  "H2CO":20,  "HOCl":21, "N2":22,   "HCN":23,  "CH3Cl":24,
    "H2O2":25, "C2H2":26,  "C2H6":27, "PH3":28,  "COF2":29, "SF6":30,
    "H2S":31,  "HCOOH":32, "HO2":33,  "O":34,    "ClONO2":35,"NO+":36,
    "HOBr":37, "C2H4":38,  "CH3OH":39,"CH3Br":40,"CH3CN":41,"CF4":42,
    "C4H2":43, "HC3N":44,  "H2":45,   "CS":46,   "SO3":47,  "C2N2":48,
    "COCl2":49,"SO":50,    "C3H4":51, "CH3":52,  "CS2":53,  "CH3I":54,
    "NF3":55,  "HC5N":61
}

# Context manager
@contextmanager
def suppress_output():
    with open(os.devnull, 'w') as devnull:
        with redirect_stdout(devnull), redirect_stderr(devnull):
            try:
                yield
            except (TimeoutError, Exception):
                raise

with suppress_output():
    try:
        from hapi import *
    except ImportError:
        pass

# Common utilities
def normalize_formula(formula):
    if not formula or not isinstance(formula, str):
        return Counter()
    formula = formula[0].upper() + formula[1:]
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    counter = Counter()
    for atom, n in tokens:
        counter[atom] += int(n) if n else 1
    return counter

def get_formula_from_name(name):
    try:
        results = pcp.get_compounds(name, 'name')
        if results:
            return results[0].molecular_formula
    except Exception:
        pass
    return None

def calculate_doppler_width(mol, nu_GHz):
    Mw = None
    if hasattr(mol, 'molar_mass'):
        if isinstance(mol.molar_mass, dict):
            Mw = (mol.molar_mass.get('chemical')
                  or mol.molar_mass.get('pubchem')
                  or mol.molar_mass.get('nist'))
        else:
            for attr in ('chemical', 'pubchem', 'nist'):
                if hasattr(mol.molar_mass, attr):
                    Mw = getattr(mol.molar_mass, attr)
                    if Mw:
                        break
    if Mw is None:
        return None
    m = (Mw / 1000) / constants.Avogadro
    nu_hz = nu_GHz * 1e9
    delta_hz = (nu_hz / constants.c) * np.sqrt(
        8 * constants.k * mol.temperature * np.log(2) / m)
    return delta_hz * 1e-6

# SpectreData helpers
def _init_spectre(mol):
    if not hasattr(mol, 'spectre') or mol.spectre is None:
        class SpectreData:
            nu_cm = []
            nu_GHz = []
            coef = []
            coef_unit = 'linear'
            hitran = None
            databank = 'none'
            doppler_width = None
            type = 'discrete'
        mol.spectre = SpectreData()

def store_spectre(mol, nu_cm, coef, databank,
                  spec_type='continuous', coef_unit='linear'):
    _init_spectre(mol)
    nu_cm = np.array(nu_cm,  dtype=float)
    coef = np.array(coef,   dtype=float)
    nu_GHz = nu_cm * constants.c * 100 * 1e-9

    if len(coef) > 0:
        if spec_type == 'discrete':
            imax = int(np.argmax(coef))
        else:
            valid = coef > 0
            if valid.any():
                imax = int(np.where(valid)[0][np.argmax(coef[valid])])
            else:
                imax = int(np.argmax(coef))
    else:
        imax = 0

    mol.spectre.nu_cm = nu_cm
    mol.spectre.nu_GHz = nu_GHz
    mol.spectre.coef = coef
    mol.spectre.coef_unit = coef_unit
    mol.spectre.hitran = nu_GHz[imax] if len(nu_GHz) > 0 else None
    mol.spectre.databank = databank
    mol.spectre.type = spec_type
    mol.spectre.doppler_width = (
        calculate_doppler_width(mol, nu_GHz[imax])
        if len(nu_GHz) > 0 else None)

def _clear_spectre(mol):
    _init_spectre(mol)
    mol.spectre.nu_cm = []
    mol.spectre.nu_GHz = []
    mol.spectre.coef = []
    mol.spectre.coef_unit = 'linear'
    mol.spectre.hitran = None
    mol.spectre.databank = 'none'
    mol.spectre.type = 'discrete'
    mol.spectre.doppler_width = None

# Thread-based timeout runner
def run_with_timeout(func, args, timeout_duration=45.0):
    local_results = {}

    def worker():
        try:
            local_results['data'] = func(*args)
            local_results['success'] = True
        except Exception as thread_error:
            local_results['success'] = False
            local_results['error'] = thread_error

    thread = threading.Thread(target=worker)
    thread.daemon = True
    thread.start()
    thread.join(timeout=timeout_duration)

    if thread.is_alive():
        thread_id = ctypes.c_long(thread.ident)
        ctypes.pythonapi.PyThreadState_SetAsyncExc(
            thread_id, ctypes.py_object(TimeoutError))
        time.sleep(0.2)
        raise TimeoutError(
            f'Execution exceeded time limit of {timeout_duration}s.')

    if not local_results.get('success', False):
        raise local_results.get('error', RuntimeError(
            'Unknown inner processing failure.'))

    return local_results.get('data')

# Pipeline 1: RADIS
def fill_from_radis(mol, mode='continuous'):
    try:
        from radis import calc_spectrum
    except ImportError:
        raise RuntimeError('[RADIS] radis is not installed')

    _init_spectre(mol)
    mol.id_hitran = None
    mol.formula_hitran = None
    mol.hitran_name = None

    working_formula = mol.formula if mol.formula is not None else mol.name
    if not working_formula or not re.match(r'^[A-Za-z]', str(working_formula)):
        pubchem_formula = get_formula_from_name(mol.name)
        working_formula = pubchem_formula if pubchem_formula else mol.name

    formula_mol = normalize_formula(working_formula)

    radis_name = None
    if formula_mol:
        for known_mol in radis_molecules:
            if normalize_formula(known_mol) == formula_mol:
                radis_name = known_mol
                break
    if radis_name is None:
        radis_name = working_formula

    else:
        try:
            for mol_id, iso_info in ISO.items():
                if mol_id[1] != 1:
                    continue
                if normalize_formula(iso_info[4]) == formula_mol:
                    mol.id_hitran = mol_id[0]
                    mol.formula_hitran = iso_info[4]
                    mol.hitran_name = iso_info[4]
                    break
        except NameError:
            pass

    fmin_cm = (mol.fmin*1e9)/(constants.c*100)
    fmax_cm = (mol.fmax*1e9)/(constants.c*100)

    if mol.resolution > 1e3:
        resolution_cm = mol.resolution /(constants.c*100)
    else:
        resolution_cm = (mol.resolution*1e9)/(constants.c*100)

    if mode == 'discrete':
        if mol.formula_hitran is not None and 'hapi' in sys.modules:
            try:
                def call_hapi_discrete():
                    with suppress_output():
                        db_begin('data')
                        fetch(mol.formula_hitran, mol.id_hitran, 1,
                              fmin_cm, fmax_cm)
                        nu_arr, sw_arr = getColumns(
                            mol.formula_hitran, ['nu', 'sw'])
                    return np.array(nu_arr), np.array(sw_arr)

                nu_cm2, sw = call_hapi_discrete()

                mask = (nu_cm2 >= fmin_cm) & (nu_cm2 <= fmax_cm)
                nu_cm2 = nu_cm2[mask]
                sw = sw[mask]

                if len(sw) == 0:
                    raise LookupError('No lines in range after masking.')

                log_ints = np.log10(np.where(sw > 0, sw, 1e-99))

                store_spectre(mol, nu_cm2, log_ints,
                              'hitran_hapi',
                              spec_type='discrete',
                              coef_unit='log10')
                print(f'[RADIS discrete]: {mol.formula_hitran} via HAPI')
                return f'[OK] [RADIS discrete]: {mol.formula_hitran} via HAPI'

            except Exception as e:
                print(f'[RADIS discrete] HAPI failed: {e}')

        _clear_spectre(mol)
        print(f'[RADIS discrete]: No discrete data for {working_formula}.')
        return f'[WARN] [RADIS discrete]: No discrete data for {working_formula}.'

    if mol.formula_hitran is not None and 'hapi' in sys.modules:
        try:
            def call_hapi_continuous():
                with suppress_output():
                    db_begin('data')
                    fetch(mol.formula_hitran, mol.id_hitran, 1, fmin_cm, fmax_cm)
                    return absorptionCoefficient_Voigt(
                        SourceTables = mol.formula_hitran,
                        Environment = {'T': mol.temperature, 'p': mol.pressure},
                        WavenumberStep = resolution_cm)

            nu_cm, coef = call_hapi_continuous()
            store_spectre(mol, nu_cm, coef, 'hitran_hapi',
                          spec_type='continuous', coef_unit='linear')
            print(f'[RADIS continuous]: {mol.formula_hitran} via HAPI')
            return f'[OK] [RADIS continuous]: {mol.formula_hitran} via HAPI'
        except Exception as e:
            print(f'[RADIS continuous] HAPI failed: {e}')

    if radis_name is not None:
        def run_fallback_databases(db_list, mol_name):
            for databank in db_list:
                try:
                    with suppress_output():
                        s = calc_spectrum(
                            fmin_cm, fmax_cm,
                            molecule = mol_name,
                            isotope = mol.isotope,
                            pressure = mol.pressure,
                            Tgas = mol.temperature,
                            mole_fraction = 1,
                            path_length = 1,
                            databank = databank,
                            verbose = False,
                            wstep = resolution_cm)
                    return s.get('crosssection', wunit='cm-1'), databank
                except Exception:
                    continue
            raise RuntimeError('All authorized databases exhausted.')

        try:
            target_dbs = ['hitran', 'hitemp', 'geisa', 'exomol']
            result_data = run_with_timeout(
                run_fallback_databases,
                [target_dbs, radis_name],
                timeout_duration=20.0)
            (nu_cm, crosssection_cm2), successful_db = result_data

            store_spectre(mol, nu_cm, crosssection_cm2,
                          successful_db,
                          spec_type='continuous',
                          coef_unit='linear')
            print(f'[RADIS continuous]: {radis_name} via {successful_db.upper()}')
            return (f'[OK] [RADIS continuous]: '
                    f'{radis_name} via {successful_db.upper()}')

        except Exception as e:
            print(f'[RADIS continuous] Targeted Fallback DBs failed: {e}')

    _clear_spectre(mol)
    print(f'[RADIS continuous]: Global failure for {working_formula}.')
    return f'[WARN] [RADIS continuous]: Global failure for {working_formula}.'

# Pipeline 2: Astroquery
def _base_formula(s):
    return str(s).split(",")[0].strip().upper()

def _is_excited(s):
    return any(m in str(s).lower()
               for m in ("v2=", "v3=", "v1=", "nu2", "nu3", "nu1",
                         "hot", " l="))

def _select_ground_state(hits):
    v0 = [r for r in hits if "v=0" in str(r["molecule"]).lower()]
    if v0:
        return v0[0]
    ground = [r for r in hits if not _is_excited(r["molecule"])]
    pool = ground if ground else list(hits)
    return min(pool, key=lambda r: len(str(r["molecule"])))

# Fallback: direct fetch of .cat file from CDMS when query_lines failed
def _fetch_cdms_raw(tag, fmin_GHz, fmax_GHz):
    import urllib.request
    url = (f"https://cdms.astro.uni-koeln.de/classic/entries/"
           f"c{int(tag):06d}.cat")
    with urllib.request.urlopen(url, timeout=20) as resp:
        content = resp.read().decode('utf-8', errors='replace')
    freqs, log_ints, elos = [], [], []
    for line in content.splitlines():
        if len(line) < 41:
            continue
        try:
            freq_ghz = float(line[0:13])/1e3
            lgint = float(line[21:29])
            elo = float(line[31:41])
        except ValueError:
            continue
        if fmin_GHz <= freq_ghz <= fmax_GHz:
            freqs.append(freq_ghz)
            log_ints.append(lgint)
            elos.append(elo)
    return np.array(freqs), np.array(log_ints), np.array(elos)

def _fetch_cdms(mol, fmin_GHz, fmax_GHz):
    import astropy.units as u
    from astroquery.linelists.cdms import CDMS

    query = mol.formula.upper().strip()
    q_counts = normalize_formula(query)

    species = CDMS.get_species_table()
    hits = [r for r in species
               if _base_formula(r["molecule"]) == query
               or _base_formula(r["Name"]) == query
               or (q_counts and normalize_formula(
                       re.sub(r'^[a-z]-', '',
                              str(r["Name"]).split(",")[0].strip())
                   ) == q_counts)]
    if not hits:
        raise LookupError('Not found.')

    sel   = _select_ground_state(hits)
    T_ref = mol.temperature

    try:
        table = CDMS.query_lines(
            fmin_GHz * u.GHz, fmax_GHz * u.GHz,
            min_strength = -6.0,
            molecule = f"{str(int(sel['tag'])).zfill(6)} {sel['Name']}")
        if len(table) == 0:
            raise LookupError('No lines.')
        freqs    = np.array(table["FREQ"].data, dtype=float) / 1e3
        log_ints = np.array(table["LGINT"].data, dtype=float)
        if getattr(mol, 'temperature', T_ref) != T_ref:
            log_ints += (np.array(table["ELO"].data, dtype=float)
                         * (1.0 / (0.695 * T_ref)
                            - 1.0 / (0.695 * mol.temperature))
                         / np.log(10))
    except ValueError:
        freqs, log_ints, elos = _fetch_cdms_raw(sel['tag'], fmin_GHz, fmax_GHz)
        if len(freqs) == 0:
            raise LookupError('No lines (raw).')
        if getattr(mol, 'temperature', T_ref) != T_ref:
            log_ints = log_ints + (elos
                         * (1.0 / (0.695 * T_ref)
                            - 1.0 / (0.695 * mol.temperature))
                         / np.log(10))

    return freqs, log_ints, f"cdms ({sel['Name']})"

def _fetch_jpl(mol, fmin_GHz, fmax_GHz):
    import astropy.units as u
    from astroquery.jplspec import JPLSpec

    query = mol.formula.upper().strip()
    species = JPLSpec.get_species_table()
    hits = [r for r in species
            if any(_base_formula(r[c]) == query
                      for c in r.colnames
                      if c.lower() in ('name', 'molecule'))]
    if not hits:
        raise LookupError('Not found.')

    sel = _select_ground_state(hits)
    tag = int(sel["TAG"]) if "TAG" in sel.colnames else int(sel["tag"])
    table = JPLSpec.query_lines(
        fmin_GHz * u.GHz, fmax_GHz * u.GHz,
        min_strength = -6.0,
        molecule = tag)

    if len(table) == 0:
        raise LookupError('No lines.')

    freqs = np.array(table["FREQ"].data, dtype=float) / 1e3
    log_ints = np.array(table["LGINT"].data, dtype=float)

    T_ref = 300.0
    if getattr(mol, 'temperature', T_ref) != T_ref:
        log_ints += (np.array(table["ELO"].data, dtype=float)
                     * (1.0 / (0.695 * T_ref)
                        - 1.0 / (0.695 * mol.temperature))
                     / np.log(10))

    return freqs, log_ints, f"jpl ({sel.get('NAME') or query})"

def _fetch_hitran(mol, fmin_GHz, fmax_GHz):
    import astropy.units as u
    from astroquery.hitran import Hitran

    query = mol.formula.strip()
    mol_key = next((k for k in HITRAN_MOLECULES
                    if k.upper() == query.upper()), None)
    if mol_key is None:
        raise LookupError('Not found.')
    mol_num = HITRAN_MOLECULES[mol_key]

    fmin_cm = (fmin_GHz * u.GHz).to(u.cm**-1, equivalencies=u.spectral())
    fmax_cm = (fmax_GHz * u.GHz).to(u.cm**-1, equivalencies=u.spectral())

    table = Hitran.query_lines(
        molecule_number = mol_num,
        isotopologue_number = 1,
        min_frequency = min(fmin_cm.value, fmax_cm.value) * u.cm**-1,
        max_frequency = max(fmin_cm.value, fmax_cm.value) * u.cm**-1)

    if len(table) == 0:
        raise LookupError('No lines.')

    freqs = (np.array(table["nu"].data, dtype=float) * u.cm**-1
             ).to(u.GHz, equivalencies=u.spectral()).value

    sw = np.array(table["sw"].data, dtype=float)
    log_ints = np.log10(np.where(sw > 0, sw, 1e-99))

    T_ref = 296.0
    if getattr(mol, 'temperature', T_ref) != T_ref:
        log_ints += (np.array(table["elower"].data, dtype=float)
                     * (1.0 / (0.695 * T_ref)
                        - 1.0 / (0.695 * mol.temperature))
                     / np.log(10))

    return freqs, log_ints, f"hitran ({mol_key})"

def fill_from_astroquery(mol):
    _init_spectre(mol)
    working = getattr(mol, 'formula', None) or getattr(mol, 'name', None)
    if not working:
        return "[WARN] [Astroquery]: No identifier."
    working = str(working).strip()

    if not re.match(r'^[A-Za-z]', working):
        pb = get_formula_from_name(working)
        if pb:
            working = pb

    w_counts = normalize_formula(working)
    if w_counts:
        mol_key = next(
            (k for k in HITRAN_MOLECULES
             if normalize_formula(k) == w_counts),
            None)
        if mol_key:
            working = mol_key

    mol.formula = working

    fmin_GHz = mol.fmin
    fmax_GHz = mol.fmax

    for src_name, fetch_f in [('CDMS',   _fetch_cdms),
                               ('JPL',    _fetch_jpl),
                               ('HITRAN', _fetch_hitran)]:
        try:
            freqs_GHz, log_ints, lbl = run_with_timeout(
                fetch_f,
                [mol, fmin_GHz, fmax_GHz],
                timeout_duration=25.0)

            nu_cm = freqs_GHz * 1e9 / (constants.c * 100)
            store_spectre(mol, nu_cm, log_ints, lbl,
                          spec_type='discrete',
                          coef_unit='log10')
            print(f'[Astroquery]: {working} via {src_name}')
            return f'[OK] [Astroquery]: {working} via {src_name}'

        except Exception:
            continue

    _clear_spectre(mol)
    print(f'[Astroquery]: No lines found for {working}.')
    return f'[WARN] [Astroquery]: No lines found for {working}.'