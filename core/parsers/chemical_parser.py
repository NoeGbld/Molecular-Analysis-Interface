# -*- coding: utf-8 -*-

from thermo import Chemical
from core.tools import convert_value

# Queries the thermo package to extract and convert native chemical data fields
def fill_from_chemical(mol):
    if not mol.CAS:
        # Return warning if the CAS registry number is unavailable
        return "[WARN] Chemical: CAS unavailable"

    try:
        results = Chemical(mol.CAS, T=mol.temperature)
        if results is None:
            return "[WARN] Chemical: Data not found for this CAS"

        # Map unique identifier
        mol.id_chemical = results.CAS

        # Extract and compute chemical properties from thermo package datasets
        mol.molar_mass.chemical = results.MW if results.MW else None    
        mol.logp.chemical = results.logP if results.logP else None
        mol.boilingtemp.chemical = convert_value(results.Tb, 'K', 'degC') if results.Tb else None
        mol.meltingtemp.chemical = convert_value(results.Tm, 'K', 'degC') if results.Tm else None
        mol.vaporpressure.chemical = results.Psat if results.Psat else None
        mol.vaporizationenthalpy.chemical = results.Hvap * results.MW / 1e6 if results.Hvap else None 
        mol.fusionenthalpy.chemical = results.Hfus * results.MW / 1e6 if results.Hfus else None
        mol.dipolemoment.chemical = results.dipole if results.dipole else None

        # Standardize temperature metadata tags for experimental property points
        mol.vaporpressure.chemical_temperature = mol.temperature if results.Psat else None
        mol.fusionenthalpy.chemical_temperature = mol.temperature if results.Hfus else None
        mol.vaporizationenthalpy.chemical_temperature = mol.temperature if results.Hvap else None
    
        print("[Chemical]: OK")
        return "[OK] Chemical data loaded"
    except Exception as e:
        raise ValueError(f"[Chemical] Parsing error: {str(e)}")