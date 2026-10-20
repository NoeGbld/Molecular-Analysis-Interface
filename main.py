# -*- coding: utf-8 -*-

from core.molecule import Molecule
from core.solver import solve
from core.plot_python import plot_absorption_spectrum
from core.tools import install_playwright_browsers

# Check dependencies before execution
install_playwright_browsers()

# Initialize molecule configuration with experimental targets and spectral ranges
mol = Molecule(
    'acetonitrile', # Target chemical name, formula or PubChem CID
    temperature=298.15, # Ambient temperature condition in Kelvin (K)
    pressure=1, # Ambient atmospheric pressure condition in atmospheres (bar)
    isotope=1, # Selected primary target isotope index
    fmin=30, # Minimum spectrum frequency boundary in GHz
    fmax=800, # Maximum spectrum frequency boundary in GHz
)

# Execute query orchestration loop
solve(mol)

# Output extracted datasets onto console
print()
print("=== Formula ==="); print(mol.formula)
print()
print("=== IUPAC Name ==="); print(mol.iupac_name)
print()
print("=== SMILES ==="); print(mol.SMILES)
print()
print("=== Molar Mass ==="); print(mol.molar_mass)
print()
print("=== logP ==="); print(mol.logp)
print()
print("=== Boiling Point ==="); print(mol.boilingtemp)
print()
print("=== Melting Point ==="); print(mol.meltingtemp)
print()
print("=== Vapor Pressure ==="); print(mol.vaporpressure)
print()
print("=== Dipole Moment ==="); print(mol.dipolemoment)
print()
print("=== Vaporization Enthalpy ==="); print(mol.vaporizationenthalpy)
print()
print("=== Fusion Enthalpy ==="); print(mol.fusionenthalpy)
print()
print("=== Henry Constant ==="); print(mol.henryconstant)
print()
print("=== Security ==="); print(mol.security)
print()
print("=== Spectrum ==="); print(mol.spectre)

# Render and present absorption plot figure
plot_absorption_spectrum(mol)