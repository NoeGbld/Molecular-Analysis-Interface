# -*- coding: utf-8 -*-

import pint
import re
import numpy as np
import subprocess
import sys

def install_playwright_browsers():
    # Forces a system installation pass on missing Chromium headless browser instances
    # This ensures that Playwright scrapers (like mVOC) have the required binaries
    try:
        subprocess.check_call([sys.executable, "-m", "playwright", "install", "chromium"])
    except Exception as e:
        print(f"[System] Playwright installation error: {e}")

# Initialize unit registry object
ureg = pint.UnitRegistry()

# Definition of compatible target units conversion maps
UNIT_TARGETS = {
    'degC': ['degC', 'degF', 'kelvin', 'K'],
    'pascal': ['Pa', 'pascal', 'kPa', 'mmHg', 'mm Hg', 'atm', 'bar'],
    'debye': ['debye', 'D', 'C.m', 'C·m'],
    'joule': ['J', 'kJ', 'cal', 'kcal'],
    'kilojoule / mole': ['kJ/mol', 'J/mol', 'kcal/mol', 'cal/mol'],
}

# Regular expression pattern to capture numbers paired with standard physical units
PATTERN = r'(-?\d+\.?\d*(?:e-?\d+)?)\s*(' \
          r'mol/\(m\^3\*Pa\)|' \
          r'atm\*m\^3/mol|' \
          r'atm\^-1|' \
          r'kJ/mol|J/mol|kcal/mol|cal/mol|' \
          r'kJ|kcal|cal|J|' \
          r'°C|°F|' \
          r'kPa|hPa|mPa|Pa|' \
          r'mm\s?Hg|mmHg|atm|bar|psi|torr|' \
          r'debye|C\.m|C·m|D|K)'


# Normalizes scientific notations, minus signs, and superscript exponents in a string
def preprocess(value_str):
    if value_str is None:
        return None
    
    # Map unicode superscript exponents to standard ASCII digits
    superscripts = {
        '⁰': '0', '¹': '1', '²': '2', '³': '3', '⁴': '4',
        '⁵': '5', '⁶': '6', '⁷': '7', '⁸': '8', '⁹': '9',
        '⁻': '-'
    }
    for sup, normal in superscripts.items():
        value_str = value_str.replace(sup, normal)
    
    value_str = value_str.replace('−', '-')
    
    # Normalize scientific multiplication sign structures
    value_str = re.sub(r'(\d+(?:\.\d+)?)\s*×\s*10\s*([+-]?\d+)', r'\1e\2', value_str)
    value_str = re.sub(r'[xX]10\^?(-?\d+)', r'e\1', value_str)
    
    return value_str


# Extracts a single value and converts its clean unit to standard predefined targets
def convert_all(value_str):
    if len(value_str) > 80:
        return None

    value_str = preprocess(value_str)
    matches = re.findall(PATTERN, value_str, re.IGNORECASE)

    if not matches:
        return None

    found = {}

    for value, unit in matches:
        unit_clean = unit.strip()
        
        # Replace notation variants to ensure Pint compatibility
        unit_clean = unit_clean.replace('mm Hg', 'mmHg')
        unit_clean = unit_clean.replace('m3', 'm^3')
        unit_clean = unit_clean.replace('°C', 'degC')
        unit_clean = unit_clean.replace('°F', 'degF')
        unit_clean = unit_clean.replace('C.m', 'coulomb * meter')
        unit_clean = unit_clean.replace('C·m', 'coulomb * meter')

        try:
            q = ureg.Quantity(float(value), unit_clean)

            # Match quantity with its compatible dimension category group
            for target, units in UNIT_TARGETS.items():
                if target not in found:
                    try:
                        found[target] = q.to(target).magnitude
                        break
                    except Exception:
                        continue
        except Exception:
            continue

    return found if found else None


# Extracts and parses multiple property-temperature pairs into a list of dictionaries
def convert_all_multiple(value_str):
    if len(value_str) > 200:
        return []
    
    value_str = preprocess(value_str)
    matches = re.findall(PATTERN, value_str, re.IGNORECASE)
    
    if not matches:
        return []
    
    results = []
    current = {}
    
    for value, unit in matches:
        unit_clean = unit.strip()
        unit_clean = unit_clean.replace('mm Hg', 'mmHg')
        unit_clean = unit_clean.replace('°C', 'degC')
        unit_clean = unit_clean.replace('°F', 'degF')
        unit_clean = unit_clean.replace('C.m', 'coulomb * meter')
        unit_clean = unit_clean.replace('C·m', 'coulomb * meter')
        
        try:
            q = ureg.Quantity(float(value), unit_clean)
            
            for target, units in UNIT_TARGETS.items():
                if target not in current:
                    try:
                        current[target] = q.to(target).magnitude
                        break
                    except Exception:
                        continue
                else:
                    # If target key conflict occurs, push previous block into list
                    if current:
                        results.append(current)
                        current = {}
                    try:
                        current[target] = q.to(target).magnitude
                        break
                    except Exception:
                        continue
        except Exception:
            continue
    
    if current:
        results.append(current)
    
    return results


# Converts a standalone numeric value from an explicit source unit to a target unit
def convert_value(value, from_unit, to_unit):
    try:
        q = ureg.Quantity(value, from_unit)
        return q.to(to_unit).magnitude
    except Exception:
        return None


# Extrapolates the reference Henry's law constant to a specific target temperature
def henry_at_T(H, dT, T):
    T_ref = 298.15
    return H * np.exp(dT * (1 / T - 1 / T_ref))


# Formats a floating number using scientific notation or standard rounding thresholds
def clean_float(value, decimals):
    if value is None:
        return None
    val = float(value)
    
    # Use standard scientific notation fallback for very small numbers
    if 0 < abs(val) < 1e-2:
        return float(f"{val:.{decimals}e}")
    return round(val, decimals)