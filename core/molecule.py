# -*- coding: utf-8 -*-

from core.tools import clean_float

class Molecule:
    
    def __init__(self, name, temperature=298.15, pressure=1, isotope=1, fmin=30, fmax=800, resolution=0.1):
        self.name = name
        self.temperature = temperature
        self.pressure = pressure/1.01325 # bar to atm
        self.isotope = isotope
        self.fmin = fmin
        self.fmax = fmax
        self.resolution = resolution
        
        # Universal molecule identifiers
        self.formula = None
        self.CAS = None
        self.SMILES = None
        self.iupac_name = None
        
        # Identifiers per database
        self.id_hitran = None
        self.id_pubchem = None # PubChem CID token
        self.id_chemical = None
        self.id_nist = None   
        
        # Sub-property objects initialization
        self.molar_mass = MolarMass()   
        self.logp = logP()
        self.boilingtemp = BoilingTemp()
        self.meltingtemp = MeltingTemp()
        self.vaporpressure = VaporPressure()
        self.dipolemoment = DipoleMoment()
        self.vaporizationenthalpy = VaporizationEnthalpy()
        self.fusionenthalpy = FusionEnthalpy()
        self.henryconstant = HenryConstant()
        self.security = Security()
        self.spectre = Spectre()


class MolarMass:
    
    def __init__(self):
        self.pubchem = None
        self.chemical = None
        self.nist = None
        self.mvoc = None

    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 3)} g/mol, "
                f"Chemical = {clean_float(self.chemical, 3)} g/mol, "
                f"NIST = {clean_float(self.nist, 3)} g/mol, "
                f"mVOC = {clean_float(self.mvoc, 3)} g/mol")


class logP:
    
    def __init__(self):
        self.chemical = None
        self.pubchem = None
        self.mvoc = None
       
    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 2)}, "
                f"Chemical = {clean_float(self.chemical, 2)}, "
                f"mVOC = {clean_float(self.mvoc, 2)}")
                

class BoilingTemp:
    
    def __init__(self):
        self.chemical = None
        self.pubchem = None
        self.nist = None
        self.mvoc = None
       
    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 2)} °C, "
                f"Chemical = {clean_float(self.chemical, 2)} °C "
                f"NIST = {clean_float(self.nist, 2)} °C, "
                f"mVOC = {clean_float(self.mvoc, 2)} °C")

class MeltingTemp:
    
    def __init__(self):
        self.chemical = None
        self.pubchem = None
        self.nist = None
        self.mvoc = None
       
    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 2)} °C, "
                f"Chemical = {clean_float(self.chemical, 2)} °C "
                f"NIST = {clean_float(self.nist, 2)} °C, "
                f"mVOC = {clean_float(self.mvoc, 2)} °C")


class VaporPressure:
    
    def __init__(self):
        self.pubchem = None
        self.pubchem_temperature = None
        self.chemical = None
        self.chemical_temperature = None
        self.nist = None
        self.nist_temperature = None
        self.mvoc = None
        self.mvoc_temperature = None
        
        self.temperature = None
       
    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 2)} Pa at {clean_float(self.pubchem_temperature, 2)} K, "
                f"Chemical = {clean_float(self.chemical, 2)} Pa at {clean_float(self.chemical_temperature, 2)} K, "
                f"NIST = {clean_float(self.nist, 2)} Pa at {clean_float(self.nist_temperature, 2)} K, "
                f"mVOC = {clean_float(self.mvoc, 2)} Pa at {clean_float(self.mvoc_temperature, 2)} K")
    

class DipoleMoment:
    
    def __init__(self):
        self.pubchem = None
        self.chemical = None
        self.nist = None
       
    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 2)} D, "
                f"Chemical = {clean_float(self.chemical, 2)} D, "
                f"NIST = {clean_float(self.nist, 2)} D")
    

class VaporizationEnthalpy:
    
    def __init__(self):
        self.pubchem = None
        self.chemical = None
        self.nist = None
        self.mvoc = None
        self.pubchem_temperature = None
        self.nist_temperature = None

    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 2)} kJ/mol at {self.pubchem_temperature} K, "
                f"Chemical = {clean_float(self.chemical, 2)} kJ/mol "
                f"NIST = {clean_float(self.nist, 2)} kJ/mol at {self.nist_temperature} K")
    
    
class FusionEnthalpy:
    
    def __init__(self):
        self.pubchem = None
        self.chemical = None
        self.nist = None
        self.pubchem_temperature = None
        self.nist_temperature = None
       
    def __repr__(self):
        return (f"PubChem = {clean_float(self.pubchem, 2)} kJ/mol at {self.pubchem_temperature} K, "
                f"Chemical = {clean_float(self.chemical, 2)} kJ/mol "
                f"NIST = {clean_float(self.nist, 2)} kJ/mol at {self.nist_temperature} K")
                

class HenryConstant:
    
    def __init__(self):
        self.sander = None
        self.sander_dT = None
        self.sander_temperature = None
        self.mvoc = None
        self.mvoc_temperature = None
        self.nist = None
        self.nist_temperature = None
        
    def __repr__(self):
        return (f"HenryConstant(Sander = {clean_float(self.sander, 2)} /atm at {self.sander_temperature} K, "
                f"mVOC = {clean_float(self.mvoc, 2)} /atm at {self.mvoc_temperature} K, "
                f"NIST = {clean_float(self.nist, 2)} /atm at {self.nist_temperature} K")


class Security:
    
    def __init__(self):
        self.signal_word = []
        self.H_codes = {}
        self.P_codes = []
    
    def __repr__(self):
        return (f"Signal = {self.signal_word}, "
                f"H = {list(self.H_codes.keys())}, "
                f"P = {self.P_codes}")


class Spectre:
    
    def __init__(self):
        self.hitran = None # Peak frequency (GHz)
        self.doppler_width = None # Doppler width FWHM (GHz)
        self.databank = None # Data source origin (hitran, exomol, etc.)

    def __repr__(self):
        return (f"Peak frequency = {clean_float(self.hitran, 2)} GHz, "
                f"Doppler width = {clean_float(self.doppler_width, 3)} MHz")