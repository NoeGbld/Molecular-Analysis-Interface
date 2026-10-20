# -*- coding: utf-8 -*-

from .parsers.pubchem_parser import fill_from_pubchem
from .parsers.astroquery_parser import fill_from_astroquery
from .parsers.chemical_parser import fill_from_chemical
from .parsers.nist_parser import fill_from_nist
from .parsers.mvoc_parser_python import fill_from_mvoc
from .parsers.sander_parser import fill_from_sander

# Sequentially executes all database parsers to fully populate a Molecule object dataset
def solve(mol):
    fill_from_pubchem(mol)
    fill_from_astroquery(mol)
    fill_from_chemical(mol)
    fill_from_nist(mol)
    fill_from_mvoc(mol)
    fill_from_sander(mol)

    return mol