#!/usr/bin/env python

from .tad_signature import compute_tad_signature, to_log_odds

from .tadmap_utils import read_Ensembl_v102_refdata, read_TADMap_from_file_or_url, retrieve_TADMap_by_species, standardize_adata_gene_names

from .tadmap_base_config import set_loglevel

__all__ = ['compute_tad_signature',
           'to_log_odds',
           'read_Ensembl_v102_refdata',
           'read_TADMap_from_file_or_url',
           'retrieve_TADMap_by_species',
           'standardize_adata_gene_names',
           'set_loglevel']


import pkgutil, pathlib, importlib

# from pkgutil import iter_modules
# from pathlib import Path
# from importlib import import_module

# https://julienharbulot.com/python-dynamical-import.html
# iterate through the modules in the current package
#
package_dir = pathlib.Path(__file__).resolve().parent
for (_, module_name, _) in pkgutil.iter_modules([package_dir]):
    if 'data' in module_name:
        module = importlib.import_module(f"{__name__}.{module_name}")
