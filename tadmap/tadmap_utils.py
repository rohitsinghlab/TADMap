#!/usr/bin/env python

###################################################################
## Primary Author:  Rohit Singh rsingh@alum.mit.edu
## License: MIT
## Repository:  http://github.io/rs239/tadmap
###################################################################

import pandas as pd
import numpy as np
import scipy, sklearn, os, sys, string, fileinput, glob, re, math, itertools, functools, copy, multiprocessing, traceback, tarfile, gzip, csv
import scipy.stats, sklearn.decomposition, sklearn.preprocessing, sklearn.covariance
from scipy.stats import describe
from scipy import sparse
import os.path
import scipy.sparse
from scipy.sparse import csr_matrix, csc_matrix
from sklearn.preprocessing import normalize
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

from . import  tadmap_base_config
dbg_print = tadmap_base_config.dbg_print


def read_Ensembl_v102_refdata(sp):
    """ Returns a dataframe consisting of gene symbology, location etc. for Ensembl v102

    :param sp: string, one of `hs` (human) or `mm` (mouse)

    :returns: A dataframe with gene name, chromosome, strad direction, TSS etc., 
"""

    assert sp in ['hs','mm']

    if sp=='mm':
        df = pd.read_csv("http://cb.csail.mit.edu/cb/tadmap/TADMap_mm10_ens102_refdata.tsv", delimiter='\t')
        df = df[~df['MGI symbol'].isnull()].reset_index(drop=True)
        df = df[df['Chromosome/scaffold name'].apply(lambda s: len(s)<=2 and s!="MT")].reset_index(drop=True)
        df["gene1"] = np.where(~df["Gene name"].isnull(), df["Gene name"], df["Gene stable ID"])
        df1 = pd.DataFrame({ "gene": df['gene1'] + "|Ensembl:" + df['Gene stable ID'] + "|MGI:" + df['MGI symbol'],
                             "chr": df['Chromosome/scaffold name'].apply(lambda s: "chr"+s if len(s)<=2 and s!="MT" else s.lower()),
                             "strand": df["Strand"].apply(lambda a: "+" if int(a)>0 else "-"),
                             "txstart": df['Gene start (bp)'], "txend": df['Gene end (bp)']})
        return df1

    if sp=='hs':
        df = pd.read_csv("http://cb.csail.mit.edu/cb/tadmap/TADMap_mm10_ens102_refdata.tsv", delimiter='\t')
        df = df[~df['HGNC symbol'].isnull()].reset_index(drop=True)
        df = df[df['Chromosome/scaffold name'].apply(lambda s: len(s)<=2 and s!="MT")].reset_index(drop=True)
        df["gene1"] = np.where(~df["Gene name"].isnull(), df["Gene name"], df["Gene stable ID"])
        df1 = pd.DataFrame({ "gene": df['Gene name'] + "|Ensembl:" + df['Gene stable ID'] + "|HGNC:" + df['HGNC symbol'],
                             "chr": df['Chromosome/scaffold name'].apply(lambda s: "chr"+s if len(s)<=2 and s!="MT" else s.lower()),
                             "strand": df["Strand"].apply(lambda a: "+" if int(a)>0 else "-"),
                             "txstart": df['Gene start (bp)'], "txend": df['Gene end (bp)']})

        return df1
        


def _do_scanpy_processing_and_filtering(adata, sample_frac=1.0, do_log1p=True):
    import scanpy as sc

    if sample_frac > 0.999:
        sample_n = int(adata.shape[0]* float(sample_frac))
        sample_idxs = np.random.choice(adata.shape[0], sample_n, replace=False)
        adata2 = adata[sample_idxs].copy()
    else:
        adata2 = adata.copy()
    dbg_print("Flag 123.10 ", adata.shape, adata2.shape, sample_frac)
    
    cinfo, _ = sc.pp.calculate_qc_metrics(adata2, inplace=False)
    c2depth_thresh = sorted(cinfo["total_counts"])[max(10,int(0.001*adata2.shape[0]))]
    c2gcnt_thresh = sorted(cinfo["n_genes_by_counts"])[max(10,int(0.001*adata2.shape[0]))]
    dbg_print("Flag 123.20 ", c2depth_thresh, c2gcnt_thresh)
    
    sc.pp.filter_cells(adata2, min_counts = c2depth_thresh+1)
    sc.pp.filter_cells(adata2, min_genes = c2gcnt_thresh+1)
    sc.pp.normalize_total(adata2, target_sum=1e6)
    if do_log1p:
        sc.pp.log1p(adata2)
    dbg_print("Flag 123.40 ", adata2.shape)
    return adata2


def _convert_adata_to_counts(adata):
    import scanpy as sc
    max1 = np.max(np.max(adata.X))
    if max1 < 15: # likely log values
        adata.X = np.expm1(adata.X)
        sc.pp.normalize_total(adata, target_sum=1e6)
        
    else:
        sc.pp.normalize_total(adata, target_sum=1e6)


        
def _parse_geneids(s):
    l = s.split("|")
    sym = l[0]

    otherid=sym
    if 'HGNC:' in s:
        try:
            otherid = [a.split(':')[2] for a in l if a.startswith("HGNC:")][0]
        except:
            otherid = sym
    if 'MGI:' in s:
        try:
            otherid = [a.split(':')[2] for a in l if a.startswith("MGI:")][0]
        except:
            otherid = sym
    try:
        ensemblid = [a.split(':')[1] for a in l if a.startswith("Ensembl:")][0]
    except:
        ensemblid = sym
        
    return (sym, ensemblid, otherid)




def read_TADMap_from_file_or_url(tad_file_or_url):
    """ Given a file that contains the TAD Map in the default format, parse it and return the results

    The file format needs to match the ones currently hosted at 
    `https://cb.csail.mit.edu/cb/tadmap/TADMap_geneset_{hs,mm}.csv`

    In particular, it is a CSV format of the form:

      <tad_name>,[semi-colon separated list of genes partially/fully contained in the TAD]

    Gene may be specified by multiple pipe-separated identifiers, with first one being the canonical,

    :param tad_file_or_url: string specifying file location or URL

    :returns: the pair (`geneset`, `tad2genelist`)
 
        - geneset: the set of all genes listed in the file. Each gene is specified as a 3-tuple: 
          (canonical_name, Ensembl_name, MGI_or_HGNC_name)

        - tad2genelist: a dictionary from tad_name -> list of genes. 
          Only the canonical_name of the gene is provided in the list
"""
    
    df = pd.read_csv(tad_file_or_url)
    dbg_print("Flag 652.01 ", df.shape)
    df = df[~df["genelist"].isnull()]
    
    df["gene_count"] = df["genelist"].apply(lambda v: len(v.split(';')))
    dbg_print("Flag 652.02 ", df.shape, df['gene_count'].describe())

    
    geneset = set()
    tad2genelist = {}
    for i in range(df.shape[0]):
        tad = df["tad"].iat[i]
        genelist = [_parse_geneids(s) for s in df['genelist'].iat[i].split(';')]
        
        #dbg_print("Flag 652.10 ", tad, len(genelist), genelist[0])        
        tad2genelist[tad] = [a[0] for a in genelist]
        geneset.update(genelist)

    return geneset, tad2genelist



def retrieve_TADMap_by_species(sp_2_letter):
    """ Retrieve the pre-computed TAD MAp from the default remote location

    The default location is     `https://cb.csail.mit.edu/cb/tadmap/TADMap_geneset_{hs,mm}.csv`

    :param sp_2_letter: string, specifying `hs` (human) or `mm` (mouse). These are the only species 
       supported currently.

    :returns:  the same output as `read_TADMap_from_file_or_url`, which is called by this function
""" 
    assert sp_2_letter in ["hs","mm"]
    server_location = {"hs": "http://cb.csail.mit.edu/cb/tadmap/TADMap_geneset_hs.csv",
                       "mm": "http://cb.csail.mit.edu/cb/tadmap/TADMap_geneset_mm.csv",}
    
    return read_TADMap_from_file_or_url(server_location[sp_2_letter])



def standardize_adata_gene_names(adata_in, sp_2_letter):
    """ Process the gene names in an AnnData object to map them to Ensembl v102 names, removing those
    that do not match.

    This is an preprocessing step to call before calling `compute_tad_signature`, 
    since the rest of this module only works with Ensembl v102 canonical names. 

    :param adata_in: the input AnnData object, with gene names in `adata_in.var_names`

    :param sp_2_letter: string specifying the species: `hs` (human) or `mm` (mouse)


    :returns: AnnData object, with var_names mapped to Ensembl v102, duplicates and non-mappable names removed. 
"""
    assert sp_2_letter in ["hs","mm"]
    geneset, _ = retrieve_TADMap_by_species(sp_2_letter)
    return _match_adata_to_geneset(adata_in, geneset)



def _match_adata_to_geneset(adata_in, geneset):
    adata = adata_in.copy()

    f_remove_ens_version = lambda s: s.split('.')[0] if s[:3]=="ENS" else s
    
    D_upper_to_standard  = {}
    for g0, g1, g2 in geneset:
        D_upper_to_standard[g0.upper()] = g0
        D_upper_to_standard[g1.upper()] = g1
        D_upper_to_standard[g2.upper()] = g2

    # remove versions and fix case to match geneset
    genes = list(adata.var_names)
    genes1 = [f_remove_ens_version(g) for g in genes]
    genes2 = [D_upper_to_standard.get(g.upper(),g) for g in genes1]
    adata.var_names = genes2
    genes = list(adata.var_names)
    
    symbol_matchcnt = len( set(genes) & set(a[0] for a in geneset))
    ensembl_matchcnt = len( set(genes) & set(a[1] for a in geneset))
    dbg_print ("Flag 342.24 ", symbol_matchcnt, ensembl_matchcnt, len(geneset), len(genes))
    
    if ensembl_matchcnt > symbol_matchcnt:
        ensembl2gene = {a[1]:a[0] for a in geneset}
        genes1 = copy.deepcopy(genes)
        genes = [ensembl2gene.get(g,g) for g in genes1]
        adata.var_names = genes
        
    geneset_symbols = set(a[0] for a in geneset)
    adata.var_names_make_unique() # because of the matching code below, all but one of the dupes will be dropped
    valid_genes = np.array([(g in geneset_symbols) for g in adata.var_names])

    dbg_print ("Flag 342.29 ", adata.shape)
    return adata[:, valid_genes].copy()



def _read_scvelo_data(rnaseq_data_path, geneset):
    from anndata import AnnData
    import scanpy as sc
    
    adata = sc.read(rnaseq_data_path)
    dbg_print ("Flag 342.20 ", adata.shape, len(geneset))
    adata.var_names = [c.upper() for c in adata.var_names]
    adata.var_names_make_unique() # because of the matching code below, all but one of the dupes will be dropped

    adata = _match_adata_to_geneset(adata, geneset)
    
    dbg_print ("Flag 342.30 ", adata.shape)
    if 'cluster' not in adata.obs.columns:
        if 'clusters_coarse' in adata.obs.columns:
            adata.obs['cluster'] = adata.obs['clusters_coarse']
            
        if 'clusters' in adata.obs.columns:
            adata.obs['cluster'] = adata.obs['clusters']
            
        
    return adata
    
