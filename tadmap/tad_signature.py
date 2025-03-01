#!/usr/bin/env python

###################################################################
## Primary Author:  Rohit Singh rsingh@alum.mit.edu
## License: MIT
## Repository:  http://github.io/rs239/tadmap
###################################################################

import pandas as pd
import numpy as np
import scipy, os, sys, string, fileinput, glob, re, math, itertools, functools, copy, multiprocessing, traceback, tarfile, gzip, csv, tqdm
import scipy.stats
from scipy.stats import describe
from scipy import sparse
import os.path
import scipy.sparse
from scipy.sparse import csr_matrix, csc_matrix
from collections import defaultdict

from . import tadmap_base_config
dbg_print = tadmap_base_config.dbg_print

from . import tadmap_utils as U

#import .tadmap_utils as U


def _score_tads_by_alphas_poisson(tEXPRs, tCapacity, lambda1, lambda2):
    #https://www.statlect.com/fundamentals-of-statistics/Poisson-distribution-maximum-likelihood
    #https://people.stat.sc.edu/Hitchcock/slides535day5spr2014.pdf
    
    f_ll = lambda L:   -tCapacity*L + np.log(L)*tEXPRs
    
    dbg_print("Flag 616.10 ", tEXPRs.shape, tCapacity.shape, lambda1, lambda2)
    
    l1, l2 = f_ll(lambda1), f_ll(lambda2)
    dbg_print("Flag 616.30 ", l1.shape, l2.shape, l1[:3,:3], l2[:3,:3])
    
    df_mle_prob = 1.0/ (np.exp(l2-l1) + 1)  # a/(a+b) = 1/ (np.exp(log(b)-log(a)) + 1)
    dbg_print("Flag 616.40 ", df_mle_prob.shape, df_mle_prob[:3,:3])
    
    return df_mle_prob




def _map_sc_to_tad_poisson(adata, tad2genelist, skip_singleton_tads=False):
    dbg_print("Flag 874.10 ", adata.shape, len(tad2genelist))
    l1 = sorted([ (t, len(v)) for t,v in tad2genelist.items()], key=lambda s:int(s[0].split('|')[0]))
    
    if skip_singleton_tads:
        l0 = [(a,b) for a,b in l1 if b>1]
    else:
        l0 = [(a,b) for a,b in l1 if b>0]
        
    dbg_print("Flag 874.15 ", len(l0), l0[:3])
    
    tadCnt = pd.Series( [a[1] for a in l0], index=[a[0] for a in l0])
    nT = len(tadCnt)
    tadNames = list(tadCnt.index)
    dbg_print("Flag 874.20 ", nT, tadNames[:3], tadCnt.head(3))
    
    gene2tadidx = defaultdict(list)
    for i, (t,_) in enumerate(l0):
        for g in tad2genelist[t]:
            if g in gene2tadidx:
                dbg_print("Flag 874.22 INFO: gene {1} seems to span multiple tads. Already seen {0}, now seeing {2}".format([tadNames[a] for a in gene2tadidx[g]], g, t))
            gene2tadidx[g].append(i)

    dbg_print("Flag 874.25 ", len(gene2tadidx))
    
    adata_geneidx2tadidx = [[] for i in range(adata.shape[1])] #-1*np.ones(adata.shape[1])
    for i, g in enumerate(adata.var_names):
        adata_geneidx2tadidx[i] = gene2tadidx[g]

    adata_tadidx2geneidx = [[] for i in range(nT)]
    for i, tidxlist in enumerate(adata_geneidx2tadidx):
        for tidx in tidxlist:
            adata_tadidx2geneidx[tidx].append(i)
    
        
    dbg_print("Flag 874.30 ", len(adata_geneidx2tadidx), len(adata_tadidx2geneidx), adata_geneidx2tadidx[:5], adata_tadidx2geneidx[:5])
        
    if scipy.sparse.isspmatrix_csr(adata.X):
        X1 = adata.X.todense()
    else:
        X1 = adata.X
        
    df_tCapacity = np.tile( tadCnt.values, (adata.shape[0],1))

    df_tEXPRs = np.full((adata.shape[0], nT),0.0)
    for tidx in range(nT):
        dbg_print("Flag 874.75 ", tidx, df_tCapacity.shape, df_tEXPRs.shape, df_tEXPRs.sum(), X1.shape)
        if len(adata_tadidx2geneidx[tidx]) == 0: continue
        gidxlist = adata_tadidx2geneidx[tidx]
        dbg_print("Flag 874.76 ", df_tEXPRs[:,tidx].shape, X1[:, gidxlist].sum(axis=1).shape)
        df_tEXPRs[:,tidx] += np.ravel(X1[:, gidxlist].sum(axis=1)) 
            
    dbg_print("Flag 874.80 ", df_tEXPRs.shape, df_tCapacity.shape, df_tEXPRs[:3,:3], df_tCapacity[:3,:3])
        
    return df_tEXPRs, df_tCapacity, tadCnt, tadNames, gene2tadidx, adata_geneidx2tadidx





#https://people.stat.sc.edu/Hitchcock/slides535day5spr2014.pdf
def _score_tads_by_EM_poisson(tEXPRs, tCapacity, gammaDistHyperParams=None):
    
    mean_umi_rate = np.mean(tEXPRs/tCapacity)
                             
    dbg_print("Flag 635.20 ", tEXPRs.shape, tCapacity.shape, mean_umi_rate)
    
    lambda1 = 2*mean_umi_rate; lambda2 = 0.2*mean_umi_rate

    gammaA1, gammaA2, gammaB = 0.5*mean_umi_rate, 0, 0.5*mean_umi_rate #3, 0, 3
    if gammaDistHyperParams is not None:
        gammaA1, gammaA2, gammaB = gammaDistHyperParams

    dbg_print("Flag 635.22 ", gammaA1, gammaA2, gammaB)
    
        
    newlambda1 = newlambda2 = -1
    itercnt = 0; EPS=0.00001
    
    # mimic with tqdm: while itercnt < 100 and (abs(newlambda1-lambda1) > EPS  or abs(newlambda2-lambda2) > EPS):
    
    tqdm_iter = tqdm.tqdm(range(100))
    for itercnt in tqdm_iter:
        if (abs(newlambda1-lambda1) <= EPS  and abs(newlambda2-lambda2) <= EPS):
            tqdm_iter.close()
            break
        
        if itercnt>0:
            lambda1 = newlambda1
            lambda2 = newlambda2
            
        df_tad_probs = _score_tads_by_alphas_poisson(tEXPRs, tCapacity, lambda1, lambda2)

        #tU1, tU2 = (tEXPRs + gammaA1), (tEXPRs + gammaA2)
        tU1, tU2 = (tEXPRs + gammaA1*tCapacity), (tEXPRs + gammaA2*tCapacity)
        tC = tCapacity + gammaB
        
        newlambda1 = np.sum(np.ravel(df_tad_probs*tU1))/np.sum(np.ravel(df_tad_probs*tC))
        newlambda2 = np.sum(np.ravel(df_tad_probs*tU2))/np.sum(np.ravel(df_tad_probs*tC))


        itercnt += 1
        dbg_print("Flag 635.40 ",df_tad_probs.shape, lambda1, lambda2, newlambda1, newlambda2, EPS, itercnt)

    if lambda1 < lambda2:
        sys.stderr.write("WARNING: EM algorithm may not have converged correctly (%g, %g)\n" % (lambda1, lambda2))
    
    return lambda1, lambda2, df_tad_probs




def _compute_tad_occupancy_by_EM_poisson( adata, tad2genelist, extra_args):

    tadmap_base_config.tadmap_info("Checking adata...")
    
    if int(extra_args.get("adata_is_logcpm_normalized",0))  < 0.5:
        U._convert_adata_to_counts(adata)
        
        #did the log transform to limit the variability of the data
        adata.X = np.round_(np.log1p(adata.X)).astype(int).astype(float)
    else:
        adata.X = np.round_(adata.X).astype(int).astype(float)

        
    dbg_print("Flag 935.20 ", adata.X.shape, adata.X[:7,:7].todense())
    
    gammaDistHyperParams = None
    if "gammaprior_hyperparams" in extra_args:
        l = [float(a) for a in extra_args.get("gammaprior_hyperparams").split(",")]
        gammaDistHyperParams = (l[0], l[1], l[2])
        dbg_print("Flag 935.625 ", gammaDistHyperParams)
        
    tadmap_base_config.tadmap_info("Mapping gene expression to TADs...")
    
    adata_tadEXPRs, adata_tadCapacity, tadCnt, tadNames, _, _ = _map_sc_to_tad_poisson(adata, tad2genelist)    
    dbg_print("Flag 935.65 ", adata_tadEXPRs.shape, adata_tadCapacity.shape)
    
    tadmap_base_config.tadmap_info("Running EM...")
    lambda1, lambda2, _ = _score_tads_by_EM_poisson(adata_tadEXPRs, adata_tadCapacity, gammaDistHyperParams)    
    dbg_print("Flag 935.67 ",  lambda1, lambda2)
    
    tad_occupancy_mtx = _score_tads_by_alphas_poisson(adata_tadEXPRs, adata_tadCapacity, lambda1, lambda2)

    return lambda1, lambda2, tadCnt, tadNames,  tad_occupancy_mtx    





def compute_tad_signature(adata, sp_2_letter):
    """Given an AnnData object and a species (`hs` or `mm`), compute a TAD activation profile for each cell

    The activation profile is computed by fitting a 2-component Poisson mixture model using the
    Expectation Maximization (EM) algorithm. One component corresponds to TADs that are transcriptionally 
    active (i.e., "ON"), while the other corresponds to "OFF" TADs. However, even "OFF" TADs can 
    have genes with active expression (e.g. isolated expression of a single gene in a 
    non-TAD-dependent fashion). For each cell, the EM algorithm computes for each TAD --- there are 
    approx 3000 of them in human and mouse --- the probability that the TAD is in "ON" state.


    :param adata:
        AnnData object (`n` cells). The gene expression matrix can be sparse or dense 
        and contain counts or log1p-transformed data--- the method will try to adapt accordingly.

    :type adata: AnnData object

    :param sp_2_letter: one of 'hs' or 'mm'
        Currently, TAD Maps are supported only for human ('hs') or mouse ('mm')

    :type sp_2_letter: string

    :returns: a pair of Pandas dataframes: the TAD signature and auxiliary information, respectively

        The first dataframe is of dimensionality n x T where n is the number of cells and T is the
        number of TADs. The algorithm will filter out TADs which had no active genes in the data so
        T may vary a little across datasets. The column names correspond to TAD names which are in
        the following format: `<numeric_id>|<chromosome>|<start>|<end>` 

        The second dataframe contains one row per (TAD,gene) pair. Some genes may span two TADs and
        will have two rows. Each row contains the TADs score dispersion, an indication of its 
        variability, similar to highly variable genes. Specifically, here are the column names of 
        this dataframe:

             - `tad_name`: see above
             - `tad_gene_count`: number of protein-coding genes partially/fully contained in the TAD
             - `lambda_ON`: lambda for the Poisson corresponding to "ON" TADs (same for all rows)
             - `lambda_OFF`: lambda for the Poisson corresponding to "OFF" TADs (same for all rows)
             - `tad_activation_mean`: the average probability score of activation for this TAD 
                 across `n` cells
             - `tad_activation_disp`: variance/mean of probability score of activation for this 
                 TAD across `n` cells. Use this as the measure for identifying highly variable TADs
             - `gene`: Ensembl v102 name of a gene contained in the TAD. 
                 There is one row for each (gene,TAD) pair


"""

    assert sp_2_letter in ["hs","mm"]
    
    extra_args = {}
    adata2 = adata.copy()
    return _compute_tad_signature(adata2, sp_2_letter, extra_args)



def _compute_tad_signature(adata, sp, extra_args):
    assert sp in ["hs","mm"]
    
    geneset, tad2genelist = U.retrieve_TADMap_by_species(sp)

    adata_genes = set(list(adata.var_names))
    l2 = set([g[0] for g in geneset])
    if len( adata_genes & l2) != len( adata_genes):
        raise ValueError("Process with standardize_adata_gene_names(...) first. Found unrecognized genes in adata.var_names")

    dbg_print("Flag 934.05 ", adata.shape, len(geneset), len(tad2genelist))
    dbg_print("Flag 934.20 ", adata.shape, adata.obs_names[:5])
    
    lambda1, lambda2, tadCnt, tadNames,  tad_occupancy_mtx = _compute_tad_occupancy_by_EM_poisson( adata, tad2genelist, extra_args)
        
    dbg_print("Flag 934.67 ",  tad_occupancy_mtx.shape, lambda1, lambda2)
    tad_occupancy_df = pd.DataFrame(tad_occupancy_mtx, columns = tadNames)
    dbg_print("Flag 934.70 ", tad_occupancy_df.shape, lambda1, lambda2)    
    
    # dbg_print the variation in TAD values
    mu_tad = tad_occupancy_df.mean(axis=0)
    disp_tad = tad_occupancy_df.var(axis=0)/(1e-12 + tad_occupancy_df.mean(axis=0))

    tad_aux = pd.DataFrame([[tadNames[i], tadCnt[i], lambda1, lambda2] for i in range(len(tadNames))],
                            columns = ["tad_name", "tad_gene_count", "lambda_ON", "lambda_OFF"])
    
    tad_aux["tad_activation_mean"] = mu_tad.values
    tad_aux["tad_activation_disp"] = disp_tad.values
    tad_aux.loc[mu_tad.values < 1e-12, "tad_activation_disp"] = np.NaN

    l3 = []
    for t, glist in tad2genelist.items():
        for g in glist:
            l3.append((t, g))
            
    tad_aux = tad_aux.merge(pd.DataFrame(l3, columns=["tad_name","gene"]), how="left").reset_index(drop=True)

    tad_occupancy_df = tad_occupancy_df.loc[:, (mu_tad.values > 1e-12)]
    tad_aux = tad_aux[ tad_aux['tad_activation_mean'] > 1e-12 ].reset_index(drop=True)

    return tad_occupancy_df, tad_aux


def to_log_odds(tad_occupancy_df):
    """ Convert probability scores `p` to log-odds, `log(p/(1-p))`

    This is a useful conversion to do before passing TAD signatures to a clustering or visualization 
    process. It widens the range of values and makes them more compatible with the Euclidean 
    distance metric, which underlies many clustering and visualization algorithms.

    :param tad_occupancy_df:  Pandas dataframe

    This is the first dataframe item in the pair of dataframes returned by `compute_tad_signature`

    :returns: Pandas dataframe, same dimensions as the input

"""
    M = np.log(np.maximum(tad_occupancy_df.values, 1e-15) / np.maximum(1-tad_occupancy_df.values, 1e-15))
    return pd.DataFrame(M, index=tad_occupancy_df.index, columns=tad_occupancy_df.columns)



# def _DEPRECATED_compute_tad_representation(sp, outdir, outsfx, tad_file, rnaseq_data_path, rnaseq_data_type, extra_args):
#     dbg_print("Flag 934.01 ", sp, outdir, outsfx, tad_file, rnaseq_data_path, rnaseq_data_type, extra_args)
    
#     # read tads. also get a list of valid genes
#     geneset, tad2genelist = U.retrieve_TADMap_by_species(sp)
#     dbg_print("Flag 934.05 ", len(geneset), len(tad2genelist))


#     # read rnaseq data to scanpy. supply list of valid genes to subset
#     if rnaseq_data_type == 'trajectorama':
#         adata = U.read_trajectorama_data(rnaseq_data_path, geneset)
#     elif rnaseq_data_path == 'scvelo':
#         adata = U.read_scvelo_data(rnaseq_data_path, geneset)
#     else:
#         assert 1==0

#     dbg_print("Flag 934.20 ", adata.shape, adata.obs_names[:5])
    
#     alpha1, alpha2, tadCnt, tadNames,  tad_occupancy_mtx = _compute_tad_occupancy_by_EM_poisson( adata, tad2genelist, extra_args)

        
#     dbg_print("Flag 934.67 ",  tad_occupancy_mtx.shape, alpha1, alpha2)
#     tad_occupancy_df = pd.DataFrame(tad_occupancy_mtx, columns = tadNames)
#     dbg_print("Flag 934.70 ", tad_occupancy_df.shape, alpha1, alpha2)
    
    
#     # dbg_print the variation in TAD values
#     mu_tad = tad_occupancy_df.mean(axis=0)
#     z_tad = np.sqrt(tad_occupancy_df.shape[0])*tad_occupancy_df.mean(axis=0)/tad_occupancy_df.std(axis=0)

#     n1 = tad_occupancy_df.shape[1]
#     l = [['alpha1']+[alpha1]*n1, ['alpha2'] +[alpha2]*n1, ['tadmemb_count'] + list(tadCnt.values),  ['mu']+list(mu_tad.values), ['z']+list(z_tad.values)]
#     tad_aux = pd.DataFrame(l, columns = ['name']+list(tad_occupancy_df.columns))

#     dbg_print("Flag 934.80 ", tad_aux.shape)

#     outfile1 = "{0}/tad_occupancy_{1}.h5".format(outdir, outsfx)
#     outfile2 = "{0}/tad_auxinfo_{1}.csv".format(outdir, outsfx)
#     tad_occupancy_df.to_hdf(outfile1, key="df", index=False)
#     tad_aux.to_csv(outfile2, index=False)
    
#     return tad_occupancy_df, tad_aux


    
        
############################################################################

if __name__ == "__main__":
    sys.path.append(os.path.join(sys.path[0],PROJ_DIR+'/src'))


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--species", help="2-letter code of species (hs or mm)", type=str, default='hs')
    parser.add_argument("--outdir", help="output directory (can set to '.')", type=str, default=PROJ_DIR+"/data/processed/")
    parser.add_argument("--outsfx", help="suffix to use when producing output files")
    parser.add_argument("--tad_file", help="the path to the TAD file", type=str)
    parser.add_argument("--rnaseq_data_path", help="the path to the rnaseq_data_path", type=str)
    parser.add_argument("--rnaseq_data_type", help="one of trajectorama|cytotrace", type=str, choices=("trajectorama","cytotrace"), default="trajectorama")
    parser.add_argument("--extra", help="put this as the LAST option and arbitrary space-separated key=val pairs after that", type=str, nargs='*')


    args = parser.parse_args()
    extra_args = dict([a.split("=") for a in args.extra]) if args.extra else {}

    tad_occupancy_df, tad_aux = run(args.species, args.outdir, args.outsfx, args.tad_file, args.rnaseq_data_path, args.rnaseq_data_type, extra_args)

