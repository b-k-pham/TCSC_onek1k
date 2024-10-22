#!/usr/bin/env python
import argparse as ap
import sys
import numpy as np
import pandas as pd
import scipy.linalg as linalg
from numpy.linalg import multi_dot as mdot
from pandas_plink import read_plink
from scipy import stats
from sklearn import linear_model as lm
import random
import warnings
import scipy
from sklearn.model_selection import KFold
from sklearn.model_selection import cross_val_score
#from sklearn.utils.testing import ignore_warnings #not loading 
from sklearn.exceptions import ConvergenceWarning
import time
import subprocess, os
import gzip
import os.path  #new
from os import path  #new
from statsmodels.regression.linear_model import OLS

mvn = stats.multivariate_normal

from pathlib import Path

def sim_geno(L, n):  #can we sim geno for 100K people using a 1Kg reference panel? we can get p from L beacuse it's how many SNPs are on chr1
    """
    removed first argument with is L (works with regular LD matrix too), added p to indicate number of snps.
    Sample genotypes from an MVN approximation.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param n: int the number of genotypes to sample

    :return: numpy.ndarray n x p centered/scaled genotype matrix
    """
    p, p = L.shape
    Z = L.dot(np.random.normal(size=(n, p)).T).T    
    #Z = np.random.normal(size=(n, p)).T #snps x people
    #Z = Z.T #I added this. people x snps 
    Z -= np.mean(Z, axis=0)
    Z /= np.std(Z, axis=0) #173 x 4 in original, ine too after transpose. 
    return Z

def sim_gwas(ngwas, Z_gwas, b_qtls, ve):
    """
    Simulate a GWAS using `ngwas` individuals such that genetics explain `var_explained` of phenotype.

    :param L: numpy.ndarray lower cholesky factor of the p x p LD matrix for the population
    :param ngwas: int the number of GWAS genotypes to sample
    :param b_qtls: numpy.ndarray latent eQTL effects for the causal gene
    :param var_explained: float the amount of phenotypic variance explained by genetic component of gene expression

    :return: (pandas.DataFrame, float) estimated GWAS beta and standard error, causal GE effect
    """
    if ve > 0:
        gwas_expr = np.dot(Z_gwas, b_qtls)
        y1 = sim_trait(gwas_expr, ve)[0]
        y2 = sim_trait(gwas_expr, ve)[0]
    y1 = np.array(y1)
    y2 = np.array(y2)
    trait1 = np.zeros([1,ngwas])
    trait1[0,:] = y1
    trait2 = np.zeros([1,ngwas])
    trait2[0,:] = y2
    return (trait1, trait2) 

def sim_trait(g, h2g): #use this twice, first simulate GE as a trait in the GE cohort, then simulate complex trait based on GE in the GWAS cohort 
    """
    Simulate a complex trait as a function of latent genetic value and env noise.

    :param g: numpy.ndarray of latent genetic values
    :param h2g: float the heritability of the trait in the population

    :return: (numpy.ndarray, float) simulated phenotype, sd of Y
    """
    n = len(g)
    if h2g > 0:
        s2g = np.var(g, ddof=1)
        s2e = s2g * ( (1.0 / h2g ) - 1 )
        e = np.random.normal(0, np.sqrt(s2e), n)
        y = g + e #appears to add noise to the genetic values, adding the same about of noise if the seed is the same 
    else:
        e = np.random.normal(0, 1, n)
        y = e
    # standardize
    y -= np.mean(y)
    y_std = np.std(y)
    y /= y_std
    y_std = y_std.item()
    return y, y_std #should we add these together to get a complex trait? 

# SIMULATE GWAS TRAIT

sim = int(sys.argv[1])
print(sim)
VE = float(sys.argv[2])
print(VE)
ngwas = int(sys.argv[3])

##note: simulate gwas cohort to be the same every time or to depend on sim, but cannot depend on gene. 

nvar = 5
nGenesTot = 1000
#ngwas = 10000 
ti = 0  
filename = "eQTLs_1000sims_Nov/Nov_sims_eQTLeffectsizes_0.75rg_SNPpolygen%s_tissue%s_sim%s.txt.gz" % (nvar,ti+1,sim) #SNP polygen = 5, will vary this: 1, 2, 3, 4, 5, was May. flag
b_qtls_df = pd.read_csv(filename, sep = "\t", header = None) #flag
b_qtls = np.array(b_qtls_df) #(position_matrix[,i],betas[,i],snp_upstream,snp_downstream) #these are 1 indexed. 

filename = "Simulated_GWAS_Cohort_%s.txt.gz" % ngwas
z_gwas = np.array(pd.read_csv(filename, sep = "\t", header = None))
print(z_gwas.shape)
print("simulated full GWAS cohort")

filename = "Nov_Y1_alpha_nCausalGenes100.txt.gz" #need to have nonzero alpha for only 100 genes , gene polygen
Y1_alpha = pd.read_csv(filename, sep = "\t", header = None)
Y1_alpha = np.array(Y1_alpha)

genenum = []
for i in range(1,nGenesTot+1):
    a = [i]*nvar
    genenum = np.concatenate((genenum,a))

#import code
#code.interact(local=locals()) #remember need to activate python 3 otherwise libraries don't work. 

for gene in range(1,nGenesTot+1): #now, we can change eqtl effect sizes to only include genes with < 900 snps nearby. 
    print(gene)
    np.random.seed(int(sim*gene))
    Y1_alpha_chunk = Y1_alpha[gene-1,sim-1]
    if ((Y1_alpha_chunk != 0) and (VE != 0)):
        print("simulating snp-gene-trait effect")
        x = np.where(genenum == gene)[0]
        min_snp = int(b_qtls[np.min(x),2]) #make 1-indexed values 0-indexed for python.
        max_snp = int(b_qtls[np.max(x),3])
        p_int = int(max_snp-(min_snp-1))
        Z_gwas = z_gwas[:,(min_snp-1):max_snp]
        b_qtls_pass = np.zeros([p_int,1]) #causal 
        pos = b_qtls[x,0]
        pos = pos.astype(int) 
        sequence = np.array(range(min_snp,max_snp+1))
        x_var = []
        for j in pos: 
            a = np.where(sequence == j)[0]
            x_var = np.concatenate((x_var,a))
        x_var = x_var.astype(int)
        b_qtls_pass[x_var,0] = b_qtls[x,1] #need new x-var where chunk is the gene. 
        keep = x_var #np.where(b_qtls_pass[:,0] != 0)[0]
        y1, y2 = sim_gwas(ngwas,Z_gwas[:,keep],b_qtls_pass[keep,0],VE)
        print(y1.shape)
        filename = "genetraits/Nov_0.75_ggcoreg/%s_vgene%s/genetraitY1_sim%s.txt.gz" % (ngwas,VE,sim) #gene traits should be correlated between Y1 and Y2 with value = variance of trait. +/- 0.2 depends on alpha sign
        Path(filename).parent.mkdir(parents=True, exist_ok=True)
        y1 = pd.DataFrame(y1)
        if gene == 1:
            y1.round(4).to_csv(filename, sep="\t", index=False, header = False)
        else:
            y1.round(4).to_csv(filename, sep="\t", index=False,mode = "a", header = False)
        filename = "genetraits/Nov_0.75_ggcoreg/%s_vgene%s/genetraitY2_sim%s.txt.gz" % (ngwas,VE,sim)
        y2 = pd.DataFrame(y2)
        if gene == 1:
            y2.round(4).to_csv(filename, sep="\t", index=False, header = False)
        else:
            y2.round(4).to_csv(filename, sep="\t", index=False,mode = "a", header = False)

