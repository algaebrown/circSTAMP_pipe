import seaborn as sns
import os
from pybedtools import BedTool
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import sys
from statsmodels.stats.multitest import fdrcorrection
from scipy.special import expit
from scipy.stats import betabinom

def get_alpha_beta(mu, rho):
    #https://dcgerard.github.io/updog/reference/betabinom.html
    
    alpha = mu*(1-rho)/rho
    beta=(1-mu)*(1-rho)/rho
    
    return alpha, beta



if __name__=='__main__':
    counts = pd.read_csv(sys.argv[1], sep = '\t')
    coef = pd.read_csv(sys.argv[2], sep = '\t')
    IP = sys.argv[3]
    IN = sys.argv[4]
    outf = sys.argv[5]
    FDR_threshold = 0.2

    # get rho
    rho = expit(coef['rho'].median()) # the inverse of logit.
    print(rho)

    r1_handle = f'BSJ-{IP}'
    r2_handle = f'BSJ-{IN}'



    total_bsj = counts[[r1_handle, r2_handle]].sum(axis = 0)

    # get beta binomial parameters
    bsj_rate = total_bsj[r1_handle]/total_bsj.sum() # mean
    alpha_bsj, beta_bsj = get_alpha_beta(bsj_rate, rho) # alpha, beta

    # size
    counts['total_bsj']=counts[r1_handle]+counts[r2_handle]

    unique_count_combinations = counts[['total_bsj', r1_handle]].drop_duplicates()
    unique_count_combinations[f'pvalue']=unique_count_combinations.apply(lambda row: 1-betabinom(row['total_bsj'], 
                                            a=alpha_bsj, 
                                            b = beta_bsj).cdf(row[r1_handle]), axis =1)
    
    counts = counts.merge(unique_count_combinations, left_on = ['total_bsj', r1_handle],
             right_on = ['total_bsj', r1_handle])

    
    # multiple hypothesis testing
    counts['perc_read']=(counts[r1_handle]/counts[r1_handle].sum())+(counts[r2_handle]/counts[r2_handle].sum())

    counts['perc_read_bins']=pd.qcut(counts['perc_read'], q = 100, duplicates = 'drop', labels = False)
    

    n_rejected = []
    for nread_bin in counts['perc_read_bins'].value_counts().sort_index(ascending = False).index:
        rejected, fdr = fdrcorrection(counts.loc[counts['perc_read_bins']>nread_bin, f'pvalue'],
                                    alpha = FDR_threshold)
        n_rejected.append([nread_bin, rejected.sum()])
    n_rejected = pd.DataFrame(n_rejected, columns = ['perc_read_bin', 'n_rejected'])
    n_rejected.set_index('perc_read_bin')['n_rejected'].plot()

    threshold = n_rejected.set_index('perc_read_bin')['n_rejected'].idxmax()

    results = counts.loc[counts['perc_read_bins']>threshold]
    results['rejected'], results[f'FDR'] = fdrcorrection(results[f'pvalue'],
                            alpha = FDR_threshold)
    

    counts[f'FDR']=results[f'FDR']
    counts[f'tested']=counts.index.isin(results.index)
    counts[f'fold change']=(counts[r1_handle]/counts['total_bsj'])/(bsj_rate)

    counts[['circ_id', 'pvalue','FDR', 'tested', 'fold change']].to_csv(outf)