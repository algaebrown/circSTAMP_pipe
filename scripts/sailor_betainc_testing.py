import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import betainc
from scipy.special import betainc
from statsmodels.stats.multitest import fdrcorrection
import sys
def process(num_reads, n_edited, alfa=0, beta=0, cov_margin=0.01, keep_all_edited=True):
    """
    Calculates, for a single line in a VCF formatted file, the
    confidence score based on depth of coverage and edit fraction %.

    :param line: string
        single vcf formatted line.
    :return confidence: float
        confidence value of the line
    :return return_string: basestring
        full vcf formatted line with confidence
    """
    
    edit_frac = n_edited / float(num_reads)

    # calc smoothed counts and confidence
    destination_smoothed = n_edited + alfa
    origin_smoothed = (num_reads-n_edited) + beta
    theta = destination_smoothed / float(destination_smoothed + origin_smoothed)

    ########  MOST IMPORTANT LINE  ########
    # calculates the confidence of theta as
    # P( theta < cov_margin | A, G) ~ Beta_theta(G, A)
    confidence = 1 - betainc(destination_smoothed, origin_smoothed, cov_margin)
    
    


    return edit_frac, theta, confidence



def call_edits(all_counts, cov_margin=0.05):
    for index, row in all_counts.iterrows():
        edit_frac, theta, confidence = process(row['n_ref']+row['n_alt'], row['n_alt'], cov_margin = cov_margin)
        all_counts.loc[index, 'edit_frac'] = edit_frac
        all_counts.loc[index, 'theta'] = theta
        all_counts.loc[index, 'confidence'] = confidence
        all_counts['confidence'] = all_counts['confidence'].fillna(0)
        all_counts['pvalue']=1-all_counts['confidence']
    
    return all_counts

if __name__ == '__main__':
    f=sys.argv[1]
    outf =sys.argv[2]

    neg = pd.read_csv(f, sep = '\t', 
                    )
    print(neg.columns)

    unique_counts = neg[['n_ref', 'n_alt']].drop_duplicates().copy()
    unique_counts = call_edits(unique_counts)

    neg = neg.merge(unique_counts, left_on = ['n_ref', 'n_alt'],right_on = ['n_ref', 'n_alt'])
    _,neg['FDR'] = fdrcorrection(neg['pvalue'])

    neg.to_csv(outf, sep = '\t')