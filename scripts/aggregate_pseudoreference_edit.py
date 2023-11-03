import numpy as np
from scipy.stats import binom
from collections import defaultdict
import pandas as pd
import sys

"""
python ~/projects/circSTAMP_pipe/scripts/aggregate_pseudoreference_edit.py /home/hsher/scratch/circ_nextera_iter13/output/edits/RBM15_STAMP.dp4.neg.vcf.tsv A /home/hsher/scratch/circ_nextera_iter13/output/edits/RBM15_STAMP.dp4.neg.combined.vcf.tsv
"""

def aggregate_counts(df):
    # the reference is repeating the sequence twice to represent BSJ
    df['length']=df['CHROM'].apply(lambda string: int(string.split('|')[1])-int(string.split('|')[0].split(':')[1])+1)
    df['POS']=df['POS']%df['length']

    ad_col = df.columns[df.columns.str.endswith('.AD')][0]

    df['n_ref'] = df[ad_col].apply(lambda s: int(s.split(',')[0]))
    df['n_alt'] = df[ad_col].apply(lambda s: int(s.split(',')[1]))
    df['pos_id']=df['CHROM']+':'+df['POS'].astype(str)

    aggregated_df = df.groupby(by = 'pos_id')[['n_ref', 'n_alt']].sum()
    aggregated_df = aggregated_df.loc[aggregated_df.sum(axis = 1)>0]
    return aggregated_df

if __name__ == '__main__':
    fname = sys.argv[1]
    alt = sys.argv[2]
    outf = sys.argv[3]
    outf_nonzero = outf.replace('.tsv', '.nonzero.tsv')



    previous_chunk = pd.DataFrame()
    i = 0
    with open(outf, 'w') as fhandle:
        with open(outf_nonzero, 'w') as fhandle_nonzero:
            for df in pd.read_csv(fname,
                        sep = '\t', chunksize = 1000000):
                
                df = df.loc[(df['ALT']=='<*>')|(df['ALT']==f'{alt},<*>')]

                last_circ = df['CHROM'].iloc[-1]
                process_next_chunk = df.loc[df['CHROM']==last_circ].copy()

                # process the last circular RNA in next chunck, and add back the ones from previous chunk
                df = df.loc[df['CHROM']!=last_circ]
                
                df = pd.concat([df,previous_chunk], axis = 0)
                print(df.shape[0])

                aggregated_df = aggregate_counts(df)
                if i == 0:
                    aggregated_df.to_csv(fhandle, sep = '\t', index = True, header = True)
                    aggregated_df.loc[aggregated_df['n_alt']>0].to_csv(
                        fhandle_nonzero, sep = '\t', index = True, header = True)
                else:
                    aggregated_df.to_csv(fhandle, sep = '\t', index = True, header = False)
                    aggregated_df.loc[aggregated_df['n_alt']>0].to_csv(
                        fhandle_nonzero, sep = '\t', index = True, header = False)
                print(f'''
                processed chunk {i}:
                found sites with edit:
                ''',
                aggregated_df.loc[aggregated_df['n_alt']>0].shape[0]
                )
                previous_chunk = process_next_chunk.copy()
                i+=1

            # process the last chunk
            aggregated_df = aggregate_counts(previous_chunk)
            aggregated_df.to_csv(fhandle, sep = '\t', index = True, header = False)
            aggregated_df.loc[aggregated_df['n_alt']>0].to_csv(
                fhandle_nonzero, sep = '\t', index = True, header = False)