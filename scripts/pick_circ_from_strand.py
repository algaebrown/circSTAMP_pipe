import pandas as pd
from pybedtools import BedTool
import sys
from statsmodels.stats.multitest import fdrcorrection

def read_ciri_gtf(fname):
    bed = BedTool(fname)
    df = bed.to_dataframe()
    # filter for non-entries
    stat = df.loc[df['seqname'].str.contains('##'), 'seqname'].str.split(': ', expand = True)
    df = df.loc[~df['seqname'].str.contains('##')].reset_index()
    
    # get attributes
    attrs = pd.DataFrame([i.attrs for i in bed])
    
    return pd.concat([df, attrs], axis = 1), stat

if __name__ == '__main__':
    neg = sys.argv[1]
    pos = sys.argv[2]
    gtf = sys.argv[3]
    outf = sys.argv[4]

    neg_df = pd.read_csv(neg, sep = '\t', index_col = 0)
    pos_df = pd.read_csv(pos, sep = '\t', index_col = 0)
    circ_anno, stat = read_ciri_gtf(gtf)

    neg_df['POS']=neg_df['pos_id'].apply(lambda s: int(s.split(':')[-1]))
    neg_df['circ_id']=neg_df['pos_id'].apply(lambda s: ':'.join(s.split(':')[:-1]))

    pos_df['POS']=pos_df['pos_id'].apply(lambda s: int(s.split(':')[-1]))
    pos_df['circ_id']=pos_df['pos_id'].apply(lambda s: ':'.join(s.split(':')[:-1]))

    neg_df['strand']='-'
    pos_df['strand']='+'

    pos_circ = circ_anno.loc[circ_anno['strand']=='+', 'circ_id']
    neg_circ = circ_anno.loc[circ_anno['strand']=='-', 'circ_id']

    pos_df = pos_df.loc[pos_df['circ_id'].isin(pos_circ)]
    neg_df = neg_df.loc[neg_df['circ_id'].isin(neg_circ)]

    edit_df = pd.concat([pos_df, neg_df])

    _,edit_df['FDR'] = fdrcorrection(edit_df['pvalue'])
    edit_df.to_csv(outf, sep = '\t')
