from collections import defaultdict
import sys
import pandas as pd

if __name__ == '__main__':
    fname = sys.argv[1] # tsv file from gatk
    alt = sys.argv[2]
    outf = sys.argv[3]

    ref_count=defaultdict(lambda:0)
    alt_count=defaultdict(lambda:0)
    
    
    for df in pd.read_csv(fname,
                sep = '\t', chunksize = 100000):
        ad_col = df.columns[df.columns.str.endswith('.AD')][0]
        for circ_id, group in df.groupby(by = 'CHROM'):

            g_sum = group[ad_col].apply(lambda s: int(s.split(',')[0])).sum()

            with_alt = group.loc[group['ALT'].str.contains(f'{alt},')]
            alt_sum = with_alt[ad_col].apply(lambda s: int(s.split(',')[1])).sum()

            ref_count[circ_id]+=g_sum
            alt_count[circ_id]+=alt_sum
    
    cnt_df = pd.concat([pd.Series(alt_count), pd.Series(ref_count)], axis = 1)
    cnt_df.columns = ['ALT', 'REF']

    cnt_df.to_csv(outf)