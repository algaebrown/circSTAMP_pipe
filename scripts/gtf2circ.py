from pybedtools import BedTool
import pandas as pd
import sys

def read_ciri_gtf(fname):
    bed = BedTool(fname)
    df = bed.to_dataframe()
    # filter for non-entries
    stat = df.loc[df['seqname'].str.contains('##'), 'seqname'].str.split(': ', expand = True)
    df = df.loc[~df['seqname'].str.contains('##')].reset_index()
    
    # get attributes
    attrs = pd.DataFrame([i.attrs for i in bed])
    
    return pd.concat([df, attrs], axis = 1), stat

if __name__=='__main__':
    input_gtf = sys.argv[1]
    output_circ = sys.argv[2]

    circ_df, stat = read_ciri_gtf(input_gtf)
    for c in ['bsj', 'start', 'end']:
        circ_df[c] = circ_df[c].astype(float).astype(int)
    circ_df[['seqname', 'start', 'end', 'bsj']].to_csv(output_circ, 
                                                    sep = '\t', header = None, index = None)
