import pandas as pd
from pybedtools import BedTool
import sys

def find_edited_windows(df, size = 50, fdr_thres = 0.2):
    ''' convert significant edits into bed file of genomic corredinates 
    expand edits, then test
    '''
    df['chrom']=df['circ_id'].apply(lambda s: s.split(':')[0])
    df['end']=df['circ_id'].apply(lambda s: int(s.split('|')[-1]))
    df['start']=df['circ_id'].apply(lambda s: int(s.split('|')[0].split(':')[1]))
    df['abs_pos']=df['start']+df['POS']

    df['edit_start']=df['abs_pos']-size
    df['edit_end']=df['abs_pos']+size

    from pybedtools import BedTool
    b = BedTool.from_dataframe(df.loc[df['FDR']<fdr_thres,
    ['chrom', 'edit_start', 'edit_end', 'pos_id', 'POS', 'strand']])
    
    return b

if __name__ == '__main__':
    df = pd.read_csv(sys.argv[1], index_col= 0, sep = '\t')
    bed = find_edited_windows(df)

    bed.saveas(sys.argv[2])