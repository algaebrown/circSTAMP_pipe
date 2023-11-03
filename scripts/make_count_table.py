from pathlib import Path
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
    data = pd.concat([df, attrs], axis = 1)
    
    data['bsj']=data['bsj'].astype(float)
    data['fsj']=data['fsj'].astype(float)
    data['junc_ratio']=data['junc_ratio'].astype(float)
    
    # stat
    
    
    return data

if __name__=='__main__':
    indir = Path(sys.argv[1])
    outdir = sys.argv[2]

    # read files
    all_files = list(indir.glob('*.gtf'))
    all_dfs = [read_ciri_gtf(f) for f in all_files]
    names = [f.name.split('.')[0] for f in all_files]

    # unify annotations
    cols = ['seqname', 'source', 'feature', 'start', 'end',
        'strand', 'frame', 'circ_id', 'circ_type','gene_id', 'gene_name', 'gene_type']
    anno_df = pd.concat([df[cols] for df in all_dfs], axis = 0).drop_duplicates()

    # BSJ
    bsjs = pd.concat([df.set_index('circ_id')['bsj'] for df in all_dfs], 
                    axis = 1).fillna(0)
    bsjs.columns = [f'BSJ-{name}' for name in names]

    # FSJ
    fsjs = pd.concat([df.set_index('circ_id')['fsj'] for df in all_dfs], 
                    axis = 1).fillna(0)
    fsjs.columns = [f'FSJ-{name}' for name in names]

    # counts
    counts = anno_df.merge(bsjs, left_on = 'circ_id', right_index = True
                        ).merge(fsjs, left_on = 'circ_id', right_index = True)

    counts.to_csv(outdir, sep = '\t')