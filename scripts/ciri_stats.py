from pybedtools import BedTool
import pandas as pd
import sys
from pathlib import Path
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
    circ_quant_output = sys.argv[1].split(' ')
    outdir = Path(sys.argv[2])

    circ_type_count = []
    names = []
    junc_ratio_tbl = []
    bsj_ratio_tbl = []
    fsj_ratio_tbl = []
    stats = []
    mega_anno = []

    features = ['seqname', 'start', 'end', 
        'strand', 'circ_type', 'gene_id', 'gene_name', 'gene_type']
    for fname in circ_quant_output:
        name = Path(fname).name.replace('.gtf', '')
        names.append(name)
        circ_df, stat = read_ciri_gtf(fname)
        circ_type_count.append(circ_df['circ_type'].value_counts())
        junc_ratio_tbl.append(circ_df.set_index('circ_id')['junc_ratio'].astype(float))
        bsj_ratio_tbl.append(circ_df.set_index('circ_id')['bsj'].astype(float))
        fsj_ratio_tbl.append(circ_df.set_index('circ_id')['fsj'].astype(float))
        stats.append(stat.set_index(0))
        mega_anno.append(circ_df.set_index('circ_id')[features])

    # every single circle found in this dataset. collapsed into one file
    mega_anno = pd.concat(mega_anno, axis = 0).drop_duplicates()

    # basic statistics
    stats_df = pd.concat(stats, axis = 1).T
    for col in ['##Total_Reads', '##Mapped_Reads', '##Circular_Reads']:
        stats_df[col] = stats_df[col].astype(int)

    stats_df = stats_df.drop_duplicates('##Sample').set_index('##Sample')
    stats_df['frac_circular'] = stats_df['##Circular_Reads']/stats_df['##Mapped_Reads']
    stats_df['frac_mapped'] = stats_df['##Mapped_Reads']/stats_df['##Total_Reads']
    
    # circle counts
    counts = pd.concat(circ_type_count, axis = 1).fillna(0)
    counts.columns = names
    counts = counts.T.reset_index().drop_duplicates('index').set_index('index')

    # BSJ, FSJ and ratio counts
    bsj = pd.concat(bsj_ratio_tbl, axis = 1)
    bsj.columns = names
    bsj=bsj.T.reset_index().drop_duplicates('index').set_index('index').T
    fsj = pd.concat(fsj_ratio_tbl, axis = 1)
    fsj.columns = names
    fsj=fsj.T.reset_index().drop_duplicates('index').set_index('index').T
    junc = pd.concat(junc_ratio_tbl, axis = 1)
    junc.columns = names
    junc=junc.T.reset_index().drop_duplicates('index').set_index('index').T

    mega_anno.to_csv(outdir / 'all_circle_annotation.csv')
    stats_df.to_csv(outdir / 'ciri_stats.csv')
    counts.to_csv(outdir / 'circ_type_counts.csv')
    bsj.to_csv(outdir / 'BSJ_counts.csv')
    fsj.to_csv(outdir / 'FSJ_counts.csv')
    junc.to_csv(outdir / 'junction_ratio.csv')






