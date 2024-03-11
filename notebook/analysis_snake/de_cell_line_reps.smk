# to run differential expression analysis
# from gtf files
"""
snakemake -s /home/hsher/projects/circSTAMP_pipe/notebook/analysis_snake/de_cell_line_reps.smk --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" -j 15 \
--use-conda --conda-prefix /home/hsher/snakeconda -np
"""

import pandas as pd
from pathlib import Path
import os

workdir: '/projects/ps-yeolab5/hsher/Tao_circ_cellline_de'

################# COMPARISON #############

try:
    os.mkdir('error_out_file')
    os.mkdir('output')
except:
    pass
################ MANIFEST ###############
# ciriquant conda env /home/hsher/snakeconda/93c88b426ffde1273db2baa3da0ab4d6

name2gtf = {'HEK_rep1': '/home/hsher/scratch/circ_nextera/output/circseq-bm-rar11-nxt.gtf',
            'HEK_rep2': '/home/hsher/scratch/circ_nextera_iter2/output/HEK_rar11.gtf'}

indir1 = '/home/hsher/scratch/circ_nextera_iter12/output/'
for f in Path(indir1).glob('*.gtf'):
    name2gtf[f.name.split('_rar11')[0]]=str(f)

name2stringtie = {}
for name in name2gtf:
    f = Path(name2gtf[name])
    folder = f.parent
    print(f.name.replace('.gtf', ''))
    name2stringtie[name] = str(folder/'gene'/ (str(f.name.replace('.gtf', '')) +'_out.gtf'))

# make manifest
for other in ['K562', 'HepG2', 'HeLa']:
    df = pd.DataFrame([
    [f'{other}_rep1', name2gtf[f'{other}_rep1'], 'T'],
    [f'{other}_rep2', name2gtf[f'{other}_rep2'], 'T'],
    [f'HEK_rep1', name2gtf[f'HEK_rep1'], 'C'],
    [f'HEK_rep2', name2gtf[f'HEK_rep2'], 'C'],
    ], columns = ['sample_name', 'path', 'group'])

    df.to_csv(f'output/{other}.HEK.lst', sep = '\t', header = False, index = False)

    df_stringtie = pd.DataFrame([
        [f'{other}_rep1', name2stringtie[f'{other}_rep1']],
        [f'{other}_rep2', name2stringtie[f'{other}_rep2']],
        [f'HEK_rep1', name2stringtie[f'HEK_rep1']],
        [f'HEK_rep2', name2stringtie[f'HEK_rep2']],
    ], columns = ['sample_name', 'path'])

    df_stringtie.to_csv(f'output/{other}.HEK.samplegene.lst', sep = '\t', header = False, index = False)



comparison = [[other, 'HEK'] for other in ['K562', 'HepG2', 'HeLa']]


############## RULES #################
rule all:
    input:
        expand("output/unadjusted_comparison/{sample_label_combination}.gtf.tsv", 
            sample_label_combination = [c[0]+'.'+c[1] for c in comparison],
            ),


rule prep_ciriquant_de:
    input:
        lst = 'output/{sample_label_control}.{sample_label_case}.lst'
    output:
        lib_info = 'output/{sample_label_control}.{sample_label_case}.csv',
        circ = 'output/{sample_label_control}.{sample_label_case}.circ.csv',
        bsj = 'output/{sample_label_control}.{sample_label_case}.bsj.csv',
        ratio = 'output/{sample_label_control}.{sample_label_case}.ratio.csv',
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "2:00:00",
        cores = "1",
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        prep_CIRIquant -i {input} \
            --lib {output.lib_info}\
            --circ {output.circ} \
            --bsj {output.bsj} \
            --ratio {output.ratio}
        """

rule prep_stringtie:
    input:
        lst = 'output/{sample_label_control}.{sample_label_case}.samplegene.lst'
    output:
        genecnt = 'output/{sample_label_control}.{sample_label_case}.genecount_matrix.csv',
        
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "2:00:00",
        cores = "1",
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        python /home/hsher/projects/circSTAMP_pipe/scripts/prepDE.py \
            -i {input} \
            -g {output.genecnt}
        """

rule run_ciri_de_compare_unadjusted:
    input:
        lib_info = rules.prep_ciriquant_de.output.lib_info,
        bsj = rules.prep_ciriquant_de.output.bsj,
        gene = rules.prep_stringtie.output.genecnt
    output:
        "output/unadjusted_comparison/{sample_label_control}.{sample_label_case}.gtf.tsv"
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "6:00:00",
        cores = "1",
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        CIRI_DE_replicate --lib {input.lib_info} \
            --bsj {input.bsj} \
            --gene {input.gene} \
            --out {output}
        """


