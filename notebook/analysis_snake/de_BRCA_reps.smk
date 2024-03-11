# to run differential expression analysis
# from gtf files
"""
snakemake -s /tscc/nfs/home/hsher/projects/circSTAMP_pipe/notebook/analysis_snake/de_BRCA_reps.smk \
--profile /tscc/nfs/home/hsher/projects/circSTAMP_pipe/profiles/tscc2_single
"""

import pandas as pd
from pathlib import Path
import os

workdir: '/tscc/projects/ps-yeolab5/hsher/Tao_circ_BRCA_de'

################# COMPARISON #############

try:
    os.mkdir('error_out_file')
except:
    pass
try:
    os.mkdir('stdout')
except:
    pass
try:
    os.mkdir('output')
except:
    pass

################ MANIFEST ###############
# ciriquant conda env /home/hsher/snakeconda/93c88b426ffde1273db2baa3da0ab4d6

name2gtf = {}

indir1 = '/tscc/projects/ps-yeolab3/tao/CIRI/TBNC_rerun/output' ### CHANGE HERE TO CONTAIN THE CIRI GTF FILES
for f in Path(indir1).glob('*.gtf'):
    name = f.name.replace('_rar11.gtf', ''
        ).replace('_4', '_1' #### WHY DO YOU NAME 231_4 instead of 231_1  !!!! I hate it
        ).replace('_5', '_2'
        ).replace('_6', '_3')
        
    name2gtf[name]=str(f)

name2stringtie = {}
for name in name2gtf:
    f = Path(name2gtf[name])
    folder = f.parent
    name2stringtie[name] = str(folder/'gene'/ (str(f.name.replace('.gtf', '')) +'_out.gtf'))

print(name2gtf)

# make manifest
ctrl = '10A'
others = ['149', '231', '436']
for other in others:
    df = pd.DataFrame([
    [f'{other}_rep1', name2gtf[f'{other}_1'], 'T'],
    [f'{other}_rep2', name2gtf[f'{other}_2'], 'T'],
    [f'{other}_rep3', name2gtf[f'{other}_3'], 'T'],
    [f'{ctrl}_rep1', name2gtf[f'{ctrl}_1'], 'C'],
    [f'{ctrl}_rep2', name2gtf[f'{ctrl}_2'], 'C'],
    [f'{ctrl}_rep3', name2gtf[f'{ctrl}_3'], 'C'],
    ], columns = ['sample_name', 'path', 'group'])

    df.to_csv(f'output/{other}.{ctrl}.lst', sep = '\t', header = False, index = False)

    df_stringtie = pd.DataFrame([
        [f'{other}_rep1', name2stringtie[f'{other}_1']],
        [f'{other}_rep2', name2stringtie[f'{other}_2']],
        [f'{other}_rep3', name2stringtie[f'{other}_3']],
        [f'{ctrl}_rep1', name2stringtie[f'{ctrl}_1']],
        [f'{ctrl}_rep2', name2stringtie[f'{ctrl}_2']],
        [f'{ctrl}_rep3', name2stringtie[f'{ctrl}_3']],
    ], columns = ['sample_name', 'path'])

    df_stringtie.to_csv(f'output/{other}.{ctrl}.samplegene.lst', sep = '\t', header = False, index = False)



comparison = [[other, ctrl] for other in others]


############## RULES #################
rule all:
    input:
        expand("output/unadjusted_comparison/{sample_label_combination}.gtf.tsv", 
            sample_label_combination = [c[0]+'.'+c[1] for c in comparison],
            ),


rule prep_ciriquant_de:
    input:
        lst = ancient('output/{sample_label_control}.{sample_label_case}.lst')
    output:
        lib_info = 'output/{sample_label_control}.{sample_label_case}.csv',
        circ = 'output/{sample_label_control}.{sample_label_case}.circ.csv',
        bsj = 'output/{sample_label_control}.{sample_label_case}.bsj.csv',
        ratio = 'output/{sample_label_control}.{sample_label_case}.ratio.csv',
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "2:00:00",
        cores = "1",
        memory = "160000"
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
        lst = ancient('output/{sample_label_control}.{sample_label_case}.samplegene.lst')
    output:
        genecnt = 'output/{sample_label_control}.{sample_label_case}.genecount_matrix.csv',
        
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "2:00:00",
        cores = "1",
        memory = "160000"
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        python /tscc/nfs/home/hsher/projects/circSTAMP_pipe/scripts/prepDE.py \
            -i {input} \
            -g {output.genecnt}
        """

# this env is hard to install as hell
rule run_ciri_de_compare_unadjusted:
    input:
        lib_info = rules.prep_ciriquant_de.output.lib_info,
        bsj = rules.prep_ciriquant_de.output.bsj,
        gene = rules.prep_stringtie.output.genecnt
    output:
        circ = "output/unadjusted_comparison/{sample_label_control}.{sample_label_case}.gtf.tsv",
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "6:00:00",
        cores = "1",
        memory = "160000"
    container:
        "docker://algaebrown/ciriquant_edger:1.1.2"
    shell:
        """
        /opt/conda/lib/python2.7/site-packages/libs/CIRI_DE.R --lib {input.lib_info} \
            --bsj {input.bsj} \
            --gene {input.gene} \
            --out {output.circ} 
        """


