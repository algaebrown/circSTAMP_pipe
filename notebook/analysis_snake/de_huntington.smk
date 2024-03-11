# to run differential expression analysis
# from gtf files
"""
snakemake -s de_huntington.smk \
--cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" \
-j 15 --use-conda --conda-prefix /home/hsher/snakeconda -n
"""
import pandas as pd
from pathlib import Path
import os

workdir: '/projects/ps-yeolab5/hsher/Tao_circ_huntington'

################# COMPARISON #############
# control, case
indir = Path('/projects/ps-yeolab3/tao/CIRI/organoids/output/')
comparison = [
    ['CVB_STO', 'HT_109_STO_2'],
    ['CVB_STO', 'HT_109_STO_1'],
]

try:
    os.mkdir('error_out_file')
except:
    pass


############## RULES #################
rule all:
    input:
        expand("output/unadjusted_comparison/{sample_label_combination}.gtf.tsv", 
            sample_label_combination = [c[0]+'.'+c[1] for c in comparison],
            ),

rule run_ciri_de_compare_unadjusted:
    input:
        gtf_control = lambda wildcards: indir / f'{wildcards.sample_label_control}.gtf',
        gtf_case = lambda wildcards: indir / f'{wildcards.sample_label_case}.gtf'
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
        CIRI_DE -n {input.gtf_control} \
        -c {input.gtf_case} \
        -o {output}
        """
