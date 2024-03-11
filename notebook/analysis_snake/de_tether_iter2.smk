# to run differential expression analysis
# from gtf files
"""
snakemake -j 16 \
    -s de_tether_iter2.smk \
    --cluster "sbatch -t {params.run_time} -e {params.error_out_file} -o stdout/{rule} -p condo -q condo -A csd792 --tasks-per-node {params.cores} --job-name {rule} --mem {params.memory}" \
    --use-conda \
    --conda-prefix /tscc/nfs/home/hsher/snakeconda \
    -n
"""

import pandas as pd
from pathlib import Path
import os

workdir: '/tscc/projects/ps-yeolab5/hsher/Tao_circ_tether_de_iter2'

try:
    os.mkdir('error_files')
    os.mkdir('stdout')
except:
    pass
# libname,RNase+,RNase-,sample_label
indir = Path('/tscc/nfs/home/hsher/scratch/circ_nextera_iter9/output/')
gtfs = list(indir.glob('*.gtf'))
sample_labels = [f.name.replace('.gtf', '') for f in gtfs]

manifest = pd.DataFrame([gtfs, sample_labels], index = ['RNase+', 'sample_label']).T

ov_trials = [s for s in manifest['sample_label'] if not s.startswith('si') and s != 'EV']

# comapre siXX to siNT, overexpress to EV(empty vector)
comparison = [['EV', case] for case in ov_trials]
print('comparisons:',comparison)

############## RULES #################
rule all:
    input:
        expand("output/unadjusted_comparison/{sample_label_combination}.gtf.tsv", 
            sample_label_combination = [c[0]+'.'+c[1] for c in comparison],
            ),
        

rule run_ciri_de_compare_unadjusted:
    input:
        gtf_control = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label_control, 'RNase+'],
        gtf_case = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label_case, 'RNase+'],
    output:
        "output/unadjusted_comparison/{sample_label_control}.{sample_label_case}.gtf.tsv"
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "8:00:00",
        cores = "1",
        memory = 320000,
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        CIRI_DE -n {input.gtf_control} \
        -c {input.gtf_case} \
        -o {output} 
        """

