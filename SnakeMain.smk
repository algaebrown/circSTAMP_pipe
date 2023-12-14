"""
snakemake -s SnakeMain.smk \
    -j 12 \
    /tscc/nfs/home/hsher/projects/circSTAMP_pipe/config/tao_nextera13.tscc2.yaml \
    --cluster "sbatch -t {params.run_time} -n {params.cores} -e {params.error_out_file} -q hotel" \
    --use-conda \
    --conda-prefix /tscc/projects/ps-yeolab5/hsher/snakeconda \
    --conda-create-envs-only

snakemake -s SnakeMain.smk \
    -j 12 \
    --configfile config/tao_trueseq2.yaml \
    --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" \
    --use-conda \
    --conda-prefix /home/hsher/snakeconda \
    -n
"""
import pandas as pd
import os
print(config.keys())
manifest = pd.read_csv(config['menifest'])

# check file existence
assert manifest['fastq1'].apply(os.path.isfile).all()
assert manifest['fastq2'].apply(os.path.isfile).all()


sample_labels = manifest['Sample'].tolist()
config['sample_labels']=sample_labels
workdir: config['workdir']
# rnase_treated_labels = manifest.loc[manifest['Rnase'],'Sample'].tolist()
# total_treated_labels = [s.replace('-R', '-A') for s in rnase_treated_labels]
try:
    if config['fit_overdispersion_from'] is None:
        config['fit_overdispersion_from'] = []
except KeyError:
    config['fit_overdispersion_from'] = []

if 'STAMP' not in config:
    config['STAMP']=[]

if 'STAMP_control' not in config:
    config['STAMP_control'] = []
if config['STAMP'] is None:
    config['STAMP'] = []

if config['STAMP_control'] is None:
    config['STAMP_control'] = []


try:
    os.mkdir('error_files')
except Exception as e:
    pass

# edit calling
NREAD=5
REF_fwd='C' # C to T
REF_rev='G' # G to A

ALT_fwd='T'
ALT_rev='A'

SCRIPT_PATH = config['SCRIPT_PATH']
CIRCRIP_PATH = config['CIRCRIP_PATH']
R_EXE = config['R_EXE']

# RIP comparisons
if config['RIP_comparison'] is not None:
    comparisons = [comp['ip']+'_vs_'+comp['in'] for comp in config['RIP_comparison'].values()]
    comparisons2 = [comp['ip']+'.'+comp['in'] for comp in config['RIP_comparison'].values()]
else:
    comparisons = []
    comparisons2 = []

print(comparisons)

rule all:
    input:
        expand("output/{sample_label}.gtf", sample_label = sample_labels),
        # expand("output/{sample_pair}.gtf", 
        #     sample_pair = [f"{total_sample_label}.{rnase_sample_label}" for total_sample_label, rnase_sample_label in zip(total_treated_labels, rnase_treated_labels)]
        #     ),
        expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai", sample_label = sample_labels),
        expand("fastQC/{sample_label}_2.Tr_fastqc/fastqc_data.txt", sample_label = sample_labels),
        #expand("output/edits/{sample_label}.dp4.{strand}.circ_level.tsv", sample_label = sample_labels, strand = ['pos', 'neg']),
        'QC/cutadapt_log.csv',
        'QC/fastQC_basic_summary.r1.csv',
        "QC/genome_mapping_stats.csv",
        "output/featureCount_output.tsv",
        "output/featureCount_output_mRNA.tsv",
        "output/circles_gc_content.nuc",
        expand("output/circRIP/{comparisons}", comparisons = comparisons),
        expand("output/circRIP/homer/{comparisons}.homer", comparisons = comparisons),
        expand("output/RIP/{comparisons}.csv", comparisons = comparisons2),
        expand("output/edits/{sample_label}.final_edit_score.fa"
        , sample_label = sample_labels),
        expand("output/edits/homer/{stamp_sample_label}.{ctrl_sample_label}.homer",
        stamp_sample_label = config['STAMP'],
        ctrl_sample_label = config['STAMP_control']
        ),
        expand("output/edits/{sample_label}.expand_edit_bed.fa"
        , sample_label = sample_labels)
    

module build_index:
    snakefile:
        "SnakeBuildIndex.smk"
    config: config

module preprocess:
    snakefile:
        "SnakePreprocess.smk"
    config: config

module ciri:
    snakefile:
        "SnakeRunCIRI.smk"
    config: config

module qc:
    snakefile:
        "SnakeQC.smk"
    config: config

module count_edit:
    snakefile:
        "SnakeCircEdit.smk"
    config: config

module rip:
    snakefile:
        "SnakeRIP.smk"
    config:
        config

# use rule * from build_index as build_*

use rule * from preprocess as pre_*

use rule * from ciri as ciri_*

use rule * from qc as qc_*

use rule * from count_edit as edit_*

use rule * from rip as rip_*