"""
snakemake -s SnakeMain.smk \
    -j 12 \
    --configfile config/tao_nextera16.yaml \
    --profile profiles/tscc2 \
    --until ciri_run_ciri_RNASE \
    -n
"""
import pandas as pd
import os
locals().update(config)

manifest = pd.read_csv(config['menifest'])

# check file existence
assert manifest['fastq1'].apply(os.path.isfile).all()
assert manifest['fastq2'].apply(os.path.isfile).all()


sample_labels = manifest['Sample'].tolist()
config['sample_labels']=sample_labels
workdir: config['workdir']

# handle exception
for key in ['STAMP', 'STAMP_control', 'fit_overdispersion_from']:
    try:
        if config[key] is None:
            config[key] = []
    except KeyError:
        config[key] = []

for key in ['external_stamp_control']:
    try:
        if config[key] is None:
            config[key] = {}
    except KeyError:
        config[key] = {}


# make error file directories
try:
    os.mkdir('error_files')
except Exception as e:
    pass

try:
    os.mkdir('stdout')
except Exception as e:
    pass

# edit calling
NREAD=5
REF_fwd='C' # C to T
REF_rev='G' # G to A

ALT_fwd='T'
ALT_rev='A'

locals().update(config)

# RIP comparisons
if config['RIP_comparison'] is not None:
    comparisons = [comp['ip']+'_vs_'+comp['in'] for comp in config['RIP_comparison'].values()]
    comparisons2 = [comp['ip']+'.'+comp['in'] for comp in config['RIP_comparison'].values()]
else:
    comparisons = []
    comparisons2 = []


rule all:
    input:
        expand("output/fastqs/{sample_label}-trimmed.log", sample_label = sample_labels),
        expand("output/{sample_label}.gtf", sample_label = sample_labels),
        expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai", sample_label = sample_labels),
        expand("fastQC/{sample_label}-trimmed-pair1_fastqc/fastqc_data.txt", sample_label = sample_labels),
        'QC/skewer_log.csv',
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
        ctrl_sample_label = config['STAMP_control']+list(config['external_stamp_control'].keys())
        ),
        expand("output/edits/homer/{stamp_sample_label}.{ctrl_sample_label}/homerResults.html",
        stamp_sample_label = config['STAMP'],
        ctrl_sample_label = config['STAMP_control']+list(config['external_stamp_control'].keys())
        ),

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

use rule * from preprocess as pre_*

use rule * from ciri as ciri_*

use rule * from qc as qc_*

use rule * from count_edit as edit_*

use rule * from rip as rip_*