import pandas as pd
import os
workdir: config['workdir']

locals().update(config)
manifest = pd.read_csv(config['menifest'])

# check file existence
assert manifest['fastq1'].apply(os.path.isfile).all()
assert manifest['fastq2'].apply(os.path.isfile).all()
assert os.path.isfile(config['CIRICONFIG'])

# check config is correct
import yaml
with open(config['CIRICONFIG'], 'r') as file:
    ciriconfig = yaml.safe_load(file)

# check if ciriconfig as the right key and path exists
for key in ['bwa', 'hisat2', 'stringtie', 'samtools']:
    assert key in ciriconfig['tools']
    assert os.path.isfile(ciriconfig['tools'][key])

assert ciriconfig['reference']['fasta']==config['GENOMEFA']
assert ciriconfig['reference']['gtf']==config['GTF']


sample_labels = manifest['Sample'].tolist()
config['sample_labels']=sample_labels


# handle exception
for key in ['STAMP', 'STAMP_control']:
    try:
        if config[key] is None:
            print(f'not a STAMP experiment, {key} is empty')
            config[key] = []
            # placeholder
            config['REF_fwd'] = 'C'
            config['ALT_fwd'] = 'T'
    except KeyError:
        config[key] = []
        # placeholder
        config['REF_fwd'] = 'C'
        config['ALT_fwd'] = 'T'

for key in ['external_stamp_control', 'fit_overdispersion_from']:
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
revcomp = {'C': 'G',
           'G': 'C',
           'A': 'T',
           'T': 'A'}

config['REF_rev'] = revcomp[config['REF_fwd']]
config['ALT_rev'] = revcomp[config['ALT_fwd']]
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
        # preprocessing
        expand("output/fastqs/{sample_label}-trimmed.log", sample_label = sample_labels),
        expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai", sample_label = sample_labels),
        expand("fastQC/{sample_label}-trimmed-pair1_fastqc/fastqc_data.txt", sample_label = sample_labels),
        'QC/skewer_log.csv',
        'QC/fastQC_basic_summary.r1.csv',
        "QC/genome_mapping_stats.csv",
        # CIRI
        expand("output/{sample_label}.gtf", sample_label = sample_labels),
        'output/circle_summary/all_circle_annotation.csv',
        'output/circle_summary/ciri_stats.csv',
        'output/circle_summary/circ_type_counts.csv',
        'output/circle_summary/BSJ_counts.csv',
        'output/circle_summary/FSJ_counts.csv',
        'output/circle_summary/junction_ratio.csv',
        
        # RIP (circRIP is not good)
        # expand("output/circRIP/{comparisons}", comparisons = comparisons) if CIRCRIP_PATH else [],
        # expand("output/circRIP/homer/{comparisons}.homer", comparisons = comparisons) if CIRCRIP_PATH else [],
        expand("output/RIP/{comparisons}.csv", comparisons = comparisons2),
        # STAMP
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