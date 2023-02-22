#snakemake -s SnakeMain.py -j 12 --configfile config/tao_nextera3.yaml --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" --use-conda --conda-prefix /home/hsher/snakeconda
#snakemake -s SnakeMain.py -j 12 --configfile config/tao_trueseq3.yaml --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" --use-conda --conda-prefix /home/hsher/snakeconda
#snakemake -s SnakeMain.py -j 12 --configfile config/tao.yaml --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" --use-conda --conda-prefix /home/hsher/snakeconda
import pandas as pd
manifest = pd.read_csv(config['menifest'])

sample_labels = manifest['Sample'].tolist()
rnase_treated_labels = manifest.loc[manifest['Rnase'],'Sample'].tolist()
total_treated_labels = [s.replace('-R', '-A') for s in rnase_treated_labels]

workdir: config['workdir']
try:
    os.mkdir('error_files')
except Exception as e:
    pass

rule all:
    input:
        expand("output/{sample_label}.gtf", sample_label = sample_labels),
        # expand("output/{sample_pair}.gtf", 
        #     sample_pair = [f"{total_sample_label}.{rnase_sample_label}" for total_sample_label, rnase_sample_label in zip(total_treated_labels, rnase_treated_labels)]
        #     ),
        expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai", sample_label = sample_labels),
        expand("fastQC/{sample_label}_2.Tr_fastqc/fastqc_data.txt", sample_label = sample_labels),
        'QC/cutadapt_log.csv',
        'QC/fastQC_basic_summary.r1.csv',
        "QC/genome_mapping_stats.csv"
    output:
        "Main_Done.txt"
    params:
        error_out_file = "error_files/all",
        run_time = "1:00:00",
        cores = "1",
        memory = "20",
        job_name = "all"
    shell:
        """
        echo done > {output}
        """

module build_index:
    snakefile:
        "SnakeBuildIndex.py"
    config: config

module preprocess:
    snakefile:
        "SnakePreprocess.py"
    config: config

module ciri:
    snakefile:
        "SnakeRunCIRI.py"
    config: config

module qc:
    snakefile:
        "SnakeQC.py"
    config: config

# use rule * from build_index as build_*

use rule * from preprocess as pre_*

use rule * from ciri as ciri_*

use rule * from qc as qc_*



