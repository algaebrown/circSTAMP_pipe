#snakemake -s SnakeQC.py -j 6 --configfile config/tao.yaml --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" --use-conda
import pandas as pd

manifest = pd.read_csv(config['menifest'])
sample_labels = manifest['Sample'].tolist()

rule all:
    input:
        'QC/cutadapt_log.csv',
        'QC/fastQC_basic_summary.r1.csv',
        "QC/genome_mapping_stats.csv"
    output:
        "QC_Done.txt"
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

rule gather_trimming_stat:
    input:
        tr1=expand("output/fastqs/{sample_label}_1.Tr.metrics", sample_label = sample_labels)
    output:
        tr1='QC/cutadapt_log.csv',
    params:
        run_time = "00:10:00",
        cores="1",
        files1=','.join(expand("output/fastqs/{sample_label}_1.Tr.metrics", sample_label = sample_labels)),
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python /home/hsher/projects/QC_tools/trimming_stat.py {params.files1} {output.tr1}
        """
rule gather_fastqc_report:
    input:
        fq1=expand("fastQC/{sample_label}_1.Tr_fastqc/fastqc_data.txt", sample_label = sample_labels),
        fq2=expand("fastQC/{sample_label}_2.Tr_fastqc/fastqc_data.txt", sample_label = sample_labels),
    output:
        basic1='QC/fastQC_basic_summary.r1.csv',
        passfail1='QC/fastQC_passfail.r1.csv',
        basic2='QC/fastQC_basic_summary.r2.csv',
        passfail2='QC/fastQC_passfail.r2.csv'
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
        files1 = ','.join(expand("fastQC/{sample_label}_1.Tr_fastqc/fastqc_data.txt", 
            sample_label = sample_labels)),
        files2 = ','.join(expand("fastQC/{sample_label}_2.Tr_fastqc/fastqc_data.txt",
            sample_label = sample_labels
        ))
    shell:
        """
        python /home/hsher/projects/QC_tools/fastqc_io.py -i {params.files1} -p {output.passfail1} -b {output.basic1}
        python /home/hsher/projects/QC_tools/fastqc_io.py -i {params.files2} -p {output.passfail2} -b {output.basic2}
        """
rule gather_GENOME_mapstat:
    input:
        #find_all_files("{libname}/bams/genome/{sample_label}.genome-mapped.Log.final.out", libs)
        expand("output/bams/{sample_label}.Log.final.out", sample_label = sample_labels)
    output:
        "QC/genome_mapping_stats.csv"
    conda:
        "envs/metadensity.yaml"
    params:
        error_out_file = "error_files/mapstat",
        run_time = "00:40:00",
        cores = "1",
        memory = "10000",
        job_name = "gather_stat",
        files = ','.join(expand("output/bams/{sample_label}.Log.final.out", sample_label = sample_labels))
    shell:
        """
        python /home/hsher/projects/QC_tools/star_mapping_stat_io.py -i {params.files} -o {output}
        """