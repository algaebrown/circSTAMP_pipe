import pandas as pd
locals().update(config)
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
        tr1=expand("output/fastqs/{sample_label}-trimmed.log", sample_label = sample_labels)
    output:
        tr1='QC/skewer_log.csv',
    params:
        run_time = "00:10:00",
        error_out_file = "error_files/trimstat",
        cores="1",
        memory = 10000,
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/trimming_stat.py "{input.tr1}" {output.tr1}
        """
rule gather_fastqc_report:
    input:
        fq1=expand("fastQC/{sample_label}-trimmed-pair1_fastqc/fastqc_data.txt", sample_label = sample_labels),
        fq2=expand("fastQC/{sample_label}-trimmed-pair2_fastqc/fastqc_data.txt", sample_label = sample_labels),
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
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/fastqc_io.py -i "{input.fq1}" -p {output.passfail1} -b {output.basic1}
        python {SCRIPT_PATH}/fastqc_io.py -i "{input.fq2}" -p {output.passfail2} -b {output.basic2}
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
    shell:
        """
        python {SCRIPT_PATH}/star_mapping_stat_io.py -i "{input}" -o {output}
        """