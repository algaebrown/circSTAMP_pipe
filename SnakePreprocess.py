import pandas as pd
import glob
try:
    manifest = pd.read_csv(config['menifest'])
except:
    pass

rule cutadapt:
    input:
        fq1=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["fastq1"].values[0]),
        fq2=lambda wildcards: glob.glob(manifest.loc[manifest.Sample == wildcards.sample_label]["fastq2"].values[0]),
    output:
        fq1=temp("output/fastqs/{sample_label}_1.Tr.fq.gz"),
        fq2=temp("output/fastqs/{sample_label}_2.Tr.fq.gz"),
        metrics = "output/fastqs/{sample_label}_1.Tr.metrics",
    params:
        error_out_file = "error_files/cutadapt.{sample_label}",
        run_time = "6:45:00",
        cores = "4",
        memory = "10000",
        job_name = "cutadapt",
        adaptor1 = config['READ1_ADAPTOR'],
        adaptor2 = config['READ2_ADAPTOR'],
    benchmark: "benchmarks/cutadapt/{sample_label}.extract.txt"
    conda:
        "envs/cutadapt.yaml"
    shell:
        """
        cutadapt -A file:{params.adaptor2} \
            -a file:{params.adaptor1} \
            -o {output.fq1} -p {output.fq2}  \
            --pair-filter=any \
            --times=2 \
            -m 20 \
            {input.fq1} {input.fq2} > {output.metrics}
        """

rule fastQC:
    input:
        fq1="output/fastqs/{sample_label}_1.Tr.fq.gz",
        fq2="output/fastqs/{sample_label}_2.Tr.fq.gz"
    output:
        fq1_zip="fastQC/{sample_label}_1.Tr_fastqc/fastqc_data.txt",
        fq2_zip="fastQC/{sample_label}_2.Tr_fastqc/fastqc_data.txt",
    params:
        outdir="fastQC/",
        error_out_file = "error_files/fastqc.{sample_label}",
        thread=6,
        run_time = "1:45:00",
        cores = "2",
        memory = "10000",
        job_name = "fastQC",
    shell:
        """
        module load fastqc
        fastqc {input.fq1} --extract --outdir {params.outdir} -t {params.thread}
        fastqc {input.fq2} --extract --outdir {params.outdir} -t {params.thread}
        """
rule align_reads:
    input:
        fq1="output/fastqs/{sample_label}_1.Tr.fq.gz",
        fq2="output/fastqs/{sample_label}_2.Tr.fq.gz"
    output:
        bam = "output/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
        stat = "output/bams/{sample_label}.Log.final.out",
    params:
        error_out_file = "error_files/star.{sample_label}",
        run_time = "03:40:00",
        cores = "4",
        memory = "10000",
        job_name = "align_reads",
        thread = "6",
        star_db=config['STAR_INDEX'],
        prefix="output/bams/{sample_label}."
    benchmark: "benchmarks/align/{sample_label}.align_reads.txt"
    shell:
        """
        module load star
        STAR --genomeDir {params.star_db} \
            --runThreadN {params.thread} \
            --readFilesIn {input.fq1} {input.fq2} \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outSAMattributes All \
            --outFilterMultimapNmax 1 \
            --readFilesCommand zcat
        """
        

rule index_bams:
    input:
        bam = "output/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
    output:
        #bam_sorted="output/bams/{sample_label}.Aligned.out.sorted.bam",
        bai = "output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai",
    params:
        error_out_file = "error_files/indexbam.{sample_label}",
        run_time = "01:40:00",
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    benchmark: "benchmarks/align/{sample_label}.index_bam.txt"
    shell:
        "module load samtools;"
        "samtools index {input.bam};"

rule quantify_linear:
    input:
        bams = expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam",sample_label = manifest['Sample']),
        bais = expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai",sample_label = manifest['Sample'])
    output:
        "output/featureCount_output.tsv"
    params:
        error_out_file = "error_files/indexbam",
        run_time = "01:40:00",
        gtf = config['GTF'],
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    shell:
        """
        module load subreadfeaturecounts
        featureCounts -p \
            -a {params.gtf} \
            -o {output} {input.bams} \
            -g gene_name
        """

rule quantify_linear_mRNA:
    input:
        bams = expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam",sample_label = manifest['Sample']),
        bais = expand("output/bams/{sample_label}.Aligned.sortedByCoord.out.bam.bai",sample_label = manifest['Sample'])
    output:
        "output/featureCount_output_mRNA.tsv"
    params:
        error_out_file = "error_files/indexbam",
        run_time = "01:40:00",
        gtf = config['GTF'],
        cores = "4",
        memory = "1000",
        job_name = "index_bam"
    shell:
        """
        module load subreadfeaturecounts
        featureCounts -p \
            -a {params.gtf} \
            -o {output} {input.bams} \
            -g gene_name \
            -t exon \
            --splitOnly
        """