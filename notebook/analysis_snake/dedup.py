import os
bams = ['/home/hsher/scratch/circ_nextera/output/align/circseq-bm-rar11-nxt.sorted.bam',
'/home/hsher/scratch/circ_nextera/output/align/circseq-bm-rar19-nxt.sorted.bam',
'/home/hsher/scratch/circ_truseq/output/align/circseq-bm-arr.sorted.bam',
'/home/hsher/scratch/circ_truseq/output/align/circseq-bm-rz.sorted.bam',
'/home/hsher/scratch/circ_truseq/output/align/circseq-bm-rar.sorted.bam']

data = dict(zip([os.path.basename(f).split('.')[0] for f in bams], bams))
print(','.join(data.keys()))

rule dedup:
    input:
        bam = lambda wildcards: data[wildcards.sample_label]
    output:
        stat = "{sample_label}.stat",
        rmdup_bam = "{sample_label}.rmdup.bam",
    params:
        error_out_file = "error_files/samtools.{sample_label}",
        run_time = "2:00:00",
        cores = "1",
    container:
        "docker://howardxu520/skipper:samtools_1.17_bedtools_2.31.0"
    shell:
        """

        samtools collate -o {wildcards.sample_label}.collate.bam {input.bam}
        samtools fixmate -m {wildcards.sample_label}.collate.bam {wildcards.sample_label}.fixmate.bam
        samtools sort -o {wildcards.sample_label}.sorted.bam {wildcards.sample_label}.fixmate.bam
        samtools markdup -r -s -f {output.stat} {wildcards.sample_label}.sorted.bam {output.rmdup_bam}
        """
#######################################
import pandas as pd
# snakemake -s dedup.py output/{circseq-bm-rar11-nxt,circseq-bm-rar19-nxt,circseq-bm-arr,circseq-bm-rz,circseq-bm-rar}.gtf --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" -j 5 -kp
config = {'CIRICONFIG': '/home/hsher/scratch/circSTAMP_pipe/ciriconfig_full.yaml',
'LIBRARY_TYPE': 2,
'BWA_INDEX': '/projects/ps-yeolab5/hsher/Tao_circSTAMP/bwa_gencode35',
'HISAT_INDEX': '/projects/ps-yeolab5/hsher/Tao_circSTAMP/hisat_gencode35',
'READ1_ADAPTOR':None,
'READ2_ADAPTOR':None,
'STAR_INDEX':'/projects/ps-yeolab5/hsher/Tao_circSTAMP/star_2_7_gencode35'
}

module ciri:
    snakefile:
        "/home/hsher/scratch/circSTAMP_pipe/SnakeRunCIRI.py"
    config: config

module pre:
    snakefile:
        "/home/hsher/scratch/circSTAMP_pipe/SnakePreprocess.py"
    config: config

workdir: '/home/hsher/scratch/dedup_ciri'
read1 = ['/home/hsher/scratch/circ_nextera/output/fastqs/circseq-bm-rar19-nxt_1.Tr.fq.gz',
'/home/hsher/scratch/circ_nextera/output/fastqs/circseq-bm-rar11-nxt_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq/output/fastqs/circseq-bm-rz_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq/output/fastqs/circseq-bm-rar_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq/output/fastqs/circseq-bm-arr_1.Tr.fq.gz']
read2 = [f.replace('_1', '_2') for f in read1]
names = [os.path.basename(f).split('_1')[0] for f in read1]

fq_data = pd.DataFrame([read1, read2, names], index = ['read1', 'read2', 'name']).T


rule fastp_dedup:
    input:
        fq1 = lambda wildcards: fq_data.loc[fq_data['name']==wildcards.sample_label, 'read1'],
        fq2 =  lambda wildcards: fq_data.loc[fq_data['name']==wildcards.sample_label, 'read2'],
    output:
        fq1 = "dedup_fq/{sample_label}.fq1.gz",
        fq2 = "dedup_fq/{sample_label}.fq2.gz",
        metric = "dedup_metric/{sample_label}.txt"
    params:
        error_out_file = "error_files/fastpdedup.{sample_label}",
        run_time = "5:00:00",
        cores = "1",
    conda:
        '/home/hsher/projects/oligoCLIP/rules/envs/fastp.yaml'
    shell:
        """
        fastp --dedup -A -Q -L -i {input.fq1} -I {input.fq2} -o {output.fq1} -O {output.fq2} > {output.metric}
        """

use rule align_reads from pre as align_start with :
    input:
        fq1="dedup_fq/{sample_label}.fq1.gz",
        fq2="dedup_fq/{sample_label}.fq2.gz",
    output:
        bam = "output/bams/{sample_label}.Aligned.sortedByCoord.out.bam",
        stat = "output/bams/{sample_label}.Log.final.out",
    

use rule run_ciri_RNASE from ciri as dedup_then_ciri with:
    input:
        bwa_index=config['BWA_INDEX']+'.amb',
        hisat_index=config['HISAT_INDEX']+'.1.ht2',
        read1 = "dedup_fq/{sample_label}.fq1.gz",
        read2 = "dedup_fq/{sample_label}.fq2.gz",
        yaml=config['CIRICONFIG']
    output:
        # "output/align/{sample_label}.bam",
        "output/{sample_label}.gtf",