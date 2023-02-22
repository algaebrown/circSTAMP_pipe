import os
import pandas as pd
# snakemake -s downsample.smk --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" -j 15 --use-conda --conda-prefix /home/hsher/snakeconda
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

workdir: '/home/hsher/scratch/downsample_ciri'
read1 = ['/home/hsher/scratch/circ_nextera/output/fastqs/circseq-bm-rar19-nxt_1.Tr.fq.gz',
'/home/hsher/scratch/circ_nextera/output/fastqs/circseq-bm-rar11-nxt_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq/output/fastqs/circseq-bm-rz_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq/output/fastqs/circseq-bm-rar_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq/output/fastqs/circseq-bm-arr_1.Tr.fq.gz',
'/home/hsher/scratch/circ_nextera_iter2/output/fastqs/HEK_JC_rar11_1.Tr.fq.gz',
'/home/hsher/scratch/circ_nextera_iter2/output/fastqs/HEK_rar11_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq_iter2/output/fastqs/HEK_JC_rar_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq_iter2/output/fastqs/HEK_rar_1.Tr.fq.gz',
'/home/hsher/scratch/circ_truseq_iter2/output/fastqs/EV_rz_1.Tr.fq.gz'
]
read2 = [f.replace('_1', '_2') for f in read1]
names = [os.path.basename(f).split('_1')[0] for f in read1]

fq_data = pd.DataFrame([read1, read2, names], index = ['read1', 'read2', 'name']).T

#to_downsample = ['circseq-bm-rar11-nxt', 'circseq-bm-rz', 'circseq-bm-rar']
#total read=337446074
target_reads = [10**8, 5*10**7, 3*10**7, 10**7, 10**6, 10**4]
rule all:
    input:
        expand("output/{sample_label}.{nread}.gtf", 
        sample_label = names,
        nread = target_reads)
    output:
        "done.txt"
    params:
        error_out_file = "error_files/all",
        run_time = "00:10:00",
        cores = 1,
    shell:
        'echo done > {output}'
    
rule downsample_reads:
    input:
        fq1 = lambda wildcards: fq_data.loc[fq_data['name']==wildcards.sample_label, 'read1'],
        fq2 =  lambda wildcards: fq_data.loc[fq_data['name']==wildcards.sample_label, 'read2'],
    output:
        fq1 = "downsample_fq/{sample_label}.{nread}.fq1.gz",
        fq2 = "downsample_fq/{sample_label}.{nread}.fq2.gz",
    params:
        error_out_file = "error_files/fastpdedup.{sample_label}",
        run_time = "5:00:00",
        cores = "1",
    conda:
        '/home/hsher/projects/oligoCLIP/rules/envs/seqtk.yaml'
    shell:
        """
        seqtk sample -s100 {input.fq1} {wildcards.nread} | gzip > {output.fq1}
        seqtk sample -s100 {input.fq2} {wildcards.nread} | gzip > {output.fq2}
        """

use rule run_ciri_RNASE from ciri as downsample_then_ciri with:
    input:
        bwa_index=config['BWA_INDEX']+'.amb',
        hisat_index=config['HISAT_INDEX']+'.1.ht2',
        read1 = "downsample_fq/{sample_label}.{nread}.fq1.gz",
        read2 = "downsample_fq/{sample_label}.{nread}.fq2.gz",
        yaml=config['CIRICONFIG']
    output:
        "output/{sample_label}.{nread}.gtf",
    params:
        name="{sample_label}.{nread}",
        error_out_file = "error_files/ciri.{sample_label}.{nread}",
        outdir='output/',
        run_time = "16:00:00",
        cores = "4",
        library_type=config['LIBRARY_TYPE'],