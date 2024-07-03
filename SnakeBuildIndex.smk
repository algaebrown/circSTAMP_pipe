"""
snakemake -s SnakeBuildIndex.smk \
    -j 4 \
    --configfile config/tao_nextera18_placental_EV.yaml \
    --profile profiles/tscc2 \
    -n
"""

import os
locals().update(config)
workdir: config['workdir']

rule all:
    input:
        config['BWA_INDEX']+'.amb',
        config['HISAT_INDEX']+'.1.ht2',
        os.path.join(config['STAR_INDEX'], 'SAindex')
    
rule get_bwa_index:
    input:
        config['GENOMEFA']
    output:
        config['BWA_INDEX']+'.amb'
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    params:
        run_time = "4:00:00",
        cores = "1",
        error_out_file = "error_files/bwabuild",
        out_file = "stdout/bwabuild",
        memory = 60000,
    benchmark: "benchmarks/bwaindex.txt"
    shell:
        """
        bwa index -a bwtsw -p {BWA_INDEX} {input}
        """
rule get_hisat_index:
    input:
        config['GENOMEFA']
    output:
       config['HISAT_INDEX']+'.1.ht2'
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    params:
        run_time = "4:00:00",
        error_out_file = "error_files/hisatbuild",
        out_file = "stdout/hisatbuild",
        memory = 60000,
        cores = "1",
    benchmark: "benchmarks/hisatindex.txt"
    shell:
        """
        hisat2-build {input} {HISAT_INDEX}
        """
        
rule build_star_index:
    input:
        genome_fasta = config['GENOMEFA'],
        gtf=config['GTF'],
    output:
        os.path.join(config['STAR_INDEX'], 'SAindex')
    params:
        outdir=config['STAR_INDEX'],
        cores = "3",
        memory = 480000,
        run_time = "5:45:00",
        error_out_file = "error_files/starbuild",
    benchmark: "benchmarks/starindex.txt"
    container:
       "docker://howardxu520/skipper:star_2.7.10b"
    shell:
        """
        STAR \
            --runMode genomeGenerate \
            --runThreadN 8 \
            --genomeDir {params.outdir} \
            --genomeFastaFiles {input.genome_fasta} \
            --sjdbGTFfile {input.gtf} \
        """