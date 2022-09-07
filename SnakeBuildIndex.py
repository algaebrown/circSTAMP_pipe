#snakemake -s SnakeBuildIndex.py -j 3 --configfile config/tao.yaml --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -q home-yeo" --use-conda --conda-frontend conda
import os

BWA_INDEX=config['BWA_INDEX']
HISAT_INDEX=config['HISAT_INDEX']
GENOMEFA=config['GENOMEFA']
GTF=config['GTF']

# rule all:
#     input:
#         BWA_INDEX+'.amb',
#         HISAT_INDEX+'.1.ht2',
#         os.path.join(config['STAR_INDEX'], 'SAindex')
#     output:
#         "BUILD_DONE"
#     params:
#         run_time = "00:03:00",
#         cores = "1",
#     shell:
#         """
#         echo DONE > {output}
#         """

rule get_bwa_index:
    input:
        config['GENOMEFA']
    output:
        BWA_INDEX+'.amb'
    conda:
        "envs/ciriquant.yaml"
    params:
        run_time = "4:00:00",
        cores = "3",
        error_out_file = "error_files/bwabuild.{sample_label}",
    shell:
        """
        bwa index -a bwtsw -p {BWA_INDEX} {input}
        """
rule get_hisat_index:
    input:
        config['GENOMEFA']
    output:
       HISAT_INDEX+'.1.ht2'
    conda:
        "envs/ciriquant.yaml"
    params:
        run_time = "4:00:00",
        error_out_file = "error_files/hisatbuild.{sample_label}",
        cores = "3",
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
        cores = "6",
        run_time = "5:45:00",
        error_out_file = "error_files/starbuild.{sample_label}",
    shell:
        """
        module load star
        STAR \
            --runMode genomeGenerate \
            --runThreadN 8 \
            --genomeDir {params.outdir} \
            --genomeFastaFiles {input.genome_fasta} \
            --sjdbGTFfile {input.gtf} \
        """