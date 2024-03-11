import pandas as pd
try:
    manifest = pd.read_csv(config['menifest'])
except Exception as e:
    print(e)

# 16 core, 4GB is sufficient to run 4 threads but takes >24hr
# 16 core, 16GB isn't sufficient to run 16 threads (Unable to flush input output)
# 2 core, 16GB mem isn't sufficient to run 8 threads (killed with signal 9)
# 4 core, 16GB mem isn't sufficient to run 8 threads (killed with signal 9)
# but when I do 8 cores, it becomes 64 CPU core available, and it gets into the forward index step, but still, signal 9
# 8 cores, 48GB mem, 8 thread, still signal 9
# 16 cores, 4GB mem, 8 thread, died at a later step, Empty fasta file: '/tscc/lustre/scratch/hsher/circ_nextera_iter13_test/output/circ/APOBEC1only_index.fa'
# no matter how much 2/4 core I requested, it always says 36 CPU cores available, using thread number
# 16 cores, 16GB mem, 8 thread
rule run_ciri_RNASE:
    input:
        bwa_index=ancient(config['BWA_INDEX']+'.amb'),
        hisat_index=ancient(config['HISAT_INDEX']+'.1.ht2'),
        read1 = "output/fastqs/{sample_label}-trimmed-pair1.fastq.gz",
        read2 = "output/fastqs/{sample_label}-trimmed-pair2.fastq.gz",
        yaml=config['CIRICONFIG']
    output:
        "output/align/{sample_label}.sorted.bam",
        "output/align/{sample_label}.sorted.bam.bai",
        "output/{sample_label}.gtf",
        "output/{sample_label}.bed",
        "output/gene/{sample_label}_cov.gtf",
        "output/gene/{sample_label}_genes.list",
        "output/gene/{sample_label}_out.gtf",
        "output/circ/{sample_label}.ciri",
        "output/circ/{sample_label}_index.fa",
        "output/circ/{sample_label}_denovo.sorted.bam",
        "output/circ/{sample_label}_denovo.sorted.bam.bai",
    resources:
        tmpdir = ''
    params:
        name="{sample_label}",
        error_out_file = "error_files/ciri.{sample_label}",
        outdir='output/',
        run_time = "24:00:00",
        cores = 64,
        cpu_threads = 16, # 16 
        memory = 512000, # 512
        library_type=config['LIBRARY_TYPE'],
        #rnase=lambda wildcards: '--RNaseR output/Rnase' if manifest.loc[manifest['Sample']==wildcards.sample_label, 'Rnase'].iloc[0] else ''
        # https://github.com/bioinfo-biols/CIRIquant/issues/4
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    benchmark: "benchmarks/ciriquant.{sample_label}.txt"
    shell:
        """
        CIRIquant \
            -1 {input.read1} \
            -2 {input.read2} \
            --config {input.yaml} \
            --library-type {params.library_type} \
            -o {params.outdir} \
            -p {params.name} \
            -t {params.cpu_threads}
        """

rule index_fa:
    input:
        "output/circ/{sample_label}_index.fa"
    output:
        "output/circ/{sample_label}_index.fa.fai"
    params:
        error_out_file = "error_files/index_fa.{sample_label}",
        run_time = "1:00:00",
        cores = "1",
        memory = 1000
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    benchmark: "benchmarks/index_fa.{sample_label}.txt"
    shell:
        """
        samtools faidx {input}
        """
# RNase correction ruin things
# rule run_ciri_polyA_rnase_coorection:
#     input:
#         bwa_index=config['BWA_INDEX']+'.amb',
#         hisat_index=config['HISAT_INDEX']+'.1.ht2',
#         read1 = "output/fastqs/{total_sample_label}_1.Tr.fq.gz",
#         read2 = "output/fastqs/{total_sample_label}_2.Tr.fq.gz",
#         yaml=config['CIRICONFIG'],
#         rnase_treated_gtf = "output/{rnase_sample_label}.gtf",
#     output:
#         "output/align/{total_sample_label}.{rnase_sample_label}.bam",
#         "output/{total_sample_label}.{rnase_sample_label}.gtf",
#     params:
#         name="{total_sample_label}.{rnase_sample_label}",
#         error_out_file = "error_files/ciri.rnasecorr.{total_sample_label}",
#         outdir='output/',
#         run_time = "24:00:00",
#         cores = "8",
#         library_type=config['LIBRARY_TYPE'],
#         #rnase=lambda wildcards: '--RNaseR output/Rnase' if manifest.loc[manifest['Sample']==wildcards.sample_label, 'Rnase'].iloc[0] else ''
#         # https://github.com/bioinfo-biols/CIRIquant/issues/4
#     conda:
#         "envs/ciriquant.yaml"
#     shell:
#         """
#         CIRIquant -t 6 \
#             -1 {input.read1} \
#             -2 {input.read2} \
#             --config {input.yaml} \
#             --library-type {params.library_type} \
#             -o {params.outdir} \
#             -p {params.name} \
#             --RNaseR {input.rnase_treated_gtf}
#         """


rule gc_cotent:
    input:
        beds = expand("output/{sample_label}.bed", sample_label = manifest['Sample'])
    output:
        megabed = "temp/all_circles.bed",
        nuc = "output/circles_gc_content.nuc"
    params:
        error_out_file = "error_files/gc_content",
        run_time = "2:00:00",
        cores = "1",
        genome=config['GENOMEFA'],
        memory = 1000,
    benchmark: "benchmarks/gc_content.txt"
    container:
        "docker://howardxu520/skipper:bedtools_2.31.0"
    shell:
        """
        cat {input.beds} | sort -k1,1 -k2,2n -k3,3n | uniq > {output.megabed}
        bedtools nuc -s -fi {params.genome} -bed {output.megabed} > {output.nuc}
        """

