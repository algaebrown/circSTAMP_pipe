import pandas as pd
try:
    manifest = pd.read_csv(config['menifest'])
except Exception as e:
    print(e)

rule run_ciri_RNASE:
    input:
        bwa_index=config['BWA_INDEX']+'.amb',
        hisat_index=config['HISAT_INDEX']+'.1.ht2',
        read1 = "output/fastqs/{sample_label}_1.Tr.fq.gz",
        read2 = "output/fastqs/{sample_label}_2.Tr.fq.gz",
        yaml=config['CIRICONFIG']
    output:
        # "output/align/{sample_label}.bam",
        "output/{sample_label}.gtf",
        "output/{sample_label}.bed",
        # "output/{sample_label}.bed",
        # "output/gene/{sample_label}_cov.gtf",
        # "output/gene/{sample_label}_genes.list",
        # "output/gene/{sample_label}_out.gtf",
        # "output/align/{sample_label}.sorted.bam",
        # "output/align/{sample_label}.sorted.bam.bai",
        # "output/circ/{sample_label}.ciri",
        #expand("output/circ/{sample_label}_index.{number}.ht2", number = list(range(20)), allow_missing=True),
        "output/circ/{sample_label}_index.fa",
        "output/circ/{sample_label}_denovo.sorted.bam",
        "output/circ/{sample_label}_denovo.sorted.bam.bai",
    params:
        name="{sample_label}",
        error_out_file = "error_files/ciri.{sample_label}",
        outdir='output/',
        run_time = "60:00:00",
        cores = "16",
        library_type=config['LIBRARY_TYPE'],
        #rnase=lambda wildcards: '--RNaseR output/Rnase' if manifest.loc[manifest['Sample']==wildcards.sample_label, 'Rnase'].iloc[0] else ''
        # https://github.com/bioinfo-biols/CIRIquant/issues/4
    conda:
        "envs/ciriquant.yaml"
    shell:
        """
        # remove incomplete files if exist. don't work, make things fail
        #rm output/gene/{wildcards.sample_label}*
        #rm output/align/{wildcards.sample_label}*
        #rm output/circ/{wildcards.sample_label}*

        CIRIquant \
            -1 {input.read1} \
            -2 {input.read2} \
            --config {input.yaml} \
            --library-type {params.library_type} \
            -o {params.outdir} \
            -p {params.name} \
            -t 4
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
    shell:
        """
        module load bedtools
        cat {input.beds} | sort -k1,1 -k2,2n -k3,3n | uniq > {output.megabed}
        bedtools nuc -s -fi {params.genome} -bed {output.megabed} > {output.nuc}
        """

