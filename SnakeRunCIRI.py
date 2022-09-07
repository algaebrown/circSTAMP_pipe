import pandas as pd
manifest = pd.read_csv(config['menifest'])

rule run_ciri_RNASE:
    input:
        bwa_index=config['BWA_INDEX']+'.amb',
        hisat_index=config['HISAT_INDEX']+'.1.ht2',
        read1 = "output/fastqs/{sample_label}_1.Tr.fq.gz",
        read2 = "output/fastqs/{sample_label}_2.Tr.fq.gz",
        yaml=config['CIRICONFIG']
    output:
        "output/align/{sample_label}.bam",
        "output/{sample_label}.gtf",
    params:
        name="{sample_label}",
        error_out_file = "error_files/ciri.{sample_label}",
        outdir='output/',
        run_time = "24:00:00",
        cores = "8",
        library_type=config['LIBRARY_TYPE'],
        #rnase=lambda wildcards: '--RNaseR output/Rnase' if manifest.loc[manifest['Sample']==wildcards.sample_label, 'Rnase'].iloc[0] else ''
        # https://github.com/bioinfo-biols/CIRIquant/issues/4
    conda:
        "envs/ciriquant.yaml"
    shell:
        """
        CIRIquant -t 4 \
            -1 {input.read1} \
            -2 {input.read2} \
            --config {input.yaml} \
            --library-type {params.library_type} \
            -o {params.outdir} \
            -p {params.name} \
            -t 12
        """

rule run_ciri_polyA_rnase_coorection:
    input:
        bwa_index=config['BWA_INDEX']+'.amb',
        hisat_index=config['HISAT_INDEX']+'.1.ht2',
        read1 = "output/fastqs/{total_sample_label}_1.Tr.fq.gz",
        read2 = "output/fastqs/{total_sample_label}_2.Tr.fq.gz",
        yaml=config['CIRICONFIG'],
        rnase_treated_gtf = "output/{rnase_sample_label}.gtf",
    output:
        "output/align/{total_sample_label}.{rnase_sample_label}.bam",
        "output/{total_sample_label}.{rnase_sample_label}.gtf",
    params:
        name="{total_sample_label}.{rnase_sample_label}",
        error_out_file = "error_files/ciri.rnasecorr.{total_sample_label}",
        outdir='output/',
        run_time = "24:00:00",
        cores = "8",
        library_type=config['LIBRARY_TYPE'],
        #rnase=lambda wildcards: '--RNaseR output/Rnase' if manifest.loc[manifest['Sample']==wildcards.sample_label, 'Rnase'].iloc[0] else ''
        # https://github.com/bioinfo-biols/CIRIquant/issues/4
    conda:
        "envs/ciriquant.yaml"
    shell:
        """
        CIRIquant -t 4 \
            -1 {input.read1} \
            -2 {input.read2} \
            --config {input.yaml} \
            --library-type {params.library_type} \
            -o {params.outdir} \
            -p {params.name} \
            -t 12 \
            --RNaseR {input.rnase_treated_gtf}
        """

