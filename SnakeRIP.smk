# to run circRIP pipeline
import pandas as pd
locals().update(config)
manifest = pd.read_csv(config['menifest'])
sample_labels = manifest['Sample'].tolist()


### CIRC-RIP: it kinda doesn't work. Rules are set for benchmarking ###
''' convert gtf file to .circ file format as required by circRIP '''
rule gtf_to_circ:
    input:
        gtf = "output/{sample_label}.gtf",
    output:
        circ = "output/{sample_label}.circ",
    params:
        error_out_file = "error_files/convert_file.{sample_label}",
        run_time = "20:00",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/gtf2circ.py {input.gtf} {output.circ}
        """

rule run_circ_RIP:
    input:
        ip_circ="output/{ip_sample_label}.circ",
        in_circ="output/{in_sample_label}.circ",
        gtf=config['GTF'],
        ip_bam="output/bams/{ip_sample_label}.Aligned.sortedByCoord.out.bam",
        ip_bai="output/bams/{ip_sample_label}.Aligned.sortedByCoord.out.bam.bai",
        in_bam="output/bams/{in_sample_label}.Aligned.sortedByCoord.out.bam",
        in_bai="output/bams/{in_sample_label}.Aligned.sortedByCoord.out.bam.bai",
        genome=config['GENOMEFA']
    output:
        "output/circRIP/{ip_sample_label}_vs_{in_sample_label}"
    params:
        error_out_file = "error_files/circRIP.{ip_sample_label}.{in_sample_label}",
        run_time = "12:00:00",
        cores = "2",
    conda:
        "envs/circRIP.yaml"
    shell:
        """
        python {CIRCRIP_PATH}/circRIP.py EnrichedcircRNA \
            -ip_circ {input.ip_circ} \
            -input_circ {input.in_circ} \
            -gtf {input.gtf} \
            -ip_bam {input.ip_bam} \
            -input_bam {input.in_bam} \
            -prefix {output} \
            -G {input.genome}
        """

rule prepare_motif_analysis:
    input:
        enrich_result = "output/circRIP/{ip_sample_label}_vs_{in_sample_label}",
        fa = "output/circ/{ip_sample_label}_index.fa"
    output:
        foreground = "output/circRIP/homer/{ip_sample_label}_vs_{in_sample_label}.enriched.fa",
        background = "output/circRIP/homer/{ip_sample_label}_vs_{in_sample_label}.background.fa",
    conda:
        "envs/metadensity.yaml"
    params:
        error_out_file = "error_files/circRIP.{ip_sample_label}.{in_sample_label}",
        run_time = "1:00:00",
        cores = "1",
        prefix = "output/circRIP/homer/{ip_sample_label}_vs_{in_sample_label}"
    shell:
        """
        python {SCRIPT_PATH}/prepare_homer.py \
            {input.enrich_result} \
            {input.fa} \
            {params.prefix}
        """

rule homer:
    input:
        foreground = "output/circRIP/homer/{ip_sample_label}_vs_{in_sample_label}.enriched.fa",
        background = "output/circRIP/homer/{ip_sample_label}_vs_{in_sample_label}.background.fa",
    output:
        "output/circRIP/homer/{ip_sample_label}_vs_{in_sample_label}.homer"
    params:
        error_out_file = "error_files/circRIP_homer.{ip_sample_label}.{in_sample_label}",
        run_time = "8:00:00",
        cores = "1",
    container:
        "docker://howardxu520/skipper:Homer_4.11"
    shell:
        """
        homer2 denovo -i {input.foreground} -b {input.background} -strand + -o {output}
        """
### My way of calling enrichment ###
rule get_counts:
    input:
        expand(expand("output/{sample_label}.gtf", sample_label = sample_labels))
    output:
        "output/count_table.tsv"
    params:
        error_out_file = "error_files/make_count_table",
        run_time = "0:30:00",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/make_count_table.py \
            output \
            {output}
        """

rule fit_overdispersion:
    input:
        "output/count_table.tsv"
    output:
        "output/RIP/coef.tsv"
    params:
        error_out_file = "error_files/fit_overdispersion",
        run_time = "0:30:00",
        cores = "1",
        reps = ','.join(config['fit_overdispersion_from'])
    container:
        "docker://howardxu520/skipper:R_4.1.3_1"
    shell:
        """
        Rscript --vanilla {SCRIPT_PATH}/fit_overdispersion.R \
            {input} \
            {params.reps} \
            {output}
        """

rule call_enriched_RIP:
    input:
        counts = "output/count_table.tsv",
        coef = "output/RIP/coef.tsv"
    output:
        "output/RIP/{ip_sample_label}.{in_sample_label}.csv"
    params:
        error_out_file = "error_files/call_RIP.{ip_sample_label}.{in_sample_label}",
        run_time = "0:30:00",
        cores = "1"
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/call_enriched_RIP.py \
            {input.counts} \
            {input.coef} \
            {wildcards.ip_sample_label} \
            {wildcards.in_sample_label} \
            {output}
        """