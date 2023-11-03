"""
snakemake -s SnakeCircEdit.smk -j 12 \
    --configfile config/tao_nextera13.yaml \
    --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" \
    --use-conda --conda-prefix /home/hsher/snakeconda -n
"""
from pathlib import Path
import pandas as pd
SCRIPT_PATH=config['SCRIPT_PATH']
# manifest = pd.read_csv(config['menifest'])
# print(manifest)

# BSJ_threshold=500

# sample_labels = manifest['Sample'].tolist()
# workdir: config['workdir']

# def get_all_circ_id(sample_label):
#     df = pd.read_csv(f'output/{sample_label}.gtf', names = ['chr', 'source', 'type','start','end',
#     'score', 'strand', '.', 'annotation'], sep = '\t', comment = '#'
#     )
#     print(df.head())
#     df['name']=df['annotation'].apply(lambda string: string.split('circ_id ')[1].split(';')[0].replace('\"', ''))
#     df['bsj']=df['annotation'].apply(lambda string: float(string.split('bsj ')[1].split(';')[0]))
#     return df.loc[df['bsj']>BSJ_threshold]['name'].tolist()

# def get_all_outputs_for_1_sample(sample_label):
#     circids = get_all_circ_id(sample_label)
#     outputs = [f"output/edits/{sample_label}.dp4.{strand}.snp_filtered/{circ_id}.tsv"
#     for strand in ['pos', 'neg'] for circ_id in circids
#     ]
#     return outputs

# def get_all_outputs(samples):
#     outputs = []
#     for s in samples:
#         outputs += get_all_outputs_for_1_sample(s)
#     return outputs

# print([len(get_all_circ_id(s)) for s in sample_labels])

# rule all:
#     input: 
#         get_all_outputs(sample_labels)

rule pileup:
    input:
        circ_ref = "output/circ/{sample_label}_index.fa",
        bam = "output/circ/{sample_label}_denovo.sorted.bam",
        bai = "output/circ/{sample_label}_denovo.sorted.bam.bai",
    output:
        vcf="output/edits/{sample_label}.dp4.vcf",
    params:
        error_out_file = "error_files/vcf.{sample_label}",
        run_time = "12:00:00",
        cores = "3"
    shell:
        """
        module load bcftools 
        
        bcftools mpileup -f {input.circ_ref} {input.bam} \
            --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR > {output.vcf}
        """


REF_fwd='C' # C to T
REF_rev='G' # G to A

ALT_fwd='T'
ALT_rev='A'

rule filter_variants:
    # need to consider strand specific thing
    input:
        vcf="output/edits/{sample_label}.dp4.vcf"
    output:
        pos_vcf = "output/edits/{sample_label}.dp4.pos.vcf",
        neg_vcf = "output/edits/{sample_label}.dp4.neg.vcf",
    params:
        error_out_file = "error_files/filter_vcf.{sample_label}",
        run_time = "12:00:00",
        cores = "3"
    shell:
        """
        module load bcftools 
        bcftools filter -i 'REF="{REF_fwd}" & TYPE!="indel"' --output {output.pos_vcf} {input.vcf}
        bcftools filter -i 'REF="{REF_rev}" & TYPE!="indel"' --output {output.neg_vcf} {input.vcf}
        """

rule gzip_index_vcf:
    input:
        "{anything}.vcf"
    output:
        gz="{anything}.vcf.gz",
        tbi="{anything}.vcf.gz.tbi"
    params:
        error_out_file = "error_files/gzip.{anything}",
        run_time = "2:00:00",
        cores = "2"
    shell:
        """
        module load bcftools
        bgzip -c {input} > {output.gz}
        tabix {output.gz}
        """

rule make_index:
    input:
        circ_ref = "output/circ/{sample_label}_index.fa",
    output:
        fadict = "output/circ/{sample_label}_index.dict"
    params:
        error_out_file = "error_files/picardindex.{sample_label}",
        run_time = "1:00:00",
        cores = "1"
    shell:
        """
        module load picard
        picard CreateSequenceDictionary R={input.circ_ref} O={output.fadict}
        """
rule make_edit_table:
    input:
        vcf = "output/edits/{sample_label}.dp4.{strand}.vcf.gz",
        circ_ref = "output/circ/{sample_label}_index.fa",
        fadict = "output/circ/{sample_label}_index.dict"
    output:
        tsv = "output/edits/{sample_label}.dp4.{strand}.vcf.tsv"
    params:
        error_out_file = "error_files/variant_table.{sample_label}.{strand}",
        run_time = "6:00:00",
        cores = "2"
    conda:
        "envs/gatk4.yaml"
    shell:
        """
        gatk VariantsToTable \
            -R {input.circ_ref} \
            -V {input.vcf} \
            -F CHROM -F POS -F REF -F ALT -F FILTER -F TYPE -F DP -GF AD -GF ADF -GF ADR  \
            --show-filtered \
            -O {output.tsv}
        """


rule aggregate_pseudoreference:
    input:
        tsv="output/edits/{sample_label}.dp4.{strand}.vcf.tsv"
    output:
        all="output/edits/{sample_label}.dp4.{strand}.vcf.aggregated.tsv",
        nonzero="output/edits/{sample_label}.dp4.{strand}.vcf.aggregated.nonzero.tsv"
    params:
        error_out_file = "error_files/aggregate.{sample_label}.{strand}",
        run_time = "02:10:00",
        cores = "2",
        alt = lambda wildcards: ALT_fwd if wildcards.strand == 'pos' else ALT_rev
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/aggregate_pseudoreference_edit.py \
            {input.tsv} {params.alt} {output.all}
        """
rule filter_snp:
    input:
        tsv = "output/edits/{sample_label}.dp4.{strand}.vcf.aggregated.nonzero.tsv"
    output:
        tsv = "output/edits/{sample_label}.dp4.{strand}.vcf.aggregated.nonzero.snpfiltered.tsv"
    params:
        error_out_file = "error_files/filter_snp.{sample_label}.{strand}",
        run_time = "01:20:00",
        cores = "1"
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/filter_snp.py {input.tsv} {output}
        """



# rule split_by_contig:
#     input:
#         vcf="output/edits/{sample_label}.dp4.{strand}.vcf.snpfiltered.tsv"
#     output:
#         by_circ= "output/edits/{sample_label}.dp4.{strand}.circ_level.tsv"
#     params:
#         error_out_file = "error_files/variant_table_sum_by_circ.{sample_label}.{strand}",
#         run_time = "2:00:00",
#         cores = "2",
#     conda:
#         "envs/metadensity.yaml"
#     shell:
#         """
#         python {SCRIPT_PATH}/edit_level_by_circ.py {input.vcf} {params.alt} {output}
#         """


rule sailor_testing:
    input:
        tsv = "output/edits/{sample_label}.dp4.{strand}.vcf.aggregated.nonzero.snpfiltered.tsv"
    output:
        "output/edits/{sample_label}.dp4.{strand}.vcf.aggregated.nonzero.snpfiltered.edit_score.tsv"
    params:
        error_out_file = "error_files/sailor_test.{sample_label}.{strand}",
        run_time = "6:00:00",
        cores = "2",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/sailor_betainc_testing.py {input.tsv} {output}
        """

rule pick_from_the_right_strand:
    input:
        pos="output/edits/{sample_label}.dp4.pos.vcf.aggregated.nonzero.snpfiltered.edit_score.tsv",
        neg="output/edits/{sample_label}.dp4.pos.vcf.aggregated.nonzero.snpfiltered.edit_score.tsv",
        gtf="output/{sample_label}.gtf"
    output:
        "output/edits/{sample_label}.final_edit_score.tsv"
    params:
        error_out_file = "error_files/select_right_strand.{sample_label}",
        run_time = "1:00:00",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/pick_circ_from_strand.py \
            {input.neg} \
            {input.pos} \
            {input.gtf} \
            {output}
        """

### Analysis
rule make_genome_size:
    ''' make genome size table for bedtool expand'''
    input:
        fai="output/circ/{sample_label}_index.fa.fai"
    output:
        genome_sizes = "output/{sample_label}.genome.sizes",
    params:
        error_out_file = "error_files/genome_sizes.{sample_label}",
        run_time = "1:00:00",
        cores = "1",
    shell:
        """
        awk '{{print $1, $2/2}}' OFS='\t' {input.fai} > {output.genome_sizes}
        """

rule fetch_sequence_around_edit:
    input:
        edits="output/edits/{sample_label}.final_edit_score.tsv",
        genome_sizes = "output/{sample_label}.genome.sizes",
        fa="output/circ/{sample_label}_index.fa",
        fai="output/circ/{sample_label}_index.fa.fai"
    output:
        bed="output/edits/{sample_label}.final_edit_score_bed",
        bed_expanded="output/edits/{sample_label}.final_edit_score_bed_expand",
        fa_expanded="output/edits/{sample_label}.final_edit_score.fa"
    params:
        error_out_file = "error_files/fetch_sequence_around.{sample_label}",
        run_time = "1:00:00",
        cores = "1",
        size = 100 # based on Kris Branan's paper
    shell:
        """
        awk 'NR>1 {{print $14,$13,$13+1,$2,$12,$15}}' OFS='\t' {input.edits} > {output.bed}
        module load bedtools
            bedtools slop -i {output.bed} -g {input.genome_sizes} -b {params.size} \
                | sort -k1,1 -k2,2n -k3,3n \
                > {output.bed_expanded}

            bedtools getfasta -s -fi {input.fa} -bed {output.bed_expanded} -name \
                > {output.fa_expanded}
        """

# method 1: significant edits, merge, fetch sequence around
rule significant_edits:
    input:
        edits="output/edits/{sample_label}.final_edit_score.tsv",
        genome_sizes = "output/{sample_label}.genome.sizes",
        fa="output/circ/{sample_label}_index.fa",
        fai="output/circ/{sample_label}_index.fa.fai"
    output:
        sig_edit="output/edits/{sample_label}.sig_edit_bed",
        edit_expanded = "output/edits/{sample_label}.expand_edit_bed",
        fa_expanded = "output/edits/{sample_label}.expand_edit_bed.fa"
    params:
        error_out_file = "error_files/collapse_edit.{sample_label}",
        run_time = "1:00:00",
        cores = "1",
    shell:
        """
            awk '$12 < 0.2 {{print $14,$13,$13+1,$2,$12,$15}}' OFS='\t' {input.edits} > {output.sig_edit}

            module load bedtools
            bedtools slop -i {output.sig_edit} -g {input.genome_sizes} -b 50 \
                | sort -k1,1 -k2,2n -k3,3n \
                | bedtools merge -s -c 4,5,6 -o  collapse,collapse,distinct > {output.edit_expanded}

            bedtools getfasta -s -fi {input.fa} -bed {output.edit_expanded} \
                > {output.fa_expanded}
        """

rule homer_expanded:
    input:
        foreground = "output/edits/{stamp_sample_label}.expand_edit_bed.fa",
        background = "output/edits/{ctrl_sample_label}.expand_edit_bed.fa"
    output:
        "output/edits/homer/{stamp_sample_label}.{ctrl_sample_label}.homer"
    params:
        error_out_file = "error_files/edits_homer.{stamp_sample_label}.{ctrl_sample_label}",
        run_time = "8:00:00",
        cores = "1",
    shell:
        """
        module load homer
        homer2 denovo -i {input.foreground} -b {input.background} -strand + -o {output}
        """

# analysis on genome coordinates
rule collapse_edits:
    input:
        "output/edits/{sample_label}.final_edit_score.tsv"
    output:
        "output/edits/expanded_bed/{sample_label}.expandbed"
    params:
        error_out_file = "error_files/collapse_edit.{sample_label}",
        run_time = "1:00:00",
        cores = "1",
    conda:
        "envs/metadensity.yaml"
    shell:
        """
        python {SCRIPT_PATH}/collapse_edits_into_bed.py \
            {input} \
            {output}
        """

rule homer_with_expanded_edits:
    input:
        stamp="output/edits/expanded_bed/{stamp_sample_label}.expandbed",
        # use either external apobec library or internal ones
        ctrl=lambda wildcards:"output/edits/expanded_bed/{ctrl_sample_label}.expandbed" if wildcards.ctrl_sample_label in config['sample_labels'] 
        else config['external_stamp_control'][wildcards.ctrl_sample_label]
    output:
        homer="output/edits/homer/{stamp_sample_label}.{ctrl_sample_label}/homerResults.html"
    params:
        error_out_file = "error_files/homer.{stamp_sample_label}.{ctrl_sample_label}",
        run_time = "1:00:00",
        cores = "1",
        outdir = lambda wildcards, output: str(Path(output.homer).parent)
    shell:
        """
        module load homer
        findMotifsGenome.pl {input.stamp} \
            hg38 \
            {params.outdir} \
            -bg {input.ctrl} \
            -size given -rna -len 5,6,7,8,9 -noknown
        """


