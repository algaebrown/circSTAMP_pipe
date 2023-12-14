# to do analysis with sailor outputs
#
# snakemake -s SnakeSailorAnalysis.py -j 36 --configfile config/tao.yaml --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" --use-conda
sample_labels, = glob_wildcards("output/sailor/{sample_label}.Aligned.sortedByCoord.out.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed")

regions = ['all', 'cds', 'distintron500', 'five_prime_utrs', 'proxintron500', 'three_prime_utrs']
rbp_labels = [s for s in sample_labels if 'APO' not in s]
apo_labels = [s.replace(config['RBP'], 'APO') for s in rbp_labels]

rule all:
    input:
        expand("output/sailor/9_edit_fraction_bedgraphs/{sample_label}.Aligned.sortedByCoord.out.bam.edit_fraction.bw",sample_label = sample_labels),
        expand("output/sailor/annotate/{sample_label}.filtered.annotate.bed",sample_label = sample_labels),
        expand("output/sailor/motif/{sample_label}.filtered.expand.svg",sample_label = sample_labels),
        expand("output/sailor/motif_apo_as_bg/{comparison}.{region}.homer",region = regions, 
        comparison = [f'{rbp_label}.vs.{apo_label}' for rbp_label, apo_label in zip(rbp_labels, apo_labels)]
        )
    output:
        "sailor_analysis.txt"
    params:
        error_out_file = "error_files/all",
        run_time = "0:15:00",
        cores = "2",
        memory = "10000",
    shell:
        """
        echo DONE > {output}
        """


# generate bigwig
rule bw:
    input:
        sailor="output/sailor/9_edit_fraction_bedgraphs/{sample_label}.Aligned.sortedByCoord.out.bam.edit_fraction.bedgraph",
    output:
        sorted="output/sailor/9_edit_fraction_bedgraphs/{sample_label}.Aligned.sortedByCoord.out.bam.edit_fraction.so.bedgraph",
        bw="output/sailor/9_edit_fraction_bedgraphs/{sample_label}.Aligned.sortedByCoord.out.bam.edit_fraction.bw"
    params:
        error_out_file = "error_files/bw.{sample_label}",
        run_time = "4:45:00",
        cores = "2",
        memory = "10000",
        chr_size = config['CHROM_SIZE']
    shell:
        """
        module load ucsctools
        sort -k1,1 -k2,2n {input.sailor} > {output.sorted}
        bedGraphToBigWig {output.sorted} {params.chr_size} {output.bw}
        """

rule filter_edits_and_expand:
    input:
        "output/sailor/{sample_label}.Aligned.sortedByCoord.out.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed",
    output:
        filtered = "output/sailor/filtered_bed/{sample_label}.filtered.bed",
        expand = "output/sailor/filtered_bed/{sample_label}.filtered.expand.bed",
    params:
        error_out_file = "error_files/filter_edit.{sample_label}",
        run_time = "1:45:00",
        cores = "2",
        memory = "10000",
        chr_size = config['CHROM_SIZE'],
        thres = config['edit_thres'],
        window_len = config['EXPAND_LEN']
    shell:
        """
        awk '{{ if ($4 > {params.thres}) print }}' {input}  > {output.filtered}
        module load bedtools;
        bedtools slop -i {output.filtered} -g {params.chr_size} -b {params.window_len} > {output.expand}
        """

rule annotate:
    input:
        peak= "output/sailor/{sample_label}.Aligned.sortedByCoord.out.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed"
    output:
         "output/sailor/annotate/{sample_label}.filtered.annotate.bed"
    params:
        run_time="6:00:00",
        error_out_file = "error_files/annotate.{sample_label}",
        cores = "1",
        gtf_db = config['GTFDB'],
        sps = config['ANNOTATOR_SPECIES']
    shell:
        """
        module load annotator
        annotator --input {input.peak} --output {output} --gtfdb {params.gtf_db} --species {params.sps}
        """

rule motif_analysis:
    input:
        peak="output/sailor/filtered_bed/{sample_label}.filtered.expand.bed"
    output:
        pickle="output/sailor/motif/{sample_label}.filtered.expand.pickle",
        svg="output/sailor/motif/{sample_label}.filtered.expand.svg",
        rbp_fa=expand("output/sailor/motif/{sample_label}.filtered.expand.homer/fasta/{sample_label}.filtered.expand.bed.{region}.real.fa",
            region = regions, sample_label = ['{sample_label}'])
    params:
        run_time="3:00:00",
        error_out_file = "error_files/motif.{sample_label}",
        cores = "1",
        fa = config['GENOMEFA'],
        sps = config['ANNOTATOR_SPECIES'],
        homer="output/sailor/motif/{sample_label}.filtered.expand.homer",
    shell:
        """
        module load eclipanalysis
        analyze_motifs --peaks {input.peak} --k 6 \
            --out_pickle_file {output.pickle} \
            --out_homer_dir {params.homer} \
            --out_file {output.svg} \
            --species {params.sps} \
            --genome_fasta {params.fa}

        """

rule motif_analysis_APO_as_bg:
    input:
        pickle="output/sailor/motif/{rbp_label}.filtered.expand.pickle",
        rbp_fa="output/sailor/motif/{rbp_label}.filtered.expand.homer/fasta/{rbp_label}.filtered.expand.bed.{region}.real.fa",
        pickle_apo="output/sailor/motif/{apo_label}.filtered.expand.pickle",
        apo_fa="output/sailor/motif/{apo_label}.filtered.expand.homer/fasta/{apo_label}.filtered.expand.bed.{region}.real.fa",
    output:
        "output/sailor/motif_apo_as_bg/{rbp_label}.vs.{apo_label}.{region}.homer"
    params:
        run_time="00:30:00",
        error_out_file = "error_files/motif_apo_as_bg.{rbp_label}",
        cores = "1",
    shell:
        '''
        module load homer
        homer2 denovo -i {input.rbp_fa} -b {input.apo_fa} -o {output}
        '''