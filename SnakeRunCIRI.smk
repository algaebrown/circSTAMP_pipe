import pandas as pd
import warnings
try:
    manifest = pd.read_csv(config['menifest'])
except Exception as e:
    print(e)
locals().update(config)

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

################## Differential Expression ##################
# This only works for a single replicate shit
rule run_ciri_de_compare_unadjusted:
    input:
        gtf_control = "output/{sample_label_control}.gtf",
        gtf_case = "output/{sample_label_case}.gtf"
    output:
        "output/ciri_de/{sample_label_control}.{sample_label_case}.gtf.tsv"
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "8:00:00",
        cores = "1",
        memory = 320000,
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        CIRI_DE -n {input.gtf_control} \
        -c {input.gtf_case} \
        -o {output} 
        """

def find_file(paths, sample_label):
    for path in paths:
        fname = Path(path)/f"output/{sample_label}.gtf"
        fname2 = Path(path)/f"output/gene/{sample_label}_out.gtf"
        print(fname, fname2)
        if fname.exists() and fname2.exists():
            return fname, fname2
    raise Exception(f'Unable to find GTF and _out.gtf for {sample_label}')
    



# All the rest makes replicate differential expressions
rule make_list:
    input:
        lambda wildcards: expand("output/{sample_label}.gtf", 
        sample_label =  [s for s in config['circ_de'][wildcards.comparison]['case']+config['circ_de'][wildcards.comparison]['ctrl']
        if s in sample_labels
        ]),
        lambda wildcards: expand("output/gene/{sample_label}_genes.list", 
        sample_label = [ s for s in config['circ_de'][wildcards.comparison]['case']+config['circ_de'][wildcards.comparison]['ctrl']
        if s in sample_labels
        ])
    output:
        lst="output/circ_de_reps/{comparison}.lst",
        sample_lst="output/circ_de_reps/{comparison}.samplegene.lst"
    params:
        error_out_file = "error_files/ciri.de.{comparison}",
        run_time = "10:00",
        cores = "1",
        memory = 1000
    run:
        import pandas as pd
        lst_df = []
        lst_sample_gene_df  = []
        for i,sample_label in enumerate(config['circ_de'][wildcards.comparison]['case']):
            fname, fname2 = find_file(config['external_de_paths']+[config['workdir']], sample_label)
            lst_df.append([f'case_rep{i+1}', fname, 'T'])
            lst_sample_gene_df.append([f'case_rep{i+1}', fname2])
                        
        for i,sample_label in enumerate(config['circ_de'][wildcards.comparison]['ctrl']):
            fname, fname2 = find_file(config['external_de_paths']+[config['workdir']], sample_label)
            lst_df.append([f'ctrl_rep{i+1}', fname, 'C'])
            lst_sample_gene_df.append([f'ctrl_rep{i+1}', fname2])
        
        df = pd.DataFrame(lst_df, columns = ['sample', 'path', 'group'])
        df_stringtie = pd.DataFrame(lst_sample_gene_df, columns = ['sample', 'path'])
        df.to_csv(output.lst, sep = '\t', header = False, index = False)
        df_stringtie.to_csv(output.sample_lst, sep = '\t', header = False, index = False)


rule prep_ciriquant_de:
    input:
        lst = 'output/circ_de_reps/{comparison}.lst'
    output:
        lib_info = 'output/circ_de_reps/{comparison}.csv',
        circ = 'output/circ_de_reps/{comparison}.circ.csv',
        bsj = 'output/circ_de_reps/{comparison}.bsj.csv',
        ratio = 'output/circ_de_reps/{comparison}.ratio.csv',
    params:
        error_out_file = "error_files/ciri_de.{comparison}",
        run_time = "2:00:00",
        cores = "1",
        memory = "160000"
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        prep_CIRIquant -i {input} \
            --lib {output.lib_info}\
            --circ {output.circ} \
            --bsj {output.bsj} \
            --ratio {output.ratio}
        """

rule prep_stringtie:
    input:
        lst = 'output/circ_de_reps/{comparison}.samplegene.lst'
    output:
        genecnt = 'output/circ_de_reps/{comparison}.genecount_matrix.csv',
    params:
        error_out_file = "error_files/ciri_de_rep_stringtie.{comparison}",
        run_time = "2:00:00",
        cores = "1",
        memory = "160000"
    container:
        "docker://mortreux/ciriquant:v1.1.2"
    shell:
        """
        python {SCRIPT_PATH}/prepDE.py \
            -i {input} \
            -g {output.genecnt}
        """

# this env is hard to install as hell
rule run_ciri_de_compare_unadjusted:
    input:
        lib_info = rules.prep_ciriquant_de.output.lib_info,
        bsj = rules.prep_ciriquant_de.output.bsj,
        gene = rules.prep_stringtie.output.genecnt
    output:
        circ = "output/circ_de_reps/{comparison}.gtf.tsv",
    params:
        error_out_file = "error_files/ciri_de_rep.{comparison}",
        run_time = "6:00:00",
        cores = "1",
        memory = "160000"
    container:
        "docker://algaebrown/ciriquant_edger:1.1.2"
    shell:
        """
        /opt/conda/lib/python2.7/site-packages/libs/CIRI_DE.R \
            --lib {input.lib_info} \
            --bsj {input.bsj} \
            --gene {input.gene} \
            --out {output.circ} 
        """

