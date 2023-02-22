# to run differential expression analysis
# from gtf files
#snakemake -s de.smk --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" -j 15 --use-conda --conda-prefix /home/hsher/snakeconda -n
import pandas as pd
from pathlib import Path
import os

workdir: '/home/hsher/scratch/circ_de'

################# COMPARISON #############
# control, case
comparison = [
    ['circseq-bm-arr', 'circseq-bm-rar'],
    ['circseq-bm-arr', 'circseq-bm-rar11-nxt'],
    ['circseq-bm-rar', 'circseq-bm-rar11-nxt'],
    ['HEK_rar11', 'EV_rar11'],
    ['HEK_rar', 'HEK_rar11'],
    ['EV_rar11', 'TOB2_rar11'],
    ['EV_rar11', 'RBM15_rar11_down'],
]

try:
    os.mkdir('error_out_file')
except:
    pass
################ MANIFEST ###############
# RNase+, RNase-
manifest = [
['/home/hsher/scratch/circ_truseq/output/circseq-bm-rar.gtf', '/home/hsher/scratch/circ_truseq/output/circseq-bm-rz.gtf'],
['/home/hsher/scratch/circ_truseq/output/circseq-bm-arr.gtf', '/home/hsher/scratch/circ_truseq/output/circseq-bm-rz.gtf'],
['/home/hsher/scratch/circ_nextera/output/circseq-bm-rar11-nxt.gtf', '/home/hsher/scratch/circ_truseq/output/circseq-bm-rz.gtf'],
['/home/hsher/scratch/circ_nextera_iter2/output/HEK_rar11.gtf', '/home/hsher/scratch/circ_truseq_iter2/output/EV_rz.gtf'],
['/home/hsher/scratch/circ_truseq_iter2/output/HEK_rar.gtf', '/home/hsher/scratch/circ_truseq_iter2/output/EV_rz.gtf'],
['/home/hsher/scratch/circ_nextera_iter2/output/EV_rar11.gtf','/home/hsher/scratch/circ_truseq_iter2/output/EV_rz.gtf'],
['/home/hsher/scratch/circ_nextera_iter2/output/TOB2_rar11.gtf','/home/hsher/scratch/circ_truseq_iter2/output/TOB2_rz.gtf'],
['/home/hsher/scratch/circ_nextera_iter2_rbm15_down/output/RBM15_rar11_down.gtf', '/home/hsher/scratch/circ_truseq_iter2/output/RBM15_rz.gtf']
]

manifest = pd.DataFrame(manifest, columns = ['RNase+', 'RNase-'])
manifest['sample_label'] = manifest['RNase+'].apply(lambda path: Path(path).name.replace('.gtf',''))

# find fq_r1
def find_fq(gtf_path, read = 1):
    uid = Path(gtf_path).name.replace('.gtf', '')
    output_dir = Path(gtf_path).parent

    if read == 1:
        fq_path = output_dir / 'fastqs' / f'{uid}_1.Tr.fq.gz'
    else:
        fq_path = output_dir / 'fastqs' / f'{uid}_2.Tr.fq.gz'

    return fq_path


manifest['fq1'] = manifest['RNase-'].apply(lambda path: find_fq(path, read = 1))
manifest['fq2'] = manifest['RNase-'].apply(lambda path: find_fq(path, read = 2))

print(manifest.applymap(os.path.isfile).sum(axis = 1))
print(manifest.applymap(os.path.isfile).sum(axis = 0))

################ CIRI CONFIG #########################
config = {'CIRICONFIG': '/home/hsher/scratch/circSTAMP_pipe/ciriconfig_full.yaml',
'LIBRARY_TYPE': 2,
'BWA_INDEX': '/projects/ps-yeolab5/hsher/Tao_circSTAMP/bwa_gencode35',
'HISAT_INDEX': '/projects/ps-yeolab5/hsher/Tao_circSTAMP/hisat_gencode35',
'READ1_ADAPTOR':None,
'READ2_ADAPTOR':None,
'STAR_INDEX':'/projects/ps-yeolab5/hsher/Tao_circSTAMP/star_2_7_gencode35'
}

############## RULES #################
rule all:
    input:
        expand("output/compare_to_rz/{sample_label}.gtf.tsv", sample_label = manifest['sample_label'].tolist()),
        expand("output/unadjusted_comparison/{sample_label_combination}.gtf.tsv", 
            sample_label_combination = [c[0]+'.'+c[1] for c in comparison],
            ),
        expand("output/adjusted_comparison/{sample_label_combination}.gtf.tsv", 
            sample_label_combination = [c[0]+'.'+c[1] for c in comparison],
            ),


rule run_ciri_polyA_rnase_coorection:
    input:
        bwa_index=config['BWA_INDEX']+'.amb',
        hisat_index=config['HISAT_INDEX']+'.1.ht2',
        total_rna_read1 = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label, 'fq1'],
        total_rna_read2 = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label, 'fq2'],
        yaml=config['CIRICONFIG'],
        rnase_treated_gtf = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label, 'RNase+'],
    output:
        "output/adjusted_gtf/{sample_label}.gtf",
    params:
        name="{sample_label}",
        error_out_file = "error_files/ciri.rnasecorr.{sample_label}",
        outdir='output/adjusted_gtf',
        run_time = "48:00:00",
        cores = "3",
        library_type=config['LIBRARY_TYPE'],
        #rnase=lambda wildcards: '--RNaseR output/Rnase' if manifest.loc[manifest['Sample']==wildcards.sample_label, 'Rnase'].iloc[0] else ''
        # https://github.com/bioinfo-biols/CIRIquant/issues/4
    conda:
        "/home/hsher/projects/circSTAMP_pipe/envs/ciriquant.yaml"
    shell:
        """
        CIRIquant -t {params.cores} \
            -1 {input.total_rna_read1} \
            -2 {input.total_rna_read2} \
            --config {input.yaml} \
            --library-type {params.library_type} \
            -o {params.outdir} \
            -p {params.name} \
            --RNaseR {input.rnase_treated_gtf}
        """

############## RULES #################
rule run_ciri_de_compare_to_rz:
    input:
        gtf_control = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label, 'RNase-'],
        gtf_case = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label, 'RNase+'],
    output:
        "output/compare_to_rz/{sample_label}.gtf.tsv"
    params:
        error_out_file = "error_files/ciri.de.{sample_label}",
        run_time = "2:00:00",
        cores = "1",
    conda:
        "/home/hsher/projects/circSTAMP_pipe/envs/ciriquant.yaml"
    shell:
        """
        CIRI_DE -n {input.gtf_control} \
        -c {input.gtf_case} \
        -o {output}
        """

rule run_ciri_de_compare_unadjusted:
    input:
        gtf_control = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label_control, 'RNase+'],
        gtf_case = lambda wildcards: manifest.loc[manifest['sample_label']==wildcards.sample_label_case, 'RNase+'],
    output:
        "output/unadjusted_comparison/{sample_label_control}.{sample_label_case}.gtf.tsv"
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "2:00:00",
        cores = "1",
    conda:
        "/home/hsher/projects/circSTAMP_pipe/envs/ciriquant.yaml"
    shell:
        """
        CIRI_DE -n {input.gtf_control} \
        -c {input.gtf_case} \
        -o {output}
        """

rule run_ciri_de_compare_adjusted:
    input:
        gtf_control = "output/adjusted_gtf/{sample_label_control}.gtf",
        gtf_case = "output/adjusted_gtf/{sample_label_case}.gtf",
    output:
        "output/adjusted_comparison/{sample_label_control}.{sample_label_case}.gtf.tsv"
    params:
        error_out_file = "error_files/ciri.de.{sample_label_control}.{sample_label_case}",
        run_time = "8:00:00",
        cores = "2",
    conda:
        "/home/hsher/projects/circSTAMP_pipe/envs/ciriquant.yaml"
    shell:
        """
        CIRI_DE -n {input.gtf_control} \
        -c {input.gtf_case} \
        -o {output} \
        -t 4
        """
