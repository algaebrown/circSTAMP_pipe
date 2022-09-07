import json
import sys
import os
# python ../../output_to_json.py /home/hsher/scratch/circSTAMP_pipe/output/bams/ /home/hsher/gencode_coords/GRCh38.primary_assembly.genome.fa /home/hsher/scratch/circSTAMP_pipe/output/sailor/ /home/hsher/scratch/circSTAMP_pipe/output/sailor/sailor.json
if __name__=='__main__':
    bam_dir = sys.argv[1]
    genome_fa = sys.argv[2]
    outdir = sys.argv[3]
    json_out = sys.argv[4]

    print(json_out)

    all_bams = [f for f in os.listdir(bam_dir) if f.endswith('Aligned.sortedByCoord.out.bam')]

    json_dict = {"samples_path": bam_dir, 
    "samples": all_bams,
    "remove_duplicates": False, 
    "junction_overhang": 10, 
    "edge_mutation": 5, 
    "non_ag": 1, 
    "reverse_stranded": True, 
    "reference_fasta": genome_fa, 
    "dp": "DP4", 
    "min_variant_coverage": 5, 
    "known_snps": "/projects/ps-yeolab3/ekofman/ReferenceData/hg38/hg38_dbSNP_knownSNPs.bed3", 
    "edit_type": "ct", 
    "edit_fraction": 0.01, 
    "alpha": 0, 
    "beta": 0, 
    "keep_all_edited": False, 
    "output_dir": outdir, 
    "snakemake_dir_path": "/projects/ps-yeolab3/ekofman/sc_STAMP_pipeline/STAMP/", 
    "config_path": json_out}

    with open(json_out, 'w') as f:
        json.dump(json_dict, f)