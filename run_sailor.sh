conda activate snakemake
cd /home/hsher/scratch/circSTAMP_iter1

snakemake -s /projects/ps-yeolab4/software/stamp/0.2.0/downloads/STAMP/workflow_sailor/Snakefile  \
    --use-singularity --singularity-args '--bind /oasis --bind /projects --bind /home/hsher ' \
    --configfile /home/hsher/projects/circSTAMP_pipe/sailor.json \
    --verbose --profile /projects/ps-yeolab4/software/stamp/0.1.0/downloads/STAMP/profiles/tscc/ \
    --latency 50 -p -j8 -n


# sailor for cirular RNA

snakemake -s /projects/ps-yeolab4/software/stamp/0.2.0/downloads/STAMP/workflow_sailor/Snakefile  \
    --use-singularity --singularity-args '--bind /oasis --bind /projects --bind /home/hsher ' \
    --configfile /home/hsher/scratch/circSTAMP_pipe/output/sailor_circ/sailor_circ.json \
    --verbose --profile /projects/ps-yeolab4/software/stamp/0.1.0/downloads/STAMP/profiles/tscc/ \
    --latency 50 -p -j8 -n