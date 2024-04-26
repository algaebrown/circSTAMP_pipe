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

# cleanup all temp directory
cd ~/projects/circSTAMP_pipe
for f in config/*.yaml
do
snakemake -s SnakeMain.smk \
    -j 12 \
    --configfile config/tao_trueseq2.yaml \
    --cluster "qsub -l walltime={params.run_time} -l nodes=1:ppn={params.cores} -e {params.error_out_file} -q home-yeo" \
    --use-conda \
    --conda-prefix /home/hsher/snakeconda \
    --delete-temp-output
done