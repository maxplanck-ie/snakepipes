main_dir_path="/home/kilpert/git/snakemake_workflows/"
indir="/data/manke/kilpert/Nicola/A277_TRR_KD/01_rna-seq-qc/"
genome="dm6"
outdir="/data/manke/kilpert/Nicola/A277_TRR_KD/12.1_rMATS"
snakefile="/home/kilpert/git/snakemake_workflows/workflows/rMATS/Snakefile"

user_configs="--config maindir=$main_dir_path indir=$indir genome=$genome subproject_names=$subproject_names"


module load snakemake

snakemake $user_configs -p --directory $outdir --snakefile $snakefile --cluster 'SlurmEasy --threads {threads} --log {cluster_logs_dir}' --jobs 999 --jobname '{rulename}.{jobid}.snakejob'
