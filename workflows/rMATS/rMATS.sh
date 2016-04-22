main_dir_path="/home/kilpert/git/snakemake_workflows/"
indir="/data/manke/kilpert/Nicola/A277_TRR_KD/01_rna-seq-qc/"
genome="dm6"
outdir="/data/manke/kilpert/Nicola/A277_TRR_KD/12.1_rMATS"
snakefile="/home/kilpert/git/snakemake_workflows/workflows/rMATS/Snakefile"

user_configs="--config maindir=$main_dir_path indir=$indir genome=$genome subproject_names=$subproject_names"

local_cores="4"
max_jobs="10"

module load snakemake

snakemake $user_configs -p --local-cores $local_cores --jobs $max_jobs --directory $outdir --snakefile $snakefile --cluster-config $main_dir_path/shared/cluster.yaml --cluster 'qsub -pe smp cluster.n -o {cluster_logs_dir} -e {cluster_logs_dir}' --jobname '{rulename}.{jobid}.snakejob'
