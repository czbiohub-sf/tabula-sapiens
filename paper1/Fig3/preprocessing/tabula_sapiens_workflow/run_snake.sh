#!/bin/bash

#name shell variables for calling in snakemake and sbatch
DATETIME=$(date "+%Y_%m_%d_%H_%M_%S")
RESTART=0

SNAKEFILE=Snakefile.smk
EMAIL=mswift2@stanford.edu

CLUSTER_CONFIG=config/slurm_config.yaml

#Snakemake config
NJOBS=200
WAIT=120

source $HOME/miniconda3/etc/profile.d/conda.sh

conda activate snakemake

mkdir -p snakemake_logs/slurm_logs/

#log file for process that calls snakemake
SBATCH_LOGFILE=snakemake_logs/cluster.$DATETIME.log
SBATCH_LOGFILE_ERR=snakemake_logs/cluster.$DATETIME.log.err

#Snakemake log file
LOGFILE=snakemake_logs/snakamaka.$DATETIME.log
LOGFILE_ERR=snakemake_logs/snakamaka.$DATETIME.log.err

# for a specific cell / output
TARGET=''


if [ $# -eq 0 ]
  then
    # Dry run snakemake
    snakemake -s $SNAKEFILE $TARGET --use-conda --keep-target-files --rerun-incomplete -n -r --quiet --keep-going
elif [ $1 = "unlock" ]
    then
        snakemake -s $SNAKEFILE $TARGET -F --rerun-incomplete --unlock --cores 1
    
elif [ $1 = "dryrun" ]
    # Dry run snakemake and print shell cmds 
    then
        snakemake \
            -s $SNAKEFILE $TARGET \
            --use-conda \
            --keep-target-files \
            --rerun-incomplete \
            -n \
            -r \
            -k \
            --printshellcmds

elif [ $1 = "dryforce" ]
    # Dry run snakemake and print shell cmds 
    then
        snakemake \
            -s $SNAKEFILE $TARGET \
            --use-conda \
            --keep-target-files \
            --rerun-incomplete \
            -n \
            -r \
            -k \
	    -F \
            --printshellcmds

elif [ $1 = "force" ]
    # Dry run snakemake and print shell cmds 
    then
        snakemake \
            -s $SNAKEFILE $TARGET \
            --use-conda \
            --keep-target-files \
            --rerun-incomplete \
            -n \
            -r \
            -k \
            --printshellcmds \
	    -F

elif [ $1 = "dag" ]
    # Dry run snakemake and print shell cmds 
    then
        snakemake \
            -s $SNAKEFILE $TARGET \
            --use-conda \
            --keep-target-files \
            --rerun-incomplete \
            -n \
            -r \
            -k \
            --printshellcmds \
	    --debug-dag

elif [ $1 = "snakemake" ]
    then
  # Run snakemake
    echo 'running snakemake'
    snakemake \
        -s $SNAKEFILE $TARGET \
        --use-conda \
        --cluster-config ${CLUSTER_CONFIG} \
        --cluster "sbatch \
                  --job-name={cluster.name} \
                  --time {cluster.time} \
                  --mem={cluster.mem} \
                  --ntasks={cluster.ntasks} \
                  --cpus-per-task={cluster.cpus-per-task} \
                  --partition={cluster.partition} \
                  --output={cluster.output} \
                  --error={cluster.error}" \
        --keep-target-files \
        --rerun-incomplete \
        -j $NJOBS \
        -w $WAIT \
        -k \
        --restart-times $RESTART \
        --keep-going

elif [ $1 = "sbatch" ]
    # Run snakemake as a job, good for long workflows
    then
    sbatch \
        --ntasks=1 \
        --cpus-per-task=1 \
        --mem=8000 \
        --mail-user=$EMAIL \
        --time 0-12 \
        -p normal,owners,quake \
        -o $SBATCH_LOGFILE \
        -e $SBATCH_LOGFILE_ERR \
        run_snake.sh snakemake

else
    echo "wrong option"
fi
