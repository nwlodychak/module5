#!/bin/bash
#SBATCH --job-name=wlodychak.s_genomeBuild                                  # name of job
#SBATCH --partition=courses                                                 # the account used for computational work
#SBATCH -N 1                                                                # number of nodes
#SBATCH -c 16                                                                # number of cpus-per-task (threads)
#SBATCH --mem 50G                                                           # memory pool for all cores
#SBATCH -t 4:00:00                                                          # time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL                                                # Get an email when the program completes or fails
#SBATCH --mail-user=wlodychak.s@northeastern.edu                             # where to send the email
#SBATCH --out=/courses/BINF6430.202510/students/wlodychak.s/logs/%x_%j.log   # captured stdout
#SBATCH --error=/courses/BINF6430.202510/students/wlodychak.s/logs/%x_%j.err # captured stdin

echo "Loading Modules"
set -e

modules=(
    "OpenJDK/19.0.1"
    "fastqc"
    "star/2.7.11a"
    "anaconda3/2021.11"
)

# load modules
for module in "${modules[@]}"; do
    module load "$module"
done

source activate BINF-12-2021
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8

echo "Modules loaded and environment set up successfully."


################################################################################################
### HOUSEKEEPING ###
################################################################################################
USER=${1}
BASE_DIR=/courses/BINF6430.202510
DATA_DIR=${BASE_DIR}/shared/gencode
REFERENCE_FILE=${DATA_DIR}/GRCh38.primary_assembly.genome.fa
GTF_FILE=${DATA_DIR}/gencode.v46.annotation.gtf
INDEX=$(basename ${REFERENCE_FILE} .fa)
GENOME_DIR=${BASE_DIR}/students/${USER}/STAR


################################################################################################
### GENOME BUILD ###
################################################################################################
mkdir -p ${GENOME_DIR}
STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
     --runMode genomeGenerate \
     --genomeDir ${GENOME_DIR} \
     --genomeFastaFiles ${REFERENCE_FILE} \
     --sjdbGTFfile ${GTF_FILE} \
     --sjdbOverhang 84