#!/bin/bash
#SBATCH --job-name=${USER}_alignment                                            # Name of the job
#SBATCH --partition=courses                                                  # the account used for computational work
#SBATCH -N 1                                                                # number of nodes
#SBATCH -c 8                                                                # number of cpus-per-task (threads)
#SBATCH --mem 32G                                                           # memory pool for all cores
#SBATCH -t 4:00:00                                                          # time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL                                                # Get an email when the program completes or fails
#SBATCH --mail-user=${USER}@northeastern.edu                             # where to send the email
#SBATCH --out=/courses/BINF6430.202510/students/${USER}/logs/%x_%j.log   # captured stdout
#SBATCH --error=/courses/BINF6430.202510/students/${USER}/logs/%x_%j.err # captured stdin

set -e
module load OpenJDK/19.0.1
module load fastqc
module load star/2.7.11a
module load bowtie/2.5.2
module load samtools/1.19.2
module load hisat2/2.2.0
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8


################################################################################################
### HOUSEKEEPING ###
################################################################################################
BASE_DIR=/courses/BINF6430.202510
DATA_DIR=${BASE_DIR}/shared/gencode
REFERENCE_FILE=${DATA_DIR}/GRCh38.primary_assembly.genome.fa
INDEX=$(basename ${REFERENCE_FILE} .fa)
RESULTS=${BASE_DIR}/students/${USER}/${SLURM_JOB_ID}_results
mkdir -p ${RESULTS}
GENOME_DIR=${BASE_DIR}/students/${USER}/STAR
SAMPLE_DIR=${BASE_DIR}/data/ReaganData/Reagan_PE85_TakaraPicoV2_HC_CM_10042022
find "$SAMPLE_DIR" -type f -name "*.gz" > ${RESULTS}/samplemanifest.txt
SAMPLE_MANIFEST=${RESULTS}/samplemanifest.txt


################################################################################################
### GENOME BUILD ###
################################################################################################
mkdir -p ${GENOME_DIR}
bowtie2-build --threads ${SLURM_CPUS_PER_TASK} ${REFERENCE_FILE} ${GENOME_DIR}/${INDEX}


################################################################################################
### ALIGNMENT PYTHON MODULE ###
################################################################################################
python STAR_alignment.py -c ${SAMPLE_MANIFEST} -g ${GENOME_DIR} -o ${RESULTS}