#!/bin/bash
#SBATCH --job-name=wlodychaks_alignment                                            # Name of the job
#SBATCH --partition=courses                                                  # the account used for computational work
#SBATCH -N 1                                                                # number of nodes
#SBATCH -c 4                                                                # number of cpus-per-task (threads)
#SBATCH --mem 32G                                                           # memory pool for all cores
#SBATCH -t 4:00:00                                                          # time (HH:MM:SS)
#SBATCH --mail-type=END,FAIL                                                # Get an email when the program completes or fails
#SBATCH --mail-user=wlodychak.s@northeastern.edu                             # where to send the email
#SBATCH --out=/courses/BINF6430.202510/students/${USER}/logs/%x_%j.log   # captured stdout
#SBATCH --error=/courses/BINF6430.202510/students/${USER}/logs/%x_%j.err # captured stdin

set -e
module load OpenJDK/19.0.1
module load fastqc
module load star/2.7.11a
module load bowtie/2.5.2
module load singularity
module load hisat2/2.2.0
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8


################################################################################################
### HOUSEKEEPING ###
################################################################################################
BASE_DIR=/courses/BINF6430.202510
DATA_DIR=${BASE_DIR}/shared/gencode
RESULTS=${BASE_DIR}/students/${USER}/results
GENOME_DIR=${BASE_DIR}/students/${USER}/STAR
SAMPLE_DIR=${BASE_DIR}/data/ReaganData/Reagan_PE85_TakaraPicoV2_HC_CM_10042022
find "$SAMPLE_DIR" -type f -name "*.gz" > ${RESULTS}/samplemanifest.txt


################################################################################################
### GENOME BUILD ###
################################################################################################
mkdir -p ${GENOME_DIR}
bowtie2-build --threads ${SLURM_CPUS_PER_TASK} {genome} {index}


################################################################################################
### STAR MODULE ###
################################################################################################
shopt -s expand_aliases
alias STAR='singularity run -B "/courses:/courses" \
			/courses/BINF6430.202510/shared/singularity_containers/star-latest.sif STAR'


################################################################################################
### ALIGNMENT PYTHON MODULE ###
################################################################################################
## TODO parsing the sample set and calling index files for genome
python STAR_alignment.py -c ${SAMPLEMANIFEST} -g ${GENOME_DIR} -t ${SLURM_CPUS_PER_TASK}

################################################################################################
### SAMTOOLS ###
################################################################################################
SAMTOOLS_PATH=/courses/BINF6430.202510/shared/samtools-latest.sif
alias SAMTOOLS="singularity run ${SAMTOOLS_PATH} samtools"
## TODO call SAMTOOLS output