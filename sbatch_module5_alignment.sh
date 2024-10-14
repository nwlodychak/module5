#!/bin/bash
#SBATCH --job-name=wlodychak.s_alignment                                            # Name of the job
#SBATCH --partition=courses                                                  # the account used for computational work
#SBATCH -N 1                                                                # number of nodes
#SBATCH -c 8                                                                # number of cpus-per-task (threads)
#SBATCH --mem 32G                                                           # memory pool for all cores
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
    "bowtie/2.5.2"
    "samtools/1.19.2"
    "anaconda3/2021.11"
)

# load modules
for module in "${modules[@]}"; do
    module load "$module"
done

source activate BINF-12-2021
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
MULTIQC_PATH=/courses/BINF6430.202510/shared/multiqc-latest.sif
alias MULTIQC="singularity run ${MULTIQC_PATH} multiqc"

echo "Modules loaded and environment set up successfully."


################################################################################################
### HOUSEKEEPING ###
################################################################################################
USER=${1}
BASE_DIR=/courses/BINF6430.202510
RESULTS=${BASE_DIR}/students/${USER}/${SLURM_JOB_ID}_results
mkdir -p ${RESULTS}
GENOME_DIR=${BASE_DIR}/students/${USER}/STAR/
SAMPLE_DIR=${BASE_DIR}/data/ReaganData/Reagan_PE85_TakaraPicoV2_HC_CM_10042022
find "$SAMPLE_DIR" -type f -name "*.gz" > ${RESULTS}/samplemanifest.txt
SAMPLE_MANIFEST=${RESULTS}/samplemanifest.txt

################################################################################################
### ALIGNMENT PYTHON MODULE ###
################################################################################################
python STAR_alignment.py -c ${SAMPLE_MANIFEST} -g ${GENOME_DIR} -o ${RESULTS}