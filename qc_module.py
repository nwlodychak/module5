import subprocess
import argparse
import os
import logging
import time


def fastqc(fq1, outdir):
    """
    Generate fastqc files
    :param fq1: full path for fastq file
    :param outdir: where you want the fastqc files to go
    :return: .html and .zip of fastq
    """
    outdir = f'{outdir}/fastqc'
    os.makedirs(outdir, exist_ok = True)
    try:
        run_command(f'fastqc {fq1} --outdir {outdir}')
    except FileNotFoundError as e:
        logging.error(f"An error occurred: {e}")
        raise


 def multiqc(outdir):
    """
    Generate multiqc files
    :param outdir: results location
    :return: .html and .zip of fastq
    """
    outdir = f'{outdir}/multiqc'
    fastqc_indir = f'{outdir}/fastqc'
    os.makedirs(outdir, exist_ok = True)
    try:
        run_command(f'multiqc -O {outdir} {fastqc_indir}')
    except FileNotFoundError as e:
        logging.error(f"An error occurred: {e}")
        raise

def worker(sample_id, read1, read2, outdir):
# fastqc block
     try:
         logging.info(f'FASTQC {sample_id}')
         fastqc(read1, outdir)
         fastqc(read2, outdir)
         time.sleep(5)
     except Exception as e:
        logging.error(f"Exception occurred during FASTQC: {e}")
