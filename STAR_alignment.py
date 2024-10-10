import argparse
import subprocess
import multiprocessing
import os
import glob
import logging
import time
import json
from os.path import isfile

import pandas as pd
from pysam import fastq


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',
                        '--input',
                        required = True)
    parser.add_argument('-c',
                        '--config',
                        default = '.',
                        required = True,
                        help = 'Config file containing samples to analyze')
    parser.add_argument('-g',
                        '--genome',
                        default = '.',
                        required = True,
                        help = 'Raw Genome file location')
    parser.add_argument('-o',
                        '--outdir',
                        default = '.',
                        required = True,
                        help = 'Where to send all files.')
    return parser.parse_args()


def run_command(command):
    """
    Run a bash command
    :param command:
    """
    subprocess.run(command, shell = True, check = True, stdout = subprocess.PIPE)


def fastqc(fq1, indir):
    """
    Generate fastqc files
    :param fq1:
    :return: .html and .zip of fastq
    """
    outdir = f'{indir}/results/fastqc'
    infile = f'{indir}/fastq/{fq1}'
    os.makedirs(outdir, exist_ok = True)
    try:
        run_command(f'fastqc {infile} --outdir {outdir}')
    except FileNotFoundError as e:
        logging.error(f"An error occurred: {e}")
        raise


def multiqc(indir, outdir):
    """
    Generate multiqc files
    :param indir: fastqc indir location
    :return: .html and .zip of fastq
    """
    outdir = f'{indir}/results/multiqc'
    fastqc_indir = f'{indir}/results/fasqc'
    os.makedirs(outdir, exist_ok = True)
    try:
        run_command(f'multiqc --outdir {outdir} {fastqc_indir}')
    except FileNotFoundError as e:
        logging.error(f"An error occurred: {e}")
        raise


def trim(sample, indir):
    """
    Standard adapter trimming function
    :param sample: sample basename
    :return: trimmed samples for R1 and R2
    """
    logging.basicConfig(level = logging.INFO, format = '%(message)s', filename = 'logs/trimming.log')
    inraw = f'{indir}/fastq'
    outtrim = f'{indir}/trimmed'
    trim_ext1 = "_R1.trimmed.fastq"
    trim_ext2 = "_R2.trimmed.fastq"

    raw1 = glob.glob(os.path.join(inraw, sample + "*_R1*"))[0]
    raw2 = glob.glob(os.path.join(inraw, sample + "*_R2*"))[0]
    trim1 = os.path.join(outtrim, sample + trim_ext1)
    trim2 = os.path.join(outtrim, sample + trim_ext2)

    os.makedirs(outtrim, exist_ok = True)

    print(raw1)
    print(raw2)

    assert raw1 is not None and raw2 is not None, \
        (logging.error(f"Read 1 is {raw1}\n")) and logging.error(f"Read 2 is {raw2}\n")

    try:
        logging.info(f"Trimming {sample}")
        command = (f"trimmomatic PE \
                    {raw1} {raw2} \
                    {trim1} {trim2} \
                    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15")
        run_command(command)
    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise

    return trim1, trim2


def align_STAR(genome, fq1, fq2, outdir, aligner):
    """
    build a genomic index for any genome
    :param genome: indexed genome of interest
    :param fq1: trimmed RNASeq read1
    :param fq2: trimmed RNASeq read2
    :param outdir: where do you want the alignment files to go?
    :param aligner: type of aligner to use
    :return: index build assignment location.
    """
    assert isfile(genome) and isfile(fq1) and isfile(fq2), f"Genome, fq1 or fq2 does not exist"
    index = os.path.basename(genome)
    sample_id = os.path.basename(fq1)
    if aligner == "STAR":
        try:
            logging.info(f"STAR alignment for {fq1} {fq2}")
            command = (f"STAR --genomeDir {genome} \
                            --runThreadN {SBATCH_cores} \
                            --readFilesIn {fq1} {fq2} \
                            --outFileNamePrefix {sample_id} \
                            --outSAMtype BAM SortedByCoordinate \
                            --outSAMunmapped Within \
                            --outSAMattributes Standard ")

        except Exception as e:
            logging.error(f"An error occurred: {e}")
            raise

    elif aligner == "HiSAT":
        try:
            logging.info(f"HiSAT alignment for {fq1} {fq2}")
            command = (f"hisat2 -f -x {genome} \
                        -1 {fq1} \
                        -2 {fq2} \
                        -S {sample_id}.sam ")

        except Exception as e:
            logging.error(f"An error occurred: {e}")
            raise
    else:
        logging.error(f"Unknown aligner {aligner}")


def worker(sample, indir):
    """
    worker function for multitprocessing algner
    :param sample:
    :param indir:
    :return:
    """

    logging.basicConfig(level = logging.INFO, format = '%(message)s', filename = f'{logs}/run.log')
    try:
        logging.info(f'Trimming {sample}')
        trim1, trim2 = trim(sample, indir)
        time.sleep(5)
        logging.info(f'Removing preprocessing files')
        remove_files = (trim1, trim2)
        for file in remove_files:
            logging.info(f'Removing {file}')
            os.remove(file)

    except Exception as e:
        logging.error(f"Exception occurred: {e}")


def get_samples(files):
    sample_dict = {}
    for file in files:
        file_name = os.path.basename(file)
        sample_name = file_name.split('_R')[0]
        read_indicator = file_name.split('_R')[-1]
        if sample_name not in sample_dict:
            sample_dict[sample_name] = {'R1': file, 'R2': None}
        if read_indicator.startswith('1'):
            sample_dict[sample_name]['R1'] = file
        elif read_indicator.startswith('2'):
            sample_dict[sample_name]['R2'] = file

    return sample_dict


def main():
    args = get_args()
    logs = f'{args.outdir}/logs'
    os.makedirs(logs, exist_ok = True)
    logging.basicConfig(level = logging.INFO, format = '%(message)s', filename = 'logs/run.log')

    # dump the sample files and locations into json
    with open (args.config, 'r') as f:
        fastqs = f.readlines().strip()
    samples = get_samples(fastqs)

    with open(f'{logs}/sample_list.json', 'w', encoding = 'utf-8') as log_file:
        json.dump(samples, log_file, indent = 4)

    # mp worker asignment
    with multiprocessing.Pool(processes = multiprocessing.cpu_count()) as pool:
        results = []
        for index, row in sample_manifest.iterrows():
            sample = row['Sample_ID']
            results.append(pool.apply_async(worker, args = (sample, indir)))

        # Wait for all processes to finish
        for result in results:
            result.get()


if __name__ == "__main__":
    main()
