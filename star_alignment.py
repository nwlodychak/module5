"""
star_alignment.py
Created on 2024.10.09
Author: Nick Wlodychak
This module is used to align sequences to a reference genome genereated with STAR.
It consists of three major parts - Sample parsing, trimming and alignment.
When running this script we should produce an alignment for samples in the manifest and a stats file
for samtools alignment.
"""
import argparse
import subprocess
import multiprocessing
import os
import sys
import logging
import time
import json


def get_args():
    """
    Parse command line arguments
    :return: args
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-c',
                        '--config',
                        required = True,
                        help = 'Config file containing samples to analyze')
    parser.add_argument('-g',
                        '--genome',
                        required = True,
                        help = 'Raw Genome file location')
    parser.add_argument('-o',
                        '--outdir',
                        required = True,
                        help = 'Where to send all files.')
    return parser.parse_args()


def run_command(command):
    """
    Run a bash command
    :param command:
    """
    process = subprocess.Popen(command,
                               shell = True,
                               stdout = subprocess.PIPE,
                               stderr = subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command failed with error: {stderr.decode().strip()}")
    return stdout.decode().strip()


def trim(sample_id, read1, read2, outdir):
    """
    Standard adapter trimming function
    :param sample_id: sample basename
    :param read1: full path to read 1
    :param read2: full path to read 2
    :param outdir: where the trimmed fastq files will go
    :return: trimmed samples for R1 and R2
    """

    # trimming variables
    outtrim = f'{outdir}/trimmed'
    trim_p_ext1 = "_R1.paired.fq.gz"
    trim_p_ext2 = "_R2.paired.fq.gz"
    trim_u_ext1 = "_R1.unpaired.fq.gz"
    trim_u_ext2 = "_R2.unpaired.fq.gz"
    trim1 = os.path.join(outtrim, sample_id + trim_p_ext1)
    trim2 = os.path.join(outtrim, sample_id + trim_p_ext2)
    utrim1 = os.path.join(outtrim, sample_id + trim_u_ext1)
    utrim2 = os.path.join(outtrim, sample_id + trim_u_ext2)

    os.makedirs(outtrim, exist_ok = True)

    assert read1 is not None and read2 is not None, \
        (logging.error('Read 1 is %s\n', read1)) and logging.error('Read 2 is %s\n', read2)

    # trimming main command
    try:
        logging.info('Trimming: %s', sample_id)
        command = f"trimmomatic PE -threads {multiprocessing.cpu_count()} \
                       {read1} {read2} \
                       {trim1} {utrim1} \
                       {trim2} {utrim2} \
                       ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 \
                       SLIDINGWINDOW:4:20 MINLEN:25"
        run_command(command)
        logging.info('Trimming complete: %s', trim1)
        logging.info('Trimming complete: %s', trim2)
    except Exception as error:
        logging.error('An error occurred: %s', error)
        raise

    return trim1, trim2


def align_reads(genome, sample_id, fq1, fq2, outdir):
    """
    align RNA-seq reads using STAR or HISAT
    :param sample_id: basename for sample
    :param genome: path to indexed genome of interest
    :param fq1: trimmed RNASeq read1
    :param fq2: trimmed RNASeq read2
    :param outdir: where do you want the alignment files to go?
    """

    # alignment main command
    try:
        logging.info('Beginning alignment: %s', sample_id)
        command = f"STAR --genomeDir {genome} \
                          --readFilesCommand zcat \
                          --runThreadN {multiprocessing.cpu_count()} \
                          --readFilesIn {fq1} {fq2} \
                          --outFileNamePrefix {outdir}/alignment_star/{sample_id} \
                          --outSAMtype BAM SortedByCoordinate \
                          --outSAMunmapped Within \
                          --outSAMattributes Standard"
        run_command(command)
        logging.info('Alignment complete: %s', sample_id)
    except Exception as error:
        logging.error('An error occurred during alignment: %s', error)


def samtools_alignment(sample_id, alignment, outdir):
    """
    Runs samtools - based on different outputs for the rna-seq alignment, we run different commands
    :param sample_id: basename of the sample
    :param alignment: BAM or SAM file, depending on the output of the aligner
    :param outdir: where do the final files go?
    :return:
    """
    os.makedirs(f'{outdir}/stats', exist_ok = True)

    # samtools main command
    try:
        command = f"samtools index {alignment}"
        run_command(command)
        command = f"samtools flagstat -O tsv {alignment} > {outdir}/stats/{sample_id}.stats.tsv"
        run_command(command)
        logging.info('Alignment complete: %s', sample_id)

    except Exception as error:
        logging.error('Samtools alignment error: %s', error)


def worker(sample_id, read1, read2, outdir):
    """
    worker function for multiprocessing aligner
    :param sample_id: base name sample id
    :param read1: read1 path location
    :param read2: read2 path location
    :param outdir: where do you want to save the output files?
    :return:
    """

    genome = os.path.normpath(get_args().genome)
    logging.info('Genome at: %s', genome)

    try:
        logging.info('Trimming: %s', sample_id)
        trim1, trim2 = trim(sample_id, read1, read2, outdir)
        time.sleep(5)
    except Exception as error:
        logging.error('Exception occurred during TRIM: %s', error)

    # alignment block
    try:
        logging.info('Aligning STAR: %s', sample_id)
        align_reads(genome, sample_id, fq1 = trim1, fq2 = trim2, outdir = outdir)
        time.sleep(5)
    except Exception as error:
        logging.error('Exception occurred during STAR: %s', error)

    # samtools block
    try:
        logging.info('Samtools: %s', sample_id)
        os.makedirs(f'{outdir}/alignment_star', exist_ok = True)
        samtools_alignment(sample_id,
                           alignment = f'{outdir}/'
                                       f'alignment_star/'
                                       f'{sample_id}Aligned.sortedByCoord.out.bam',
                           outdir = outdir)
        time.sleep(5)
    except Exception as error:
        logging.error('Exception occurred during STAR: %s', error)

    # logging.info(f'Removing preprocessing files')
    # remove_files = [trim1, trim2]
    # for file in remove_files:
    #     logging.info(f'Removing {file}')
    #     os.remove(file)


def get_samples(files):
    """
    pairs up r1 and r2 samples from a sample manifest into a dict
    :param files: list of file locations
    :return: dict of sample R1 and R2 locations
    """
    sample_dict = {}

    # sample parser from file names
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
    """
    Main function for the script, map to workers set variables
    :return:
    """
    # environ variables
    cpus = int(os.environ.get('SLURM_CPUS_PER_TASK',
                              multiprocessing.cpu_count()))
    job_id = os.environ.get('SLURM_JOB_ID', 'default_value')
    args = get_args()

    # logging assignment
    outdir = os.path.normpath(args.outdir)
    os.makedirs(outdir, exist_ok = True)
    logs = f'{args.outdir}/logs'
    os.makedirs(logs, exist_ok = True)
    logging.basicConfig(handlers = [logging.FileHandler(filename = f'{logs}/{job_id}_run.log',
                                                        encoding = 'utf-8', ),
                                    logging.StreamHandler(sys.stdout)],
                        format = '%(levelname)-8s %(message)s',
                        level = logging.INFO)

    # dump the sample files and locations into json
    with open(args.config, 'r', encoding = 'utf-8') as file:
        fqs = [line.strip() for line in file.readlines()]

    samples = get_samples(fqs)

    with open(f'{logs}/sample_list.json', 'w', encoding = 'utf-8') as log_file:
        json.dump(samples, log_file, indent = 4)

    # mp worker assignment
    logging.info('Assigning workers for job - %d.', job_id)
    logging.info('Running on %d cores.', cpus)
    logging.info('Samples accounted for %d.', len(samples.keys()))
    with multiprocessing.Pool(processes = cpus) as pool:
        results = []
        for sample_id, paths in samples.items():
            read1 = paths['R1']
            read2 = paths['R2']
            assert os.path.exists(read1), f'{read1} does not exist'  # ensure the files exist
            assert os.path.exists(read2), f'{read2} does not exist'
            results.append(pool.apply_async(worker,
                                            args = (sample_id, read1, read2, outdir)))

        # wait for all processes to finish
        for result in results:
            result.get()


if __name__ == "__main__":
    main()
