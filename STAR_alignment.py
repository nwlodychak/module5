import argparse
import subprocess
import multiprocessing
import os
import logging
import time
import json


def get_args():
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
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command failed with error: {stderr.decode().strip()}")
    return stdout.decode().strip()


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


def trim(sample_id, read1, read2, outdir):
    """
    Standard adapter trimming function
    :param sample_id: sample basename
    :param read1: full path to read 1
    :param read2: full path to read 2
    :param outdir: where the trimmed fastq files will go
    :return: trimmed samples for R1 and R2
    """
    outtrim = f'{outdir}/trimmed'
    trim_ext1 = "_R1.trimmed.fastq"
    trim_ext2 = "_R2.trimmed.fastq"
    trim1 = os.path.join(outtrim, sample_id + trim_ext1)
    trim2 = os.path.join(outtrim, sample_id + trim_ext2)

    os.makedirs(outtrim, exist_ok = True)

    assert read1 is not None and read2 is not None, \
        (logging.error(f"Read 1 is {read1}\n")) and logging.error(f"Read 2 is {read2}\n")

    try:
        logging.info(f"Trimming {sample_id}")
        command = (f"trimmomatic PE \
                    {read1} {read2} \
                    {trim1} {trim2} \
                    SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15")
        run_command(command)
        logging.info(f"Trimming complete - {trim1}")
        logging.info(f"Trimming complete - {trim2}")
    except Exception as e:
        logging.error(f"An error occurred: {e}")
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

    try:
        logging.info(f"Beginning alignment {sample_id}")
        command = (f"STAR --genomeDir {genome} \
                          --readFilesCommand zcat \
                          --runThreadN {multiprocessing.cpu_count()} \
                          --readFilesIn {fq1} {fq2} \
                          --outFileNamePrefix {outdir}/alignment_star/{sample_id} \
                          --outSAMtype BAM SortedByCoordinate \
                          --outSAMunmapped Within \
                          --outSAMattributes Standard")
        run_command(command)
        logging.info(f"Alignment complete - {sample_id}")
    except Exception as e:
        logging.error(f"An error occurred during alignment: {e}")


def samtools_alignment(sample_id, alignment, outdir):
    """
    Runs samtools - based on different outputs for the rna-seq alignment, we run different commands
    :param sample_id: basename of the sample
    :param alignment: BAM or SAM file, depending on the output of the aligner
    :param outdir: where do the final files go?
    :return:
    """
    try:
        command = f"samtools index {alignment}"
        run_command(command)
        command = f"samtools flagstat -O tsv {alignment} > {outdir}/stats/{sample_id}.stats.tsv"
        run_command(command)
        logging.info(f"Alignment complete {sample_id}")

    except Exception as e:
        logging.error(f"Samtools alignment error: {e}")


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
    logging.info(f"Genome at - {genome}")

    # fastqc block
    # try:
    #     logging.info(f'FASTQC {sample_id}')
    #     fastqc(read1, outdir)
    #     fastqc(read2, outdir)
    #     time.sleep(5)
    # except Exception as e:
    #    logging.error(f"Exception occurred during FASTQC: {e}")

    # trimming block
    try:
        logging.info(f'Trimming {sample_id}')
        trim1, trim2 = trim(sample_id, read1, read2, outdir)
        time.sleep(5)
        logging.info(f'Trimming - {trim1}')
        logging.info(f'Trimming - {trim2}')
    except Exception as e:
        logging.error(f"Exception occurred during TRIM: {e}")

    # alignment block
    try:
        logging.info(f'Aligning STAR - {sample_id}')
        align_reads(genome, sample_id, fq1 = trim1, fq2 = trim2, outdir = outdir)
        time.sleep(5)
    except Exception as e:
        logging.error(f"Exception occurred during STAR: {e}")

    # samtools block
    try:
        logging.info(f'Samtools {sample_id}')
        os.makedirs(f'{outdir}/alignment_star', exist_ok = True)
        samtools_alignment(sample_id,
                           alignment = f'{outdir}/alignment_star/{sample_id}Aligned.sortedByCoord.out.bam',
                           outdir = outdir)
        time.sleep(5)
    except Exception as e:
        logging.error(f"Exception occurred during STAR: {e}")

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
    logging.basicConfig(level = logging.INFO, format = '%(message)s', filename = f'{logs}/{job_id}_run.log')

    # dump the sample files and locations into json
    with open(args.config, 'r', encoding = 'utf-8') as f:
        fqs = [line.strip() for line in f.readlines()]

    samples = get_samples(fqs)

    with open(f'{logs}/sample_list.json', 'w', encoding = 'utf-8') as log_file:
        json.dump(samples, log_file, indent = 4)

    # mp worker assignment
    logging.info(f'Assigning workers for job {job_id}.')
    logging.info(f'Running on {cpus} cores.')
    logging.info(f'Samples accounted for - {len(samples.keys())}.')
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

    # final quality assessment
    multiqc(outdir)


if __name__ == "__main__":
    main()
