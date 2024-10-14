# `star_alignment.py`

## Overview
The `star_alignment.py` script is designed to align RNA-seq sequences to a reference genome using the STAR aligner. 
It consists of three major components: sample parsing, trimming, and alignment. 
the script generates trims the samples using `trimmomatic`. It aligns using `STAR` aligner. Finally it generates `flagstat` using `samtools`

## Requirements

    Python 3.x
    Required Python packages: argparse, subprocess, multiprocessing, os, sys, logging, time, json
    External tools:
        STAR
        Trimmomatic
        Samtools

## Instructions

Install Dependencies: Ensure that all required `python` packages and external tools are installed and accessible in your `PATH`.
We loaded them as modules on a HPC cluster
### Prepare Configuration:
Create a configuration file listing the sample files to be analyzed. Each line should contain the path to a sample file.
This was generated via:
```shell
find "ABS_PATH_TO_DIR" -type f -name "*.gz" > ${RESULTS}/samplemanifest.txt
```
### Run the Script:
Use the command line to execute the script with the necessary arguments:
```shell
python star_alignment.py -c <config_file> -g <genome_directory> -o <output_directory>
```
Replace `<config_file>` with the path to your `samplemanifest.txt` file (this is generated as a list of found fastq in a directory), 
`<genome_directory>` with the path to your `STAR` indexed genome directory, and `<output_directory>` with the desired output directory path.

## Output
Upon completion there will be a directory with various files. In the `alignment_star` folder, all of the alignment files will be found.
In the logs file we will find statistics about the run. in the `stats` folder we will find the metrics from `samtools` and finally the trimmed `fq.gz` filles will be in `trimmed`
```shell
.
|-- alignment_star
|   |-- <sample_name>.*
|-- logs
|   |-- <run_id>.log
|-- stats
|   |-- <sample_name>.flagstat.tsv
|-- trimmed
|   |-- <sample_name>.R1.unpaired.fq.gz
|   |-- <sample_name>.R2.unpaired.fq.gz
|   |-- <sample_name>.R1.paired.fq.gz
|   |-- <sample_name>.R2.paired.fq.gz
```