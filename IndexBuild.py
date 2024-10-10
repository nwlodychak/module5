def build_index(genome, indir, outdir):
    """
    build a genomic index for any genome
    :param genome: raw genome of interest
    :param indir: where do your genome files live
    :param outdir: where do you want the index files to go?
    :return: index build assignment location.
    """
    assert isfile(genome), f"{genome} does not exist"
    index = os.path.basename(genome)
    try:
        logging.info(f"Building index for {genome}")
        command = (f"bowtie-build -p ${SBATCH_THREADS} {genome} {index}")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise