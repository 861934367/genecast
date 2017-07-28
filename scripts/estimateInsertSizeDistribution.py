import matplotlib as mpl
mpl.use('Agg')
import pysam
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import argparse as apa
import multiprocessing


def isPaired(bamfile, alignments=1000):
    '''check if a *bamfile* contains paired end data

    The method reads at most the first *alignments* and returns
    True if any of the alignments are paired.
    '''
    samfile = pysam.Samfile(bamfile)
    n = 0
    for read in samfile:
        if read.is_paired:
            break
        n += 1
        if n == alignments:
            break
    samfile.close()
    return n != alignments

def estimateInsertSizeDistribution(bamfile, alignments=1000):
    '''
    estimate insert size from first alignments in bam file.
    returns mean and stddev of insert sizes.
    '''
    assert isPaired(bamfile), \
        'can only estimate insert size from' \
        'paired bam files'
    samfile = pysam.Samfile(bamfile)
    # only get positive to avoid double counting
    inserts = np.array([read.tlen for read in samfile if read.is_proper_pair and read.tlen > 0])
    return inserts


def plot(list, name):
    sns.distplot(list)
    plt.savefig(name + "_insert_size.png")
    plt.close()


def plot_insert_size(bam_file):
    inserts = estimateInsertSizeDistribution(bam_file)
    plot(inserts, bam_file.rstrip(".bam"))


if __name__ == "__main__":
    parser = apa.ArgumentParser(prog="convert")
    parser.add_argument("bams", nargs='*', help="*.bam")
    parser.add_argument('-p', '--progress', required=False, default=1, type=int,
                        help='parallel')
    args = parser.parse_args()
    pool = multiprocessing.Pool(processes=args.progress)
    for file in args.bams:
        pool.apply_async(plot_insert_size, (file,))
    pool.close()
    pool.join()