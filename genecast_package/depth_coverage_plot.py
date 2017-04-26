import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import multiprocessing
import sh
import pysam
from collections import defaultdict


class TypeException(Exception):
    pass


def bin(group, n):
    depth = []
    group = np.array(group)
    for i in range(0, len(group), n):
        depth.append(np.median(group[i:i + n]))
    return np.log10(depth)


def plot(data, args=None, file=None):
    fig, axs = plt.subplots(nrows=2, ncols=1,figsize=(15,12))
    average_depth = data[args.type].mean()
    percentage20 = len(data.loc[data[args.type] > average_depth * 0.2]) / len(data)
    reads = bin(data[args.type], 150)
    ax1 = axs[0]
    ax1.bar(np.arange(len(reads)), reads, 1, color="slateblue")
    ax1.set_ylabel("$Log10(%s)$" % args.type, size=20)
    ax1.set_xlabel("panel_loction(bin=%d)" % 150, size=20)
    if average_depth >= args.depth:
        ax1.set_title("Uniformity of Coverage (Average Coverage = %d percentage20 = %0.3f)" % (average_depth, percentage20), size=20)
    else:
        ax1.set_title("Uniformity of Coverage (Average Coverage = %d)" % (average_depth), size=20)
    ax2 = axs[1]
    reads.sort()
    ax2.bar(np.arange(len(reads)), reads, 1, color="slateblue")
    ax2.set_ylabel("$Log10(%s)$" % args.type, size=20)
    ax2.set_xticks([])
    ax2.set_xlabel("sort_depth(bin=%d)" % 150, size=20)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace = 0.2, wspace = 0.3)
    plt.savefig("%s_uniformity_coverage.png" % file.split(".")[0], dpi=600)
    plt.close()


def multiprocessing_plot(file, args):
    if args.type == "reads":
        sam = pysam.AlignmentFile(file)
        data = pd.read_table("/home/zhout/data/%s_for_qc.bed" % args.panel)
        data["reads"] = [sam.count(chr, pos, pos + 1) for chr, pos in zip(data["chr"], data["pos"])]
    elif args.type == "base":
        re = sh.samtools("depth", file, "-b", "/home/zhout/data/%s.bed" % args.panel)
        f = open(file.strip(".") + ".depth", "wb")
        f.write(re.stdout); f.close()
        data = pd.read_table(file.strip(".") + ".depth", names=["chr", "pos", "base"])
    else:
        raise TypeException("data type is wrong")
    plot(data, args=args, file=file)
    return np.log10(data[args.type])


def plot_coverage(args=None):
    pool = multiprocessing.Pool(processes=args.progress)
    box_data = {}
    for file in args.bams:
        box_data[file.split(".")[0]] = pool.apply_async(multiprocessing_plot, (file, args))
    pool.close(); pool.join()
    box_data = {key: value.get() for key, value in box_data.items()}
    data = pd.DataFrame(box_data)
    fig, ax1 = plt.subplots(figsize=(len(args.bams), 12))
    sns.boxplot(data=data, ax=ax1, width=0.2, linewidth=.5)
    ax1.set_title("Uniformity of overage")
    ax1.set_ylabel("$Log10(%s)$" % args.type)
    ax1.set_xticklabels(ax1.xaxis.get_majorticklabels(), rotation=90)
    fig.autofmt_xdate(ha='center', rotation=0)
    fig.savefig(r'%s_Boxplot.png' % ("Uniformity"), dpi=600)
    plt.close()

