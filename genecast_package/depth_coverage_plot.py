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
    if args.type == "base": reads = bin(data[args.type], args.n)
    else: reads = np.log10(data[args.type])
    ax1, ax2 = axs[0], axs[1]
    ax1.bar(np.arange(len(reads)), reads, 1, color="slateblue")
    ax1.set_ylabel("$Log10(%s)$" % args.type, size=20)
    
    #ax1.set_title("Uniformity of Coverage (Average Coverage = %d)" % (average_depth), size=20)
    reads.sort()
    ax2.bar(np.arange(len(reads)), reads, 1, color="slateblue")
    ax2.set_ylabel("$Log10(%s)$" % args.type, size=20)
    ax2.set_xticks([])
    if args.type == "base":
        ax1.set_xlabel("panel_loction(bin=%d)" % args.n, size=20)
        ax2.set_xlabel("sort_depth(bin=%d)" % args.n, size=20)
        ax1.set_title("Uniformity of Coverage (Average Coverage = %d percentage20 = %0.3f)" % (average_depth, percentage20), size=20)
    else:
        ax1.set_xlabel("panel_loction(bin=panel_region)", size=20)
        ax2.set_xlabel("sort_depth(bin=panel_region)", size=20)
        ax1.set_title("Uniformity of Coverage (Average reads of panel_region = %d percentage20 = %0.3f)" % (average_depth, percentage20), size=20)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.9, top=0.9, hspace = 0.2, wspace = 0.3)
    plt.savefig("%s_uniformity_coverage.png" % file.split(".")[0], dpi=600)
    plt.close()


def multiprocessing_plot(file, args):
    if args.type == "reads":
        sam = pysam.AlignmentFile(file)
        data = pd.read_table("/home/zhout/data/%s.bed" % args.panel, names=["chr", "start", "end", "gene", "transcript"])
        data["reads"] = [sam.count(chr, start, end) / (end - start) for chr, start, end in zip(data["chr"], data["start"], data["end"])]
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
    ax1.set_title("Uniformity of Coverage")
    ax1.set_ylabel("$Log10(%s)$" % args.type)
    fig.autofmt_xdate(ha='center', rotation=0)
    plt.xticks(rotation=90)
    fig.savefig(r'%s_Uniformity_Boxplot.png' % (args.out), dpi=600)
    plt.close()

