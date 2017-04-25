import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import multiprocessing
import sh


def bin(group, n):
    depth = []
    group = np.array(group)
    for i in range(0, len(group), n):
        depth.append(np.median(group[i:i + n]))
    return depth


def multiprocessing_plot(file, args):
    re = sh.samtools("depth", file, "-b", "/home/zhout/data/panel6.bed")
    f = open(file.strip(".") + ".depth", "wb")
    f.write(re.stdout); f.close()
    fig, ax = plt.subplots()
    data = pd.read_table(file.strip(".") + ".depth", names=["chr", "pos", "depth"])
    ax.plot(bin(data["depth"], args.n), linewidth=0.5, label=file.split(".")[0])
    plt.ylabel("depth")
    plt.xlabel("panel_loction(bin=%d)" % args.n)
    plt.legend(loc="upper left", ncol=2)
    plt.title("Coverage Depth")
    plt.savefig("%s_depth_coverage.png" % file.split(".")[0], dpi=600)
    plt.close()
    fig, ax = plt.subplots()
    ax.plot(bin(data["depth"].sort_values(), args.n), linewidth=2, label=file.split(".")[0])
    plt.ylabel("depth")
    plt.xticks([])
    plt.xlabel("sort_depth(bin=%d)" % args.n)
    plt.legend(loc="upper left", ncol=2)
    plt.title("Uniformity")
    plt.savefig("%s_uniformity.png" % file.split(".")[0], dpi=600)
    plt.close()
    


def plot_coverage(args=None):
    pool = multiprocessing.Pool(processes=args.progress)
    for file in args.bams:
        pool.apply_async(multiprocessing_plot, (file, args))
    pool.close(); pool.join()
    
