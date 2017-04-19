## this tool is the core function of cnv and snv analysis
## author: taozhou
## email: zhou.tao@genecast.com.cn

import matplotlib as mpl
mpl.use('Agg')
import warnings
warnings.filterwarnings("ignore")
import itertools
import seaborn as sns
import matplotlib.pylab as plt
import matplotlib.colors as mc
from genecast_package.svm_analysis import feature_select, evaluate_model
from sklearn.decomposition import PCA
from collections import OrderedDict
from collections import defaultdict
import datetime
import pandas as pd
import os
import sh
import warnings
warnings.filterwarnings("ignore")


def z_score(data, axis):
    if axis == 1:
        z_scored = data
    else:
        z_scored = data.T
    z_scored = (z_scored - z_scored.mean()) / z_scored.std()

    if axis == 1:
        return z_scored
    else:
        return z_scored.T


def pheatmap(data, length, col_cluster=True, xticklabels=True, yticklabels=True, save="pdf", color=None, name=None):
    data = z_score(data, axis=0)
    if len(data.columns) > 30:
        xticklabels = False
    if len(data) > 80:
        yticklabels = False
    vmin, vmax = data.unstack().quantile([.01, .99])
    re = sns.clustermap(data, cmap="bwr", row_cluster=True, col_cluster=col_cluster, figsize=(13, 10), \
                        xticklabels=True, yticklabels=yticklabels, vmin=vmin, vmax=vmax, col_colors=color)
    re.ax_heatmap.set_xticklabels(re.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
    re.ax_heatmap.set_yticklabels(re.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    if col_cluster == False:
        for group, number in length.items():
            re.ax_col_colors.text((number[0] + number[1])/2 - len(group)/2, 1.1, group, size=30)
        re.savefig(name + "." + save)
    else:
        re.savefig(name + "_col_cluster." + save)
    plt.close()


def make_col_color_heatmap(group_dic):
    common_color = ["blue", "red", "green", "grey"]
    color = {}; length = {}
    temp = 0
    i = 0
    for name, group in group_dic.items():
        length[name] = [temp, temp + len(group)]
        temp += len(group)
        for sample in group:
            color[sample] = common_color[i]
        i += 1
    color = pd.Series(color)
    color.name = "group"
    return color, length


def pca(data, group_dic, n=None):
    pca = PCA(n_components=2)
    group = []
    length = OrderedDict()
    temp = 0
    for name, g in group_dic.items():
        length[name] = [temp, temp + len(g)]
        temp += len(g)
        group += g
    data = data[group]
    newData = pca.fit_transform(data.T)
    colors = ["blue", "red", "green", 'turquoise', "grey"]
    i = 0
    for name, number in length.items():
        plt.scatter(newData[number[0]:number[1], 0], newData[number[0]:number[1], 1], label=name, color=colors[i])
        i += 1
    plt.title("PCA analysis")
    pc1 = 100*pca.explained_variance_ratio_[0]
    pc2 = 100*pca.explained_variance_ratio_[1]
    plt.xlabel("PC1(%.1f)" % pc1)
    plt.ylabel("PC1(%.1f)" % pc2)
    plt.legend()
    plt.savefig("PCA_%s.png" % n)
    plt.close()


def plot_box(data, which, outname, palette, regulation, group):
    fig, ax1 = plt.subplots(figsize=(8,12))
    box_data = defaultdict(list)
    if which == "cnv":
        how = "mean"
        for name, g in group.items():
            box_data[name] = data[g].mean(0)
    else:
        how = "sum"
        for name, g in group.items():
            box_data[name] = data[g].sum(0)
    data.to_csv(outname + "_box_data_%s_%s" % (regulation, how) + ".txt", sep="\t")
    sns.boxplot(data=pd.DataFrame(box_data), ax=ax1, width=0.2, linewidth=.5, palette=palette)
    ax1.set_title(outname)
    ax1.set_ylabel('%s value(%s)' % (which, how))
    fig.autofmt_xdate(ha='center', rotation=0)
    fig.savefig(r'%s_box_data_%s_%s_Boxplot.png' % (outname, regulation, how), dpi=600, size=0.5)
    plt.close()


def databox(raw, which, outname=None, group=None):
    palette = {}
    up = []; down = []
    group1_data = raw[list(group.values())[0]]
    group2_data = raw[list(group.values())[1]]
    color = ["red", "green", "blue"]
    for gene in raw.index:
        if group1_data.ix[gene].sum() - group2_data.ix[gene].sum() >= 0:
            up.append(gene)
        else:
            down.append(gene)
    for i, (name, g) in enumerate(group.items()):
        palette[name] = color[i]
    plot_box(raw.ix[up], which, outname, palette, "up", group)
    plot_box(raw.ix[down], which, outname, palette, "down", group)


def save_data_pdf(data, name, length, color, group_dic, which):
    data.to_csv("%s.txt" % name, sep="\t")
    length = {key.split("/")[-1]: value for key, value in length.items()}
    group_dic = {key.split("/")[-1]: value for key, value in group_dic.items()}
    pheatmap(data, length, col_cluster=True, color=color, name=name, save="png")
    pheatmap(data, length, col_cluster=False, color=color, name=name, save="png")
    pca(data, group_dic, n=name)
    databox(data, which, outname=name, group=group_dic)


def save_parameters(args=None, which="cnv"):
    pass


def make_result_folder(args=None, which="cnv", fun=None):
    feature_genes = []; gene_lists = {}; color_length = {}
    os.chdir(args.outdir)
    i = datetime.datetime.now()
    for two_group in itertools.combinations([args.group1, args.group2], 2):
        target = two_group[0].split("/")[-1] + "_VS_" + two_group[1].split("/")[-1] + "_%s%s%s_%s%s" % (i.year, i.month, i.day, i.hour, i.minute)
        try:
            os.mkdir(target)
        except FileExistsError:
            sh.rm("-rf",target)
            os.mkdir(target)
        if which == "cnv":
            name = "cnv_median_" + args.data_type
            gene_list, a_group, b_group = fun(args.host_gene, two_group[0], two_group[1], data_type=args.data_type)
        else:
            if args.cal_type == "num":
                name = "snv_number"
            else:
                name = "snv_mean"
            gene_list, a_group, b_group = fun(args.host_gene, two_group[0], two_group[1], args.cal_type, which)
        feature_gene = feature_select(gene_list, a_group, b_group, pval=args.pval, method=args.feature_selection_method,\
                                      criterion=args.criterion, penalty=args.penalty, C=args.C, threshold=args.threshold)
        feature_genes.append(feature_gene)
        gene_lists[two_group[0]] = gene_list[a_group]; gene_lists[two_group[1]] = gene_list[b_group]
        os.chdir(target)
        save_parameters(args=args, which=which)
        group_dic = {two_group[0]: a_group, two_group[1]: b_group}
        color_length[two_group[0]] = a_group; color_length[two_group[1]] = b_group
        color, length = make_col_color_heatmap(group_dic)
        save_data_pdf(gene_list, "host_gene_%s" % name, length, color, group_dic, which)
        pd.DataFrame({"gene":feature_gene}).to_csv("feature_gene_pval%0.2f.txt" % args.pval, sep="\t", index=False)
        feature_gene_cnv = gene_list.ix[feature_gene]
        evaluate_model(gene_list, a_group, b_group, feature_gene, name="feature_gene_%s" % name, method=args.prediction_method, C=args.C, n_folds=args.n_folds)
        save_data_pdf(feature_gene_cnv, "feature_gene_%s" % name, length, color, group_dic, which)
        os.chdir(args.outdir)
    if len([args.group1, args.group2]) > 2:
        try:
            os.mkdir("intersection")
        except FileExistsError:
            pass
        os.chdir("intersection")
        color, length = make_col_color_heatmap(color_length)
        intersection_feature_gene = list(set(feature_genes[0]).intersection(*feature_genes[1:]))
        intersection_feature_gene_cnv = pd.concat([data.ix[intersection_feature_gene] for [args.group1, args.group2], data in gene_lists.items()], axis=1)
        try:
            save_data_pdf(intersection_feature_gene_cnv, "intersection", length, color, color_length)
        except Exception:
            print("no intersection\njob finish...")
        os.chdir(args.outdir)