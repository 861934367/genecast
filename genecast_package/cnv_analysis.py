## this tool is for cnv analysis
## author: taozhou
## email: zhou.tao@genecast.com.cn

import pandas as pd
from glob import glob
import numpy as np
import os
from genecast_package.core import make_result_folder
from genecast_package.snv_analysis import get_host_gene
import warnings
warnings.filterwarnings("ignore")


class MethodException(Exception):
    pass


def split_gene_data(data, data_type):
    new_data = {"gene": [], data_type: []}
    for genes, value in zip(data["gene"], data[data_type]):
        for gene in genes.split(";"):
             new_data["gene"].append(gene)
             new_data[data_type].append(value)
    data = pd.DataFrame(new_data)
    return data


def parser_cnr(file, args=None):
    data = pd.read_table(file, usecols=["gene", "log2"])
    data = data.loc[data["gene"] != "Background"]
    data = split_gene_data(data, args.data_type)
    groups = pd.DataFrame(data.groupby(data["gene"]).median())
    groups.columns = [file.split("/")[-1].split(".")[0]]
    return groups


def parser_call_cnr(file, args=None):
    data = pd.read_table(file, usecols=["gene", "log2", "cn"])
    data = data.loc[data["gene"] != "Background"]
    data = split_gene_data(data, args.data_type)
    groups = pd.DataFrame(data[args.data_type].groupby(data["gene"]).median())
    groups.columns = [file.split("/")[-1].split(".")[0]]
    return groups


def get_host_gene_cnv(args=None):
    gene_list = get_host_gene(args=args)
    if args.data_type == "log2": fun = parser_cnr; fillna_num = 0
    else: fun = parser_call_cnr; fillna_num = 2
    a_group = []
    for file in args.group1:
        a_group.append(file.split("/")[-1].split(".")[0])
        gene_list = pd.merge(gene_list, fun(file, args=args), left_on="gene", right_index=True, how="left").fillna(fillna_num)
    b_group = []
    for file in args.group2:
        b_group.append(file.split("/")[-1].split(".")[0])
        gene_list = pd.merge(gene_list, fun(file, args=args), left_on="gene", right_index=True, how="left").fillna(fillna_num)
    gene_list.index = gene_list["gene"]
    del gene_list["gene"]
    # if 0 in gene_list.dropna(how="all").fillna(0):
        # data = gene_list.dropna(how="all").fillna(0).drop(0, axis=0)
    # else:
        # data = gene_list.dropna(how="all").fillna(0)
    if args.data_type == "log2": data = gene_list.loc[~(gene_list.T==0).all()]
    else: data = gene_list.loc[~(gene_list.T==2).all()]
    return data, a_group, b_group


def cnv(args=None):
    make_result_folder(args=args, fun=get_host_gene_cnv, which="cnv")


if __name__ == "__main__":
    host_gene_file = "target_gene.txt"
    groups = ["CESC", "OV", "UCEC"]
    p = 0.05
    root_dir = os.getcwd()
    make_result_folder(host_gene_file, groups, p, root_dir, fun=get_host_gene_cnv, which="cnv", \
                       prediction_method="LinearSVC", C=1, n_folds=5, criterion='aic', penalty="l2", alpha=0.025, threshold=0)






