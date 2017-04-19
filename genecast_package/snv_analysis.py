## this tool is for snv analysis
## author: taozhou
## email: zhou.tao@genecast.com.cn

import pandas as pd
from glob import glob
import numpy as np
import os
from genecast_package.core import make_result_folder


class MethodException(Exception):
    pass


def replace(file, indel):
    base_name = file.split("/")[-1].replace("snp", indel)
    return "/".join(file.split("/")[:-1]) + "/" + base_name


def merge_snp_indel(file, which):
    title = ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.refGene', 'Gene.refGene',
             'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene',
             'ljb2_pp2hdiv', 'ljb2_pp2hvar', 'ExAC_ALL', 'ExAC_AFR', 'ExAC_AMR',
             'ExAC_EAS', 'ExAC_FIN', 'ExAC_NFE', 'ExAC_OTH', 'ExAC_SAS',"CLINSIG","CLNDBN","CLNACC","CLNDSDB","CLNDSDBID",
             "gnomAD_exome_ALL", "gnomAD_exome_AFR", "gnomAD_exome_AMR", "gnomAD_exome_ASJ", "gnomAD_exome_EAS",
             "gnomAD_exome_FIN", "gnomAD_exome_NFE", "gnomAD_exome_OTH", "gnomAD_exome_SAS",
             'Otherinfo', ".", "_", "_", "_", "_", "_", "_", "_", "_", "_", "_", "_", "ratio"]
    if which == "snp":
        return pd.read_table(file, skiprows=1, names=title)
    elif which == "indel":
        indel_file = replace(file, "indel")
        pd.read_table(indel_file, skiprows=1, names=title)
    else:
        data_snp = pd.read_table(file, skiprows=1, names=title)
        indel_file = replace(file, "indel")
        data_indel = pd.read_table(indel_file, skiprows=1, names=title)
        data_snp_indel = pd.concat([data_snp, data_indel])
        #data_snp_indel.to_csv(replace(file, "snp_indel"), index=False, sep="\t")
        return data_snp_indel


def filter(data, file_name, cal, circos=False):
    target_ExonicFunc = ["frameshift", "nonframeshift", "nonsynonymous", "stopgain", "stoploss"]
    data = data.loc[data["Func.refGene"] == "exonic"]
    data["ExonicFunc"] = [i.split(" ")[0] for i in data["ExonicFunc.refGene"]]
    data = data.loc[data["ExonicFunc"].isin(target_ExonicFunc)]
    B = []
    for ljb2_pp2hdiv, ljb2_pp2hvar in zip(data["ljb2_pp2hdiv"], data["ljb2_pp2hvar"]):
        if "B" in ljb2_pp2hdiv and "B" in ljb2_pp2hvar:
            B.append("B")
        else:
            B.append("ok")
    data["B"] = B
    data = data.loc[data["B"] != "B"]
    ExAC_columns = [i for i in data.columns if "gnomAD" in i]
    data["gnomAD_max"] = data[ExAC_columns].max(1)
    ratio = []
    for i in data["ratio"]:
        try:
            ratio.append(i.split(":")[5])
        except AttributeError:
            ratio.append(None)
    data["ratio"] = ratio
    data["max"] = [False if i != "." and float(i) >= 0.001 else True for i in data["gnomAD_max"]]
    data["ratio"] = [float(i.rstrip("%")) for i in data["ratio"]]
    if circos:
        data = data[['Chr', 'Start', 'End', 'Gene.refGene', "ratio"]]
        return data
    else:
        data = data.loc[data["max"] == True][["Gene.refGene", "ratio"]]
        groups = data.groupby(data["Gene.refGene"])
        data = pd.merge(groups.count(), groups.mean(), left_index=True, right_index=True, how="inner")
        data.columns = ["num", "mean"]
        data = pd.DataFrame({file_name:data[cal]}, index=data.index)
        return data


def get_host_gene(host_gene_file):
    if len(pd.read_table(host_gene_file).T) == 1:
        try:
            gene_list = pd.read_table(host_gene_file, usecols=["gene"]).drop_duplicates()
        except Exception:
            raise MethodException('the format of your target file is wrong, '
                              'please make sure it only contain one colomn and its title must be gene,'
                              'or panal bed file also be ok')
    else:
        try:
            gene_list = pd.read_table(host_gene_file, names=["chr", "start", "end", "gene", "trans"])[["gene"]].drop_duplicates()
        except Exception:
            raise MethodException('the format of your target file is wrong, '
                              'please make sure it only contain one colomn and its title must be gene,'
                              'or panal bed file also be ok')
    return gene_list


def _get_group_data(gene_data, g, cal, which, circos=False):
    group = []
    for file in glob(g + "/*snp*hg19_multianno.txt"):
        file_name = file.split("/")[-1].split(".")[0]
        group.append(file_name)
        if circos:
            gene_data = pd.merge(gene_data, filter(merge_snp_indel(file, which), file_name, cal, circos=circos), on=["chr", "start", "end", "gene"], how="outer")
        else:
            gene_data = pd.merge(gene_data, filter(merge_snp_indel(file, which), file_name, cal), left_on="gene", right_index=True, how="left")
    return group, gene_data


def get_host_gene_snv(host_gene_file, a, b, cal, which):
    gene_list = get_host_gene(host_gene_file)
    a_group, a_gene_data = _get_group_data(gene_list, a, cal, which)
    b_group, b_gene_data = _get_group_data(gene_list, b, cal, which)
    gene_data = pd.merge(a_gene_data, b_gene_data, on="gene", how="outer")
    gene_data.index = gene_data["gene"]
    del gene_data["gene"]
    if 0 in gene_data.dropna(how="all").fillna(0):
        data = gene_data.dropna(how="all").fillna(0).drop(0, axis=0)
    else:
        data = gene_data.dropna(how="all").fillna(0)
    return data, a_group, b_group
#    return gene_data.dropna(how="all").fillna(0).drop(0, axis=0), a_group, b_group


def make_karyotype(gene_list, unit):
    f = open("karyotype_hg19.txt", "w")
    for i, gene in enumerate(gene_list):
        f.write(
            "chr - " + gene + " " + gene + " " + str(0) + " " + str(unit) + " " + "chr" + str(i + 1) + "\n")


def fill_0(gene_group, unit):
    data = {"gene": [], "start": [], "end": [], "ratio": []}
    temp = 0
    for gene, start, ratio in zip(gene_group["gene"], gene_group["start"], gene_group["ratio"]):
        data["gene"].append(gene)
        if start != temp:
            data["start"].append(temp)
            data["end"].append(start)
            data["ratio"].append(0)
        data["start"].append(start)
        data["end"].append(start + 1)
        data["ratio"].append(ratio)
        temp = start + 1
    data["start"].append(temp)
    data["end"].append(unit)
    data["ratio"].append(0)
    return pd.DataFrame(data)


def circos_data(host_gene_file, root_dir, gene_location, groups, which="snv", cal="mean", unit=100000):
    try:
        os.mkdir("circos")
    except FileExistsError:
        pass
    gene_list = get_host_gene(host_gene_file)
    gene_location = pd.read_table(gene_location, names=["Chr", "Start", "End", "Gene.refGene"])
    gene_list = gene_location.loc[gene_location["Gene.refGene"].isin(gene_list)]
    norm_loc = {}
    for gene, start, end in zip(gene_list["Gene.refGene"], gene_list["Start"],
                                gene_list["End"]):
        norm_loc[gene] = [gene, start, (end - start) / unit]
    for g in groups:
        group, data = _get_group_data(gene_list, g, cal, which, circos=True)
        ratio = [i for i in data.columns if "ratio" in i]
        print(data)
        last_data = data[["gene", "start", "end"]]
        last_data["ratio"] = data[ratio].fillna(0).sum(1)
        last_group = []
        gene_groups = last_data.groupby(last_data["gene"])
        for gene, gene_group in gene_groups:
            gene_group["start"] = np.round((gene_group["start"] - norm_loc[gene][1]) / norm_loc[gene][2])
            gene_group = gene_group.sort_values("start")
            last_group.append(fill_0(gene_group, unit))
        last_group = pd.concat(last_group, axis=0)
        last_group.to_csv(g.split("/")[-1] + ".txt", sep="\t", index=False, header=None)
    make_karyotype(gene_list, unit)


def snv(args=None):
    make_result_folder(args=args, fun=get_host_gene_snv, which=args.data_type)


if __name__ == "__main__":
    host_gene_file = "panel6.bed"
    gene_location = "USCS_gene_location.txt"
    groups = ["snp_A_1_cf", "snp_B_1_cf"]
    p = 0.05
    root_dir = os.getcwd()
    #circos_data(host_gene_file, root_dir, gene_location, groups, which="snp", cal="mean", unit=100000)
    ## cal: mean num
    ## which: snv snp or indel
    make_result_folder(host_gene_file, groups, p, root_dir, which="snp", cal="mean", fun=get_host_gene_snv)