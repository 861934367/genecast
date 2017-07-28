import pandas as pd
import os
from genecast_package.snv_analysis import *


def get_geneid_stand():
    pass


def save_parameters(args=None):
    f = open("parameters.txt", "w")
    for arg in dir(args):
        if not arg.startswith("_"):
            f.write(arg + ": " + str(getattr(args, arg)) + "\n")


def prepare_for_maf(data, file_name, host_gene=None, args=None):
    data = data[["Gene.refGene", "Chr", "Start", "End", "ExonicFunc", "Ref", "Alt", "ratio", "AAChange.refGene"]]
    tran, aa, Variant_Type = [], [], []
    for AAChange in data["AAChange.refGene"]:
        AAChange = AAChange.split(",")[0]
        tran.append(AAChange.split(":")[1])
        aa.append(AAChange.split(":")[-1])
    for ref, alt in zip(data["Ref"], data["Alt"]):
        if "-" in ref:
            Variant_Type.append("INS")
        elif "-" in alt:
            Variant_Type.append("DEL")
        else:
            Variant_Type.append("SNP")
    data["Chr"] = [i.strip("chr") for i in data["Chr"]]
    data["i_transcript_name"] = tran
    data["Protein_Change"] = aa
    data["Tumor_Sample_Barcode"] = [file_name] * len(aa)
    data["NCBI_Build"] = ["hg19"] * len(aa)
    data["Center"] = ["genecast"] * len(aa)
    data["Tumor_Seq_Allele1"] = data["Ref"]
    data["Variant_Type"] = Variant_Type
    data["Entrez_Gene_Id"] = [None] * len(aa)
    data["Strand"] = [None] * len(aa)
    data = data[["Gene.refGene", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chr", "Start", "End", "Strand", "ExonicFunc",
                 "Variant_Type", "Ref", "Tumor_Seq_Allele1", "Alt", "Tumor_Sample_Barcode", "Protein_Change", "ratio",
                 "i_transcript_name"]]
    data.columns = ["Hugo_Symbol", "Entrez_Gene_Id", "Center", "NCBI_Build", "Chromosome", "Start_Position", "End_position",
                    "Strand", "Variant_Classification", "Variant_Type", "Reference_Allele", "Tumor_Seq_Allele1",
                    "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "Protein_Change", "i_TumorVAF_WU", "i_transcript_name"]
    if host_gene:
        data = pd.merge(host_gene, data, left_on="gene", right_on="Hugo_Symbol", how="left")
        del data["gene"]
    return data


def _get_group_data(gene_data, g, group_name, args=None):
    group = {"Tumor_Sample_Barcode": [], "FAB_classification": []}
    if len(g) == 0:
        raise FileNoExist('no file get!!!, please check your file name!')
    data = []
    for i, file in enumerate(g):
        file_name = file.split("/")[-1].split(".")[0]
        group["Tumor_Sample_Barcode"].append(file_name)
        group["FAB_classification"].append(group_name)
        data1 = filter(merge_snp_indel(file, args=args), file_name, args=args)
        #data.append(prepare_for_maf(data1, file_name, host_gene=None, args=None))
        data.append(data1)
    data = pd.concat(data)
    group = pd.DataFrame(group)
    group = group[["Tumor_Sample_Barcode", "FAB_classification"]]
    return group, data


def maf_main(args=None):
    target = args.group1[0].split("/")[-2] + "_VS_" + args.group2[0].split("/")[-2] + "_maf"
    try:
        os.mkdir(target)
    except FileExistsError:
        sh.rm("-rf", target)
        os.mkdir(target)
    gene_list = get_host_gene(args=args)
    a_group, a_gene_data = _get_group_data(gene_list, args.group1, args.group1[0].split("/")[-2], args=args)
    b_group, b_gene_data = _get_group_data(gene_list, args.group2, args.group2[0].split("/")[-2], args=args)
    data = pd.concat([a_gene_data, b_gene_data])
    group = pd.concat([a_group, b_group])
    os.chdir(target)
    save_parameters(args=args)
    data = data.astype(str)
    data.to_csv("laml.txt", sep="\t", index=False)
    group.to_csv("laml_fab_annotation.txt", sep="\t", index=False)