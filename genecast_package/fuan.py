## use gff anno the fusion result of breakdancer
import pandas as pd
from glob import glob
import re
from collections import defaultdict
import multiprocessing


def find_gene(pos, data):
    target = set()
    temp_start, temp_end, temp_gene, temp_exon, temp_strand  = 0, 1000000000, " ", " ", " "
    for start, end, gene, exon, strand in data:
        if pos - int(temp_end) >= 0 and pos - int(start) <= 0:
            if temp_strand == "-" and strand == "-":
                target.add(temp_gene + "," + temp_exon + "," + temp_strand)
            else:
                target.add(gene + "," + str(exon) + "," + strand)
        if pos >= int(start) and pos <= int(end):
            return [gene + "," + str(exon) + "," + strand]
        temp_start, temp_end, temp_gene, temp_exon, temp_strand = start, end, gene, str(exon), strand
    return list(target)


def find_strand(pos, data):
    temp_start, temp_end, temp_gene, temp_exon, temp_strand = 0, 1000000000, " ", " ", " "
    for start, end, gene, exon, strand in data:
        if pos - int(temp_end) >= 0 and pos - int(start) <= 0:
            return strand
        if pos >= start and pos <= end:
            return strand
        temp_start, temp_end, temp_gene, temp_exon, temp_strand = start, end, gene, str(exon), strand


def anno_fusion(file, dic_data):
    title = ["chrom","pos", "orientation","chrom2","pos2","orientation2","type","size"
                                                    ,"score", "num_reads","num_reads_lib","tumor_1a.sorted.bam", "_"]
    fusion = pd.read_table(file, skiprows=5, names=title)
    fusion = fusion.loc[fusion["type"].isin(["ITX", "CTX"])]
    target_chr = ["chr" + str(i) for i in range(1, 22)] + ["chrX", "chrY"]
    fusion = fusion.loc[(fusion["chrom"].isin(target_chr)) & (fusion["chrom2"].isin(target_chr))]
    strand1, strand2, freq = [], [], []
    for chr1, pos1, chr2, pos2, orientation, orientation2, supports\
            in zip(fusion["chrom"], fusion["pos"], fusion["chrom2"], fusion["pos2"],
                   fusion["orientation"], fusion["orientation2"], fusion["num_reads"]):
        data1 = dic_data[chr1]
        if chr1 == chr2: data2 = data1
        else:
            data2 = dic_data[chr2]
        re = find_strand(int(pos1), data1)
        strand1.append(re)
        re = find_strand(int(pos2), data2)
        strand2.append(re)
        depth1 = int(orientation.strip("-").split("+")[0]) + int(orientation.strip("-").split("+")[1])
        depth2 = int(orientation2.strip("-").split("+")[0]) + int(orientation2.strip("-").split("+")[1])
        freq.append(supports * 2 / (depth1 + depth2))
    fusion["strand"], fusion["strand2"], fusion["freq"] = strand1, strand2, freq
    fusion = fusion[["chrom","pos", "strand", "chrom2", "pos2", "strand2", "num_reads", "orientation", "orientation2", "freq"]]
    fusion.columns = ["chrom", "pos", "strand", "chrom2", "pos2", "strand2", "supports", "depth1", "depth2", "freq"]
    fusion.to_csv(file + "_anno.txt", sep="\t", index=False)


def parser_gff(gff_file):
    gff_data = pd.read_table(gff_file, skiprows=6, names=["chr", "source", "type", "start", "end", "."
                                                          ,"strand", ".", "info"])
    chr = []
    for i in gff_data["chr"]:
        try:
            int(i)
            chr.append(True)
        except:
            chr.append(False)
    data = gff_data.loc[(gff_data["type"] == "exon") & chr]
    data = data[["chr", "start", "end","strand", "info"]]
    data["chr"] = data["chr"].astype(int)
    data["start"] = data["start"].astype(int)
    data["end"] = data["end"].astype(int)
    groups = data.groupby("chr")
    dic_data = defaultdict(list)
    for group, data in groups:
        data = data.sort_values(by="start")
        gene_exon = defaultdict(int)
        for start, end, strand, infomation in zip(data["start"], data["end"], data["strand"], data["info"]):
            p_gene = re.compile(r'gene_name=(.*?);')
            gene = p_gene.findall(infomation)[0]
            gene_exon[gene] += 1
            p_type = re.compile(r'gene_type=(.*?);')
            #if p_type.findall(infomation)[0] == "protein_coding":
            dic_data["chr" + str(group)].append([start, end, gene, gene_exon[gene], strand])
    return dic_data


def fuanno(args=None):
    dic_data = parser_gff(args.gff)
    pool = multiprocessing.Pool(processes=args.progress)
    for file in args.fusion_files:
        pool.apply_async(anno_fusion, (file, dic_data,))
    pool.close(); pool.join()


if __name__ == "__main__":
    files = glob("*.fusion")
    dic_data = parser_gff("hg37.gff")
    for file in files:
        re = anno_fusion(file, dic_data)
        re.to_csv(file + "_anno.txt", sep="\t", index=False)


