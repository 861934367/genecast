## use gff anno the fusion result of breakdancer
import pandas as pd
from glob import glob
import re
from collections import defaultdict
import multiprocessing


def find_gene(pos, data):
    for start, end, gene, exon in data:
        if pos >= int(start) and pos <= int(end):
            return (gene, exon)


def anno_fusion(file, dic_data):
    title = ["chr1","pos1", "orientation1","chr2","pos2","orientation2","type","size"
                                                    ,"score	num_reads","num_reads_lib","tumor_1a.sorted.bam", "_"]
    fusion = pd.read_table(file, skiprows=5, names=title)
    gene1, exon1, gene2, exon2 = [], [], [], []
    for chr1, pos1, chr2, pos2 in zip(fusion["chr1"], fusion["pos1"], fusion["chr2"], fusion["pos2"]):
        data1 = dic_data[chr1]
        if chr1 == chr2: data2 = data1
        else:
            data2 = dic_data[chr2]
        re = find_gene(pos1, data1)
        if re == None: gene1.append(None); exon1.append(None)
        else:gene1.append(re[0]); exon1.append("exon" + re[1])
        re = find_gene(pos2, data2)
        if re == None:
            gene2.append(None); exon2.append(None)
        else:
            gene2.append(re[0]); exon2.append("exon" + re[1])
    fusion["gene1"], fusion["exon1"], fusion["gene2"], fusion["exon2"] = gene1, exon1, gene2, exon2
    fusion.to_csv(file + "_anno.txt", sep="\t", index=False)


def parser_gff(gff_file):
    gff_data = pd.read_table(gff_file, skiprows=6, names=["chr", "source", "type", "start", "end", "."
                                                          ,"stand", ".", "info"])
    data = gff_data.loc[gff_data["type"] == "exon"]
    data = data[["chr", "start", "end", "info"]]
    dic_data = defaultdict(list)
    for chr, start, end, infomation in zip(data["chr"], data["start"], data["end"], data["info"]):
        p_gene = re.compile(r'gene_name=(.*?);')
        p_exon = re.compile(r'exon_number=(.*?);')
        dic_data["chr" + str(chr)].append((start, end, p_gene.findall(infomation)[0], p_exon.findall(infomation)[0]))
    return dic_data

#
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


