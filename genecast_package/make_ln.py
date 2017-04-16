## this tool is for making ln of bam or mpileup of GeneticTest or research folder
## author: taozhou
## email: zhou.tao@genecast.com.cn

import pandas as pd
import os
import sh
from glob import glob


def make_ln_GeneticTest(groups, search_dir, target_dir, which="bam"):
    if which == "bam":
        source_name = "/mapped/%s_1a.sorted.bam"
        target_name = "sorted.bam"
    else:
        source_name = "/mpileup/%s_1a.mpileup.gz"
        target_name = "mpileup.gz"
    for group, group_info in groups:
        try:
            os.mkdir(group + "_%s" % which)
        except FileExistsError:
            pass
        for id in group_info["patient"]:
            for file in search_dir:
                if id == file.split("/")[-1]:
                    target_file = file
            sh.ln("-s", target_file + source_name % "blood", target_dir + "/" + group + "_%s" % which + id + "_bc" + target_name)
            sh.ln("-s", target_file + source_name % "plasma", target_dir + "/" + group + "_%s" % which + id + "_cf" + target_name)


def make_ln_research(groups, search_dir, target_dir, which="bam"):
    if which == "bam":
        source_name = ".sorted.filtered.bam"
        source_folder = "/mapped/"
        target_name = ".bam"
    else:
        source_name = ".mpileup.gz"
        source_folder = "/mpileup/"
        target_name = ".mpileup.gz"
    for group, group_info in groups:
        try:
            os.mkdir(group + "_%s" % which)
        except FileExistsError:
            pass
        i = 0
        for s, id, line, date, type in zip(group_info["patient"], group_info["sampleid"], group_info["line"],
                                           group_info["data"], group_info["type"]):
            for file in search_dir:
                if date.lower() in file or date in file:
                    target_file = file
                    sample_name = "S" + str(line) + "_" + id
            sh.ln("-s", target_file + source_folder + sample_name + source_name,
                  target_dir + "/" + group + "_%s" % which + "/" + type + "_" + s + target_name)


def get_file_ln_info(file, research):
    ## file: use_sample_train.txt
    target_dir = os.getcwd()
    sample_info = pd.read_table(file)
    groups = sample_info.groupby(sample_info["group"])
    if research != "yes":
        search_dir = glob("/work-z/shared/GeneticTest/*")
        make_ln_GeneticTest(groups, search_dir, target_dir, which="bam")
        make_ln_GeneticTest(groups, search_dir, target_dir, which="mpileup")
    else:
        search_dir = glob("/work-z/shared/Mapping/*")
        make_ln_research(groups, search_dir, target_dir, which="bam")
        make_ln_research(groups, search_dir, target_dir, which="mpileup")
    return sample_info["group"]


def ln(sf=None, research="yes"):
    get_file_ln_info(sf, research)


if __name__ == "__main__":
    get_file_ln_info("use_sample_train1.txt", research)