#!/home/zhout/app/anaconda3/bin/python3
'''
main functional
'''

import time
from argparse import ArgumentParser
import os

__version__ = '0.1.3'


def main():
    parser = pre_parser()
    args = parser.parse_args()
    subCommand = args.subcommand
    if subCommand == 'cnv':
        from genecast_package.cnv_analysis import cnv
        cnv(args=args)
    elif subCommand == "snv":
        from genecast_package.snv_analysis import snv
        snv(args=args)
    elif subCommand == "ln":
        from genecast_package.make_ln import ln
        ln(args=args)
    elif subCommand == "fuanno":
        from genecast_package.fuan import fuanno
        fuanno(args=args)
    elif subCommand == "prediction":
        pass

def pre_parser():
    '''
    准备命令行解析
    :return:
    '''

    des = 'genecast analysis pipeline'
    epilog = '每个子命令的说明如下 genecast subcommond -h'

    parser = ArgumentParser(description=des, epilog=epilog)
    parser.add_argument("--version", action='version',
                        version="genecast: " + __version__)

    subparser = parser.add_subparsers(dest='subcommand')
    add_anno_fusion(subparser)
    add_make_ln(subparser)
    add_cnv(subparser)
    add_snv(subparser)
    add_circos(subparser)
    add_prediction(subparser)

    return parser


def add_output(parser):
    '''output option'''
    parser.add_argument('-o', '--outdir', required=False, default=os.getcwd(),
                        help='the dir of output, default is a_vs_b')


def add_make_ln(subparser):
    parser = subparser.add_parser('ln', help='make ln module, the template file please contact the author')
    parser.add_argument('-sf', '--sample_file', required=True, help='the dir of group1', type=str)
    parser.add_argument('-r', '--research', choices=("yes", "no"), required=False, default="yes",
                        help='if you sample is scientific research cooperation please choose yes,else choose no')


def add_anno_fusion(subparser):
    parser = subparser.add_parser('fuanno', help='anno fusion resulut of breakdancer module, the fusion resulut please contact the author')
    parser.add_argument('fusion_files', nargs='*', help="fusion resulut of breakdancer files (.fusion)")
    parser.add_argument('-gff', required=True, help='the gff file', type=str)
    parser.add_argument('-p', '--progress', required=False, default=1, type=int,
                        help='parallel anno fusion result')


def add_circos(subparser):
    parser = subparser.add_parser('circos',
                                  help='arranged the data for plot circos')
    # prediction参数
    pass


def add_prediction(subparser):
    parser = subparser.add_parser('prediction',
                                  help='accept a dataframe that row is sample, columns is genes, a group file, '
                                       'feature gene selection and prediction new sample')
    # prediction参数
    parser.add_argument("-td1", "--train_data", required=True, help='the data file name', type=str)
    parser.add_argument("-td2", "--test_data", required=True, help='the data file name', type=str)
    parser.add_argument("-g", '--group', default='group.txt', type=str,
                        help='the group1 of the group file must be group_a, the group2 of the group file must be group_b'
                             'such as A  B\n'
                             '        s1 s2\n'
                             '        s3 s4\n'
                             '        .. ..\n'
                             'table split')
    add_common_parameter(parser)


def add_snv(subparser):
    parser = subparser.add_parser('snv', help='snv analysis module')
    parser.add_argument("-a", '--group1', nargs='*', help="snv resulut of annovar files (*snp*.hg19_multianno.txt)")
    parser.add_argument("-b", '--group2', nargs='*', help="snv resulut of annovar files (*snp*.hg19_multianno.txt)")
    # snv参数
    parser.add_argument("-r", '--ratio', default=0, type=float,
                        help='filter somatic ratio  default=0')
    parser.add_argument("-os", '--one_strand', default=1, type=int,
                        help='strand preference filter if two strand both greater than default,'
                             'wo believe in it is a somatic default=1')
    parser.add_argument("-ts", '--two_strand', default=5, type=int,
                        help='if one of the two strand greater than defailt'
                             ',whatever wo believe in it is a somatic default=10')
    parser.add_argument("-cal", '--cal_type', default='num', type=str,
                        choices=('num', 'mean'),
                        help='which type to cal  default=num')
    parser.add_argument("-dt", '--data_type', default='snv', type=str,
                        choices=('snv', 'snp', "indel"),
                        help='snv contain snp and indel  default=snv')
    parser.add_argument("-locus", default=False, type=bool,
                        choices=(False, True),
                        help='if True, all calculate will base on locus level')
    parser.add_argument("-somatic", default='Y', type=str,
                        choices=('y', 'n'),
                        help='if you only have Blood cell and cfdna you should choose n else you have Blood cell and cfdna you should choose y; default=y')
    parser.add_argument('-circos', default=False, type=bool,
                        choices=(True, False), help="if you need plot circos, please choose True, depault=False")
    add_common_parameter(parser)


def add_cnv(subparser):
    parser = subparser.add_parser('cnv', help='cnv analysis module')
    parser.add_argument("-a", '--group1', nargs='*', help="cnv resulut of cnvkit files (*.cnr)")
    parser.add_argument("-b", '--group2', nargs='*', help="cnv resulut of cnvkit files (*.cnr)")
    # cnv参数
    parser.add_argument("-dt", '--data_type', default='log2', type=str,
                        choices=('log2', 'cn'),
                        help='if you choose cn, please must use cnvkit call the cnr file and the call result filename '
                             'must be *call, if you choose, the target filename must be *cnr, default is log2')
    add_common_parameter(parser)


def add_common_parameter(parser):
    #input

    parser.add_argument('-hg', '--host_gene', required=True, type=str,
                        help='host gene file, please make sure it only contain one colomn and its title must be gene,'
                             'or panal bed file also be ok')
    parser.add_argument('-z', '--z_score', required=False, type=int,
                        choices=(0, 1, 3), default=0,
                        help='z_score for heatmap, z = (x - mean)/std if 0, z for row if 1 z for column if 3 do nothing'
                             'default=0')
    parser.add_argument("-cmp", default="bwr", type=str,
                        help='colormaps for heatmap, default=bwr, refer http://matplotlib.org/examples/color/colormaps_reference.html')
    parser.add_argument('-fsm', "--feature_selection_method", required=False,
                        choices=('wilcox', "fisher", 'pearsonr', "Lasso", "logistic", "RandomizedLasso", "RandomForest", "Wrapper", "variance"),
                        type=str, default="logistic",
                        help='please choose a feature selection method, default is logistic')
    parser.add_argument('-threshold', required=False, type=float, default=0,
                        help='the threshold weight coefficient of choose feature gene, default=0')
    parser.add_argument('-criterion', required=False, type=str, default="aic",
                        choices=('bic', 'aic'),
                        help=''''bic' | 'aic' The type of criterion to use, default=aic''')
    parser.add_argument('-pm', "--prediction_method", required=False,
                        choices=('LinearSVC', 'SVC', "logistic"),
                        type=str, default="SVC", help='please choose a prediction method, default is SVC')
    parser.add_argument('-p', '--pval', required=False, default=0.05, type=float,
                        help='if the method of feature selection is statistical test, please provide this parameter,'
                             'default is 0.05')
    parser.add_argument('-penalty', required=False, type=str, default="l2",
                        help='''str, 'l1' or 'l2', default: 'l2' Used to specify the norm used in the penalization. The 'newton-cg',
                                 'sag' and 'lbfgs' solvers support only l2 penalties''')
    parser.add_argument('-C', required=False, type=float, default=1,
                        help='smaller values specify stronger regularization, default=1')
    parser.add_argument('-n', '--n_folds', required=False, type=int, default=3,
                        help='int, default=3 Number of folds. Must be at least 2')
    parser.add_argument('-alpha', required=False, type=float, default=0.25,
                        help='float, default=0.25 parameter for lasso.')
    parser.add_argument('-ne', '--n_estimators', required=False, type=int, default=20,
                        help='int, default=20 number of tree parameter for RandomForestRegressor.')
    parser.add_argument('-nf', '--n_feature', required=False, type=int, default=15,
                        help='int, default=15 number of feature')
    parser.add_argument('-step', required=False, type=int, default=200,
                        help='int, default=200 parameter for Wrapper_LogisticRegression')
    parser.add_argument('-cluster_method', required=False, type=str, default="complete",
                        choices=('complete', 'average', "weighted", "single"),
                        help='str, default=complete parameter for cluster_method')
    parser.add_argument('-save', required=False, type=str, default="png",
                        choices=('png', 'pdf'),
                        help='str, default=png parameter for save file format')
    add_output(parser)

if __name__ == '__main__':
    main()

