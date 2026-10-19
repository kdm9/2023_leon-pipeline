import argparse
import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

def get_args():
    """
    parses arguments from the command line
    :return: list of command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-s',
                        '--sqanti',
                        type=str,
                        help='sqanti classification',
                        required=True)
    parser.add_argument('-r',
                        '--reads',
                        type=str,
                        help='read support',
                        required=True)
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='Output file',
                        required=False)
    parser.add_argument('-a',
                        '--acc',
                        type=str,
                        help='accession',
                        required=False)

    args = parser.parse_args()
    return args


def main():

    args = get_args()
    # -------------1-------------
    # prepare sqanti output
    sqanti = pd.read_csv(args.sqanti, sep="\t")    
    sub = sqanti[["isoform","exons","length","structural_category","associated_gene","associated_transcript","subcategory","ORF_seq"]]
    sub["gene"] = sub["isoform"].str.split(".").str[0]
    # -------------3-------------
    # prepare read support
    reads=pd.read_csv(args.reads, sep="\t")
    r_sub = reads[["merge_trans_id","trans_read_count","gene_read_count"]]
    r_sub.columns = ["isoform", "isoform_reads", "gene_reads"]
    # -------------4-------------
    # add accession ID and merge all information
    merged_all = pd.merge(left=sub, right=r_sub, how='left', left_on='isoform', right_on='isoform')
    merged_all['accession'] = args.acc
    merged_all.to_pickle(args.output)

if __name__ == '__main__':
    main()
