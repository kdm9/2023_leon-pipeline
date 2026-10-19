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
                        help='Output file',
                        required=False)

    args = parser.parse_args()
    return args


def main():

    args = get_args()
    
    # -------------1-------------
    # prepare sqanti output
    # sqanti = pd.read_csv("/ebio/abt6_projects7/dliso/dlis/leon_pipeline/output/v2/collapse_pb/{}_liftoff_sqanti/{}-pb-vs-liftoff_classification.txt".format(acc,acc), sep="\t")
    sqanti = pd.read_csv(args.sqanti, sep="\t")
    
    
    sub = sqanti[["isoform","exons","length","structural_category","associated_gene","associated_transcript","subcategory","ORF_seq"]]
    sub["gene"] = sub["isoform"].str.split(".").str[0]+"."+sub["isoform"].str.split(".").str[1]

    # -------------3-------------
    # prepare read support
    # reads=pd.read_csv("/ebio/abt6_projects7/dliso/dlis/leon_pipeline/output/v2/collapse_pb/at9900.collapsed.pb.v2.abundance.txt",sep="\t", skiprows=3)
    reads=pd.read_csv(args.reads, sep="\t", skiprows=3)
    r_sub = reads[["pbid","fl_assoc"]]
    r_sub.columns = ["isoform", "isoform_reads"]
    r_sub["gene"]=r_sub["isoform"].str.split(".").str[0]+r_sub["isoform"].str.split(".").str[1]
    for i in r_sub.groupby("gene"):
        r_sub.loc[r_sub["gene"]==i[0], "gene_reads"] = i[1]["isoform_reads"].sum()
    r_sub["gene_reads"] = r_sub["gene_reads"].astype("int")

    # -------------4-------------
    # add accession ID and merge all information
    merged_all = pd.merge(left=sub, right=r_sub[["isoform","isoform_reads","gene_reads"]], how='left', left_on='isoform', right_on='isoform')
    merged_all['accession'] = args.acc

    # merged_all.to_pickle("./sqanti/{}_pb_liftoff_merged_{}.pkl".format(acc, run))
    merged_all.to_pickle(args.output)

if __name__ == '__main__':
    main()
