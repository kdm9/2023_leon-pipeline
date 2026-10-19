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
    parser.add_argument('-m',
                        '--merged',
                        type=str,
                        help='sqanti classification',
                        required=True)
    parser.add_argument('-o',
                        '--output',
                        type=str,
                        help='Output file',
                        required=False)

    args = parser.parse_args()
    return args


def main():

    args = get_args()

    sqanti = pd.read_csv(args.sqanti, sep="\t")
    # merged_all = pd.read_pickle("./sqanti/{}_tama_liftoff_merged_{}.pkl".format(acc, run))
    merged_all = pd.read_pickle(args.merged)

    sub = sqanti[["isoform","exons","length","structural_category","associated_gene","associated_transcript","subcategory","ORF_seq"]]
    sub["gene"] = sub["isoform"].str.split(".").str[0]

    # -------------5-------------
    # Generate list of subdataframes (Sqanti ORF analysis)
    genes = list(sub.groupby("gene"))

    gene2ORF = {}
    for i in genes:
        gene2ORF[i[0]]={}
        orfs = i[1]["ORF_seq"].value_counts().index.tolist()
        for orf in orfs:
            gene2ORF[i[0]][orf] = i[1].loc[i[1]["ORF_seq"] == orf]["isoform"].tolist()


    ORF2iso = pd.DataFrame()

    for gene in gene2ORF.keys():
        # get the isoforms per ORF + calculate the number of isoforms per ORF
        isoforms=[]
        n_isoforms=[]
        for orf in gene2ORF[gene].keys():
            isoforms.append(gene2ORF[gene][orf])
            n_isoforms.append(len(gene2ORF[gene][orf]))
        # get read support and sqanti subcategory per isoform
        read_support=[]
        subcategories=[]
        for l in isoforms:
            read_supp_list=[]
            category_list = []
            for iso in l:
                # read support
                rs=merged_all.loc[merged_all["isoform"] == iso]["isoform_reads"]
                read_supp_list.append(int(rs.to_string().split("    ")[1]))
                # subcategory
                subcat = merged_all.loc[merged_all["isoform"] == iso]["subcategory"]
                category_list.append(subcat.to_string().split("    ")[1])
            read_support.append(read_supp_list)
            subcategories.append(category_list)
        gene={"Gene":gene,
              "n ORFs":len(gene2ORF[gene]),
              "n Isoforms": sum(n_isoforms),
              "Sqanti ORFs":list(gene2ORF[gene].keys()),
              "Isoforms":isoforms,
              "n_Isoforms_ORF":[len(x) for x in isoforms],
              "Subcategories":subcategories,
              "Read_Support":read_support,
              "Read_Support_ORF":[sum(x) for x in read_support]}
        ORF2iso=pd.concat([ORF2iso, pd.DataFrame([gene])])


    # pb_liftoff_sqanti.to_pickle("./sqanti/{}_pb_liftoff_sqanti_orf_{}.pkl".format(acc, run)
    ORF2iso.to_pickle(args.output)

if __name__ == '__main__':
    main()
