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
    parser.add_argument('-f',
                        '--fasta',
                        type=str,
                        help='orf fastas',
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

def read(file_name):
    try:
        fasta = open(file_name, 'r')
        fasta_pairs = {}
        header = ''
        # check every line for a new header
        # if so, set as key and set non header lines as value
        for line in fasta:
            if line[0] == '>':
                header = line.strip('\n')
                fasta_pairs[header] = ''
            else:
                fasta_pairs[header] += line.strip('\n')
        # return dict of fasta entries
        return fasta_pairs
    except IOError:
        print("File \"{}\" could not be opened.".format(file_name))

def main():

    args = get_args()

    # merged_all = pd.read_pickle("./sqanti/{}_tama_liftoff_merged_{}.pkl".format(acc, run))
    merged_all = pd.read_pickle(args.merged)

    # -------------5-------------
    # Generate list of subdataframes (Transdecoder ORF analysis)
    # fasta=read("/ebio/abt6_projects7/dliso/dlis/leon_pipeline/output/v2/collapse_pb/{}_transdecoder/{}_transcripts.fasta.transdecoder.genome.sorted.protein.fasta".format(acc,acc))
    fasta=read(args.fasta)

    trans_df=pd.DataFrame()
    for key in fasta.keys():
        gene = key.split(".")[0].replace(">","")
        iso = key.split(" ")[0].replace(">","").rsplit('.', 1)[0]
        p = key.split(" ")[0].replace(">","").rsplit('.', 1)[1]    
        row={"Gene":gene,
            "isoform":iso,
            "ORF_variant":p,
            "Transdecoder ORF":fasta[key].replace("*","")}
        trans_df=pd.concat([trans_df,pd.DataFrame([row])])

    trans_genes = list(trans_df.groupby("Gene"))
    gene2ORF_trans = {}
    for i in trans_genes:
        gene2ORF_trans[i[0]]={}
        orfs = i[1]["Transdecoder ORF"].value_counts().index.tolist()
        for orf in orfs:
            gene2ORF_trans[i[0]][orf] = i[1].loc[i[1]["Transdecoder ORF"] == orf]["isoform"].tolist()

    ORF_trans=pd.DataFrame()
    for gene in gene2ORF_trans.keys():
        # get the isoforms per ORF + calculate the number of isoforms per ORF
        isoforms=[]
        n_isoforms=[]
        for orf in gene2ORF_trans[gene].keys():
            isoforms.append(gene2ORF_trans[gene][orf])
            n_isoforms.append(len(gene2ORF_trans[gene][orf]))
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
              "n ORFs":len(gene2ORF_trans[gene]),
              "n Isoforms": sum(n_isoforms),
              "Transdecoder ORFs":list(gene2ORF_trans[gene].keys()),
              "Isoforms":isoforms,
              "n_Isoforms_ORF":[len(x) for x in isoforms],
              "Subcategories":subcategories,
              "Read_Support":read_support,
              "Read_Support_ORF":[sum(x) for x in read_support]}
        ORF_trans=pd.concat([ORF_trans, pd.DataFrame([gene])])

    # ORF_trans.to_pickle("./orf/{}_pb_liftoff_transdecoder_orf_{}.pkl".format(acc, run))
    ORF_trans.to_pickle(args.output)

if __name__ == '__main__':
    main()
