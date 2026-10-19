import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

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

accessions= ['at6923','at6929','at7143','at8285','at9104','at9336','at9503','at9578','at9744','at9762','at9806','at9830','at9847','at9852','at9879','at9883','at9900','col0']

#runs=["full","iso2","gene2"]
runs=["full"]


for run in tqdm(runs):
    pb_liftoff_sqanti = pd.DataFrame()
    pb_liftoff_trans = pd.DataFrame()
    for acc in tqdm(accessions):
        print(acc)
        # -------------0-------------
        # make directories
        os.makedirs("{}/{}/isoform".format(run,acc), exist_ok=True)
        os.makedirs("{}/{}/gene".format(run,acc), exist_ok=True)
        # -------------1-------------
        # prepare sqanti output
        sqanti = pd.read_csv("/ebio/abt6_projects7/dliso/dlis/leon_pipeline/output/v2/collapse_pb/{}_liftoff_sqanti/{}-pb-vs-liftoff_classification.txt".format(acc,acc), sep="\t")
        sub = sqanti[["isoform","exons","length","structural_category","associated_gene","associated_transcript","subcategory","ORF_seq"]]
        sub["gene"] = sub["isoform"].str.split(".").str[0]+sub["isoform"].str.split(".").str[1]

        # -------------3-------------
        # prepare read support
        reads=pd.read_csv("/ebio/abt6_projects7/dliso/dlis/leon_pipeline/output/v2/collapse_pb/{acc}.collapsed.pb.v2.abundance.txt".format(acc),sep="\t", skiprows=3)
        r_sub = reads[["pbid","fl_assoc"]]
        r_sub.columns = ["isoform", "isoform_reads"]
        r_sub["gene"]=r_sub["isoform"].str.split(".").str[0]+r_sub["isoform"].str.split(".").str[1]
        for i in r_sub.groupby("gene"):
            r_sub.loc[r_sub["gene"]==i[0], "gene_reads"] = i[1]["isoform_reads"].sum()
        r_sub["gene_reads"] = r_sub["gene_reads"].astype("int")

        # -------------4-------------
        # add accession ID and merge all information
        merged_all = pd.merge(left=sub, right=r_sub[["isoform","isoform_reads","gene_reads"]], how='left', left_on='isoform', right_on='isoform')
        merged_all['accession'] = acc

        # set dataset to work with
        if run == "iso2":
            merged_all=merged_all[merged_all["isoform_reads"]>1]
        elif run == "gene2":
            merged_all=merged_all[merged_all["gene_reads"]>1]

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

        if run == "full":
            merged_all.to_pickle("./sqanti/{}_pb_liftoff_merged_{}.pkl".format(acc, run))
            pb_liftoff_sqanti.to_pickle("./sqanti/{}_pb_liftoff_sqanti_orf_{}.pkl".format(acc, run))
            pn_liftoff_sqanti = pd.concat([pb_liftoff_sqanti, ORF2iso])

    if run == "full":
        pn_liftoff_sqanti.to_pickle("./pb_liftoff_sqanti_orf_{}.pkl".format(run))
        pb_liftoff_sqanti.to_csv("./pb_liftoff_sqanti_orf_{}.tsv".format(run), sep="\t")

    print("\n##############################Tama-Sqanti-ORF Done!##############################")

    # ORF transdecoder
    for acc in tqdm(accessions):
        merged_all = pd.read_pickle("./sqanti/{}_tama_liftoff_merged_{}.pkl".format(acc, run))
        # -------------5-------------
        # Generate list of subdataframes (Transdecoder ORF analysis)
        fasta=read("/ebio/abt6_projects7/dliso/dlis/leon_pipeline/output/v2/collapse_pb/{}_transdecoder/{}_transcripts.fasta.transdecoder.genome.sorted.protein.fasta".format(acc,acc))
        
        trans_df=pd.DataFrame()
        for key in fasta.keys():
            gene = key.rsplit(".",1)[0].replace(">","").rsplit('.', 2)[0]
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
                  "Sqanti ORFs":list(gene2ORF_trans[gene].keys()),
                  "Isoforms":isoforms,
                  "n_Isoforms_ORF":[len(x) for x in isoforms],
                  "Subcategories":subcategories,
                  "Read_Support":read_support,
                  "Read_Support_ORF":[sum(x) for x in read_support]}
            ORF_trans=pd.concat([ORF_trans, pd.DataFrame([gene])])
        
        if run == "full":     
            ORF_trans.to_pickle("./orf/{}_pb_liftoff_transdecoder_orf_{}.pkl".format(acc, run))
            pb_liftoff_trans = pd.concat([pb_liftoff_trans, ORF_trans])

    if run == "full":
        pb_liftoff_trans.to_pickle("./pb_liftoff_transdecoder_orf_{}.pkl".format(run))
        pb_liftoff_trans.to_csv("./pb_liftoff_transdecoder_orf_{}.tsv".format(run), sep="\t")
    
    print("\n##############################Tama-Transdecoder-ORF Done!##############################")
