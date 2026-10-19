import argparse
import subprocess
from os import listdir
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import os.path
import glob
from rg import main

def get_uniques(sqanti, ref, acc, tool, version, name):
    # df reduction
    sqanti_genes = sqanti[['isoform', 'associated_transcript', 'structural_category']]
    sqanti_old = sqanti_genes.loc[sqanti_genes['associated_transcript'] != 'novel']
    sqanti_ref = sqanti_old.loc[sqanti_old['associated_transcript'].isin(ref['genes'])]
    sqanti_ref_full = sqanti_ref.loc[sqanti_ref['structural_category'] == 'full-splice_match']
    # merged = sqanti_ref_full['isoform'].str.split('.').str[0] + '_' + sqanti_ref_full['associated_transcript'] + '.' + acc + '.' + tool + '.' + version
    merged = sqanti_ref_full['isoform'].str.split('.').str[0] + '_' + sqanti_ref_full[
        'associated_transcript']

    merged_values = merged.value_counts().to_frame()
    merged_values.columns = ['values']
    merged_values['accession'] = acc
    merged_values['tool'] = tool
    merged_values['version'] = version
    merged_values['reference'] = name
    return merged_values

# TODO: function to generate all violin plots in one

# TODO: function to create all pairwise plots


def main():
    ref_names = ['Busco', 'NLR', 'Housekeeping']
    references = ['/home/leon/Documents/master_thesis/Analysis/sqanti/ref_genes/busco_Araport11_20221221.txt',
          '/home/leon/Documents/master_thesis/Analysis/sqanti/ref_genes/nlrs_in_orthogroups_20221221.txt',
          '/home/leon/Documents/master_thesis/Analysis/sqanti/ref_genes/arat_housekeeping.txt']

    # files = ["/home/leon/Documents/master_thesis/Analysis/sqanti/at9852-tama-vs-liftoff_v2_classification.txt", "/home/leon/Documents/master_thesis/Analysis/sqanti/at9852-tama-vs-liftoff_v3_classification.txt"]
    files = glob.glob("/home/leon/Documents/master_thesis/Analysis/sqanti/v*/*")

    df_all = pd.DataFrame()
    l = 0
    # for everey reference, do all counting and append to df
    for r_path, r_name in zip(references, ref_names):
        ref = pd.read_csv(r_path, header=None, names=['genes'])
        for path in files:
            df = pd.read_csv(path, sep="\t")
            basename = os.path.basename(path)
            # v = basename.split('_')[1]
            v = path.split('/')[6]
            acc = basename.split('-')[0]
            tool = basename.split('-')[1]
            values = get_uniques(df, ref, acc, tool, v, r_name)
            l += len(values) #check if concat worked as inteded
            df_all = pd.concat([df_all, values])

    df_all['log_val'] = np.log(df_all['values'])
    df_all['bin'] = df_all['values']
    df_all['bin'][df_all['bin'] >= 50] = 50
    print(df_all)
    print(l)
    df_all.to_csv('/home/leon/Documents/master_thesis/Analysis/sqanti/all.tsv', sep='\t', header=True)

    # ref_specific dfs
    busco = df_all[df_all['reference'] == 'Busco']
    nlr = df_all[df_all['reference'] == 'NLR']
    house = df_all[df_all['reference'] == 'Housekeeping']

    # ref_comparisons
    busco_nlr = df_all[df_all['reference'] == 'Busco' and df_all['reference'] == 'NLR']
    busco_nlr_v2 = busco_nlr[busco_nlr['version'] == 'v2']
    busco_nlr_v3 = busco_nlr[busco_nlr['version'] == 'v3']
    print(busco_nlr_v2)
    busco_house = df_all[df_all['reference'] == 'Busco' and df_all['reference'] == 'Housekeeping']
    nlr_house = df_all[df_all['reference'] == 'NLR' and df_all['reference'] == 'Housekeeping']


    fig = plt.figure()
    sns.violinplot(data=df_all, y='bin', x='accession', hue='version', split=True)
    plt.show()
    # fig.savefig('/home/leon/Documents/master_thesis/Analysis/sqanti/normal.png', dpi=fig.dpi)


if __name__ == '__main__':
    main()