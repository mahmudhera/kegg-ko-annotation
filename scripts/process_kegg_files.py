#!/usr/bin/env python
import argparse
import pandas as pd

def read_organism_table(org_table_filename):
    # read dataframe
    df = pd.read_csv(org_table_filename, delimiter='\t')

    # keep only bacteria
    df['lineage'] = df['lineage'].str.lower()
    df = df[ df['lineage'].str.contains('bacteria') ]

    # take interesting things
    gb_ids = df['gb_ncbi_seq_id'].tolist()
    org_codes = df['org_code'].tolist()
    names = df['name'].tolist()

    # keep only non-empty ncbi things here
    ret_dic = {}
    for gb_id, org_code, name in list( zip(gb_ids, org_codes, names) ):
        if gb_id == '[]':
            continue
        else:
            gb_id = gb_id[1:-1].replace('\'', '').split(',')
            ret_dic[org_code] = gb_id

    return ret_dic

def read_organism_table_for_single_chr_organisms(org_table_filename):
    # read dataframe
    df = pd.read_csv(org_table_filename, delimiter='\t')

    # keep only bacteria
    df['lineage'] = df['lineage'].str.lower()
    df = df[ df['lineage'].str.contains('bacteria') ]

    # take interesting things
    gb_ids = df['gb_ncbi_seq_id'].tolist()
    org_codes = df['org_code'].tolist()
    names = df['name'].tolist()

    # keep only non-empty ncbi things here
    ret_dic = {}
    for gb_id, org_code, name in list( zip(gb_ids, org_codes, names) ):
        if gb_id == '[]':
            continue
        else:
            gb_id = gb_id[1:-1].replace('\'', '').split(',')
            if len(gb_id) > 1:
                continue
            ret_dic[org_code] = gb_id

    return ret_dic

def read_gene_ids(kegg_gene_filename):
    f = open(kegg_gene_filename, 'r')
    lines = f.readlines()[1:]
    f.close()
    return [ line.split(' ')[0].split('\t')[0] for line in lines ]


def parse_args():
    parser = argparse.ArgumentParser(description="This preprocesses kegg_db data files. It opens the organisms table, locates the genomes that has GenBank or RefSeq entries available.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--org_table_file', type=str, help="The organisms table file.")
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = parse_args()
    org_table_filename = args.org_table_file
    org_code_to_ncbi_ids = read_organism_table(org_table_filename)
    print( len( org_code_to_ncbi_ids.keys() ) )

    org_code_to_ncbi_ids_single_chr = read_organism_table(org_table_filename)
    print( len( org_code_to_ncbi_ids_single_chr.keys() ) )

    sample_gene_filename = 'data/sce_kegg_genes.txt'
    gene_ids = read_gene_ids(sample_gene_filename)
    #print(gene_ids)

    ## We now know a lot
    ## The same exact coords to be used
    ## But, how to locate the contig??
