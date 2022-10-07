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

def read_gene_id_with_ko(kegg_gene_filename):
    df = pd.read_csv(kegg_gene_filename, delimiter='\t')
    df = df[ ~pd.isnull(df['koid']) ]
    return(df['kegg_gene_id'].tolist())

def parse_args(): # pragma: no cover
    parser = argparse.ArgumentParser(description="This preprocesses kegg_db data files. It opens the organisms table, locates the genomes that has GenBank or RefSeq entries available.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--org_table_file', type=str, help="The organisms table file.")
    args = parser.parse_args()
    return args

def main(): # pragma: no cover
    args = parse_args()
    org_table_filename = args.org_table_file

    # get list of all bacteria organisms
    org_code_to_ncbi_ids = read_organism_table(org_table_filename)
    print('Num of all bacteria organism:')
    print( len( org_code_to_ncbi_ids.keys() ) )

    # get all bacteria with single chromosome
    print('Num of all bacteria organism with single chromosome:')
    org_code_to_ncbi_ids_single_chr = read_organism_table(org_table_filename)
    print( len( org_code_to_ncbi_ids_single_chr.keys() ) )

    # for all these organisms
        # if kegg gene file not existing, then skip

    # report reduced num of organisms

    # for all these organisms:
        # look in the directory and find all genes with KOs present

    # report num of all genes here

    # for all these organisms
        # get a list of all the genes with KOs present
        # make a directory for this organism
        # write fasta file in that

        # for all the genes:
            # look up KEGG website
            # get the position of that gene
            # create an entry in the mapping file

        # write mapping file

    # print number of genes

    sample_gene_filename = 'data/sce_kegg_genes.txt'
    gene_ids = read_gene_ids(sample_gene_filename)
    print(len(gene_ids))

    sample_gene_filename = 'data/sce_kegg_genes.txt'
    gene_ids = read_gene_id_with_ko(sample_gene_filename)
    print(len(gene_ids))

if __name__ == '__main__': # pragma: no cover
    main()

    ## We now know a lot
    ## The same exact coords to be used
    ## But, how to locate the contig??
