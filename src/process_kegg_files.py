#!/usr/bin/env python
import argparse
import pandas as pd
import requests
from bs4 import BeautifulSoup

def read_organism_table(org_table_filename):
    """
    Takes as argument the organism table filename with full path.
    Returns a dictionary indexed by 3-character organism code, mapping to a list
    of contig ids.
    :param str org_table_filename: full path to organism filename
    :return: dictionary indexed by 3-character organism code, mapping to a list
    of contig ids.
    :rtype: dic
    """
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
    """
    Takes as argument the organism table filename with full path.
    Returns a dictionary indexed by 3-character organism code, mapping to a list
    of contig ids. Only those bacteria with single chromosome are kept.
    :param str org_table_filename: full path to organism filename
    :return: dictionary indexed by 3-character organism code, mapping to a list
    of contig ids.
    :rtype: dic
    """
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
    """
    Given a kegg gene filename, returns list of all genes
    :param str kegg_gene_filename: filename with full path
    :return: all gene ids in that file
    :rtype: list
    """
    f = open(kegg_gene_filename, 'r')
    lines = f.readlines()[1:]
    f.close()
    return [ line.split(' ')[0].split('\t')[0] for line in lines ]

def read_gene_id_with_ko(kegg_gene_filename):
    """
    Given a kegg gene filename, returns list of all genes with KO ids available
    :param str kegg_gene_filename: filename with full path
    :return: all gene ids in that file
    :rtype: list
    """
    df = pd.read_csv(kegg_gene_filename, delimiter='\t')
    df = df[ ~pd.isnull(df['koid']) ]
    return(df['kegg_gene_id'].tolist())

def read_gene_id_with_kos_labeled(kegg_gene_filename):
    """
    Given a kegg gene filename, returns list of all genes with KO ids available
    :param str kegg_gene_filename: filename with full path
    :return: all gene ids in that file, along with their assigned KO ids
    :rtype: list of 2-tuples, tuple[0] is the gene id, tuple[1] is the KOid
    """
    df = pd.read_csv(kegg_gene_filename, delimiter='\t')
    df = df[ ~pd.isnull(df['koid']) ]
    return(  list( zip( df['kegg_gene_id'].tolist(), df['koid'].tolist() ) )  )

def download_page(url):
    while True:
        response = requests.get(url)
        response.raise_for_status()
        if response.status_code == 200:
            break
    return response.text

def read_gene_start_and_end_positions(kegg_gene_id):
    """
    Given kegg gene id, looks up KEGG website, reads response and parses it
    to obtain start and end position of that gene.
    :param str kegg_gene_id: id of the gene in kegg database
    :return: 2-tuple, tuple[0] is start, tuple[1] is end position
    :rtype: tuple of int
    """
    url = "https://www.genome.jp/entry/" + kegg_gene_id
    content = download_page(url)
    soup = BeautifulSoup(content, 'html.parser')

    for row in soup.table.find_all('tr'):
        row_header = row.th
        row_cell = row.td
        if row_header is None:
            continue
        if row_header.get_text() == "Position":
            full_text = row_cell.get_text()
            full_text = full_text.replace('complement(', '')
            full_text = full_text.replace(')', '')
            start_end_merged = full_text.split('\n')[0].split(':')[1]
            start = int( start_end_merged.split('..')[0] )
            end = int( start_end_merged.split('..')[1] )
            return start, end

    return (-1, -1)

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

    print(gene_ids[:10])

    sample_gene_filename = 'data/sce_kegg_genes.txt'
    gene_ids = read_gene_id_with_ko(sample_gene_filename)
    print(len(gene_ids))

    print(gene_ids[:10])

if __name__ == '__main__': # pragma: no cover
    main()

    ## We now know a lot
    ## The same exact coords to be used
    ## But, how to locate the contig??
