#!/usr/bin/env python
import argparse
import pandas as pd
import requests
from bs4 import BeautifulSoup
from os.path import exists
from Bio import SeqIO
import os
import subprocess
from tqdm import tqdm
import random

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

def read_organism_table_for_single_chr_bacteria(org_table_filename):
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
    :rtype: list of 4-tuples, tuple[0] is the gene id, tuple[1] is the KOid,  [2] is ntseq, [3] is aaseq
    """
    df = pd.read_csv(kegg_gene_filename, delimiter='\t')
    df = df[ ~pd.isnull(df['koid']) ]
    return(  list( zip( df['kegg_gene_id'].tolist(), df['koid'].tolist(), df['ntseq'].tolist(), df['aaseq'] ) )  )

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
    :return: 3-tuple, tuple[0] is start, tuple[1] is end position, tuple[2] is the strand
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
            if 'complement' in full_text:
                strand = '-'
            else:
                strand = '+'
            full_text = full_text.replace('complement(', '')
            full_text = full_text.replace(')', '')
            start_end_merged = full_text.split('\n')[0].split(':')[-1]
            start = int( start_end_merged.split('..')[0] )
            end = int( start_end_merged.split('..')[1] )
            return start, end, strand

    return (-1, -1)

def make_kegg_gene_file_name(org_code):
    """
    Using org_code (3-character code for organisms), returns kegg gene filename
    :param str org_code: 3-character code for organisms, internal to kegg database
    :return: kegg gene filename
    """
    return org_code + '_kegg_genes.txt'

def parse_args(): # pragma: no cover
    parser = argparse.ArgumentParser(description="This preprocesses kegg_db data files. It opens the organisms table, locates the genomes that has GenBank or RefSeq entries available.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', '--org_table_file', type=str, help="The organisms table file.")
    parser.add_argument('-O', '--out_dir', type=str, help="Full path to the output directory", default='./extracted_genomes')
    parser.add_argument('-s', '--seed', type=int, default=0, help="Random seed for shuffling.")
    args = parser.parse_args()
    return args

def find_key(all_keys, ncbi_id):
    """
    Given a list of all keys, find the key corresponding to an NCBI id.
    The keys in all_keys are the full contig ids. The given ncbi_id should be
    a substring of one of these keys in the list all_keys.
    """
    for key in all_keys:
        if ncbi_id in key:
            return key
    return None

def main(): # pragma: no cover
    args = parse_args()
    org_table_filename = args.org_table_file
    out_dir = args.out_dir
    seed = args.seed

    random.seed(seed)

    print(out_dir)
    if os.path.exists(out_dir) and not os.path.isfile(out_dir):
        print('Out directory exists.')
    else:
        print('Out directory does not exist, creating...')
        subprocess.call(['mkdir', out_dir])

    # get list of all bacteria organisms
    org_code_to_ncbi_ids = read_organism_table(org_table_filename)
    print('Num of all bacteria organism:')
    print( len( org_code_to_ncbi_ids.keys() ) )

    # get all bacteria with single chromosome
    print('Num of all bacteria organism with single chromosome:')
    org_code_to_ncbi_ids_single_chr = read_organism_table_for_single_chr_bacteria(org_table_filename)
    print( len( org_code_to_ncbi_ids_single_chr.keys() ) )

    directory_with_kegg_gene_files = '/data/shared_data/KEGG_data/organisms/kegg_gene_info'
    list_bacteria_single_chr_with_existing_gene_file = []
    for org_code in org_code_to_ncbi_ids_single_chr.keys():
        gene_filename = make_kegg_gene_file_name(org_code)
        file_with_path = directory_with_kegg_gene_files + '/' + gene_filename
        if exists(file_with_path):
            list_bacteria_single_chr_with_existing_gene_file.append(org_code)

    print('Number of bacterial organisms with single chromosome and existing gene file:')
    print( len(list_bacteria_single_chr_with_existing_gene_file) )

    all_genes_and_kos = []
    selected_organisms = []
    org_code_to_gene_and_ko = {}
    for org_code in tqdm(list_bacteria_single_chr_with_existing_gene_file):
        gene_filename = make_kegg_gene_file_name(org_code)
        gene_file_with_path = directory_with_kegg_gene_files + '/' + gene_filename
        gene_and_ko_list = read_gene_id_with_kos_labeled(gene_file_with_path)
        org_code_to_gene_and_ko[org_code] = gene_and_ko_list
        all_genes_and_kos = all_genes_and_kos + gene_and_ko_list
        selected_organisms.append(org_code)
        if len(all_genes_and_kos) > 500000:
            break

    print('Number of selected organisms:')
    print(len(selected_organisms))
    print('Number of total genes in these organisms:')
    print(len(all_genes_and_kos))

    print('Constructing database now...')

    kegg_db_path = '/data/shared_data/KEGG_data/organisms/'
    long_fasta_filename = 'gb_ncbi_organism.fasta'

    print('Indexing the complete fasta file...')
    record_dict = SeqIO.index(kegg_db_path+long_fasta_filename, "fasta")

    print('Randomizing the organisms:')
    random.shuffle(selected_organisms)

    all_keys = list(record_dict.keys())
    org_count = 1
    total_orgs = len(selected_organisms)
    for org_code in selected_organisms:
        print('Working with organism ' + str(org_count) + '/' + str(total_orgs))

        ncbi_ids = org_code_to_ncbi_ids[org_code]
        ncbi_id = ncbi_ids[0]
        key_in_dict = find_key(all_keys, ncbi_id)

        if key_in_dict is None:
            print(org_code, ncbi_ids, ncbi_id)
            print(list(all_keys))
            print('Problem occurred here!')
            exit(-1)

        genome_dir = os.path.join(out_dir, org_code)
        subprocess.call(['rm', '-rf', genome_dir])
        subprocess.call(['mkdir', genome_dir])

        fasta_filename = os.path.join(genome_dir, org_code+'.fasta')
        f = open(fasta_filename, 'w')
        f.write('> ' + record_dict[key_in_dict].description + '\n')
        f.write(str(record_dict[key_in_dict].seq))
        f.close()

        print('Completed writing fasta file. Now handling the genes...')

        mapping_records = []
        mapping_filename = os.path.join(genome_dir, org_code+'_mapping.csv')
        for gene_name, ko_id, nt_seq, aa_seq in tqdm(org_code_to_gene_and_ko[org_code][:10]):
            start_pos, end_pos, strand = read_gene_start_and_end_positions(gene_name)
            genome_name = org_code
            contig_id = key_in_dict
            assembly_id = org_code + '_' + contig_id
            gene_name = gene_name
            protein_id = gene_name
            mapping_records.append( (genome_name, assembly_id, gene_name, protein_id, contig_id, start_pos, end_pos, strand, aa_seq, nt_seq) )

        df = pd.DataFrame(mapping_records, columns=['genome_name', 'assembly_id', 'gene_name', 'protein_id', 'contig_id', 'start_position', 'end_position', 'strand', 'aa_sequence', 'nt_sequence'])
        df.to_csv(mapping_filename)

if __name__ == '__main__': # pragma: no cover
    main()
