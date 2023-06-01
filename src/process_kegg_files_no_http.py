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
import time

def org_code_to_T_number(org_table_filename):
    """
    Takes as argument the organism table filename with full path.
    Returns a dictionary indexed by 3-character organism code, mapping to a list
    of T-numbers. T-numbers are unique identifiers of organisms, that is used in
    KEGG database,
    :param str org_table_filename: full path to organism filename
    :return: dictionary indexed by 3-character organism code, mapping to a list
    of T-numbers.
    :rtype: dic
    """
    # read dataframe
    df = pd.read_csv(org_table_filename, delimiter='\t')

    # keep only bacteria
    df['lineage'] = df['lineage'].str.lower()
    df = df[ df['lineage'].str.contains('bacteria') ]

    # take interesting things
    org_codes = df['org_code'].tolist()
    t_numbers = df['T_number'].tolist()

    # keep only non-empty ncbi things here
    ret_dic = {}
    for org_code, t_number in list( zip(org_codes, t_numbers) ):
        if org_code is None:
            continue
        ret_dic[org_code] = t_number

    return ret_dic

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

    # get list of all bacterial org_code to t_number mapping
    org_code_to_t_number = org_code_to_T_number(org_table_filename)

    # get all bacteria with single chromosome
    print('Num of all bacteria organism with single chromosome:')
    org_code_to_ncbi_ids_single_chr = read_organism_table_for_single_chr_bacteria(org_table_filename)
    print( len( org_code_to_ncbi_ids_single_chr.keys() ) )

    # for these organisms, identify those with gene info available
    # TODO: change here: these are hardcoded!
    directory_with_kegg_gene_files = '/scratch/shared_data/archive/KEGG_data/organisms/kegg_gene_info'
    directory_with_kegg_kff_files = '/scratch/mbr5797/all_kff_files_from_kegg'
    list_bacteria_single_chr_with_existing_gene_file = []
    for org_code in org_code_to_ncbi_ids_single_chr.keys():
        gene_filename = make_kegg_gene_file_name(org_code)
        file_with_path = directory_with_kegg_gene_files + '/' + gene_filename
        if not exists(file_with_path):
            print(f'{file_with_path} was not found!')
            continue
        t_number = org_code_to_t_number[org_code]
        kff_filename = t_number + '.kff'
        file_with_path = directory_with_kegg_kff_files + '/' + kff_filename
        if not exists(file_with_path):
            print(f'{file_with_path} was not found!')
            continue
        list_bacteria_single_chr_with_existing_gene_file.append(org_code)

    # Now, list_bacteria_single_chr_with_existing_gene_file has org codes of all bacteria
    # that we want to process!

    print('Number of bacterial organisms with single chromosome and existing gene & kff file:')
    print( len(list_bacteria_single_chr_with_existing_gene_file) )

    # now, select only a few so that we have 500K genes
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

    print( f'Num of total organisms: {len(selected_organisms)}' )

    print('Number of selected organisms:')
    print(len(selected_organisms))
    print('Number of total genes in these organisms:')
    print(len(all_genes_and_kos))

    print('Constructing database now...')
    # TODO: change things here -- these are hardcoded!
    kegg_db_path = '/scratch/shared_data/archive/KEGG_data/organisms/'
    long_fasta_filename = 'gb_ncbi_organism.fasta'

    print('Indexing the complete fasta file...')
    start_time = time.time()
    record_dict = SeqIO.index(kegg_db_path+long_fasta_filename, "fasta")
    end_time = time.time()
    print(f'Indexing complete. It took {end_time-start_time} seconds')

    print('Randomizing the organisms...')
    random.shuffle(selected_organisms)

    all_keys = list(record_dict.keys())
    num_problematic_genes = 0
    for org_code in tqdm(selected_organisms):
        ncbi_ids = org_code_to_ncbi_ids[org_code]
        ncbi_id = ncbi_ids[0]
        key_in_dict = find_key(all_keys, ncbi_id)

        if key_in_dict is None:
            print(org_code, ncbi_ids, ncbi_id)
            print(list(all_keys))
            print('Problem occurred here!')
            exit(-1)

        # prepare a directory with org_code
        genome_dir = os.path.join(out_dir, org_code)
        subprocess.call(['rm', '-rf', genome_dir])
        subprocess.call(['mkdir', genome_dir])

        # write the fasta file here
        fasta_filename = os.path.join(genome_dir, org_code+'.fasta')
        f = open(fasta_filename, 'w')
        f.write('> ' + record_dict[key_in_dict].description + '\n')
        f.write(str(record_dict[key_in_dict].seq))
        f.close()

        # list for keeping mapping records, and the mapping filename
        mapping_records = []
        mapping_filename = os.path.join(genome_dir, org_code+'_mapping.csv')

        # find kff file for this org_code
        t_number = org_code_to_t_number[org_code]
        kff_filename = t_number + '.kff'
        file_with_path = directory_with_kegg_kff_files + '/' + kff_filename
        #read in a dataframe
        df = pd.read_csv(file_with_path, delimiter='\t', header = None)
        gene_ids = df[0].tolist()
        position_strings = df[4].tolist()

        #indexed by gene_id -> maps to positions info
        gene_id_to_position_string = {}
        for gene_id, position_string in list( zip(gene_ids, position_strings) ):
            gene_id_to_position_string[gene_id] = position_string

        f = open('problematic_genes.log')
        for gene_name, ko_id, nt_seq, aa_seq in org_code_to_gene_and_ko[org_code]:
            try:
                # get the position string from the dict using gene_name/gene_id
                position_string = gene_id_to_position_string[gene_name]
                # extract the positions
                if 'complement' in position_string:
                    strand = '-'
                else:
                    strand = '+'
                full_text = position_string.replace('complement(', '')
                full_text = full_text.replace(')', '')
                full_text = full_text.replace('<', '')
                full_text = full_text.replace('>', '')
                start_end_merged = full_text.split('\n')[0].split(':')[-1]
                start_pos = int( start_end_merged.split('..')[0] )
                end_pos = int( start_end_merged.split('..')[1] )
            # record how many times this fails
            except:
                f.write(f'{gene_name}\t{position_string}\n')
                num_problematic_genes += 1
                continue
                #print('Problem with gene: ' + str(gene_name))

            # list other required info for the mapping record
            genome_name = org_code
            contig_id = key_in_dict
            assembly_id = org_code + '_' + contig_id
            gene_name = gene_name
            protein_id = gene_name

            # add to the records
            mapping_records.append( (genome_name, assembly_id, gene_name, protein_id, contig_id, start_pos, end_pos, strand, aa_seq, nt_seq) )

        # write to file
        df = pd.DataFrame(mapping_records, columns=['genome_name', 'assembly_id', 'gene_name', 'protein_id', 'contig_id', 'start_position', 'end_position', 'strand', 'aa_sequence', 'nt_sequence'])
        df.to_csv(mapping_filename)

    print('Found problems in ' + str(num_problematic_genes) + ' genes.')

if __name__ == '__main__': # pragma: no cover
    main()
