import pytest
import src.process_kegg_files as proc

def test_read_gene_ids():
    organism_filename = 'data/organism_table.txt'
    sample_gene_filename = 'data/sce_kegg_genes.txt'

    gene_list = proc.read_gene_ids(sample_gene_filename)
    gene_list_with_kos = proc.read_gene_id_with_ko(sample_gene_filename)

    assert len(gene_list) >= len(gene_list_with_kos)

def test_read_organism_table():
    organism_filename = 'data/organism_table.txt'
    org_to_contig_id_dic = proc.read_organism_table(organism_filename)
    org_to_contig_id_dic_single_chr = proc.read_organism_table_for_single_chr_organisms(organism_filename)
    assert len( org_to_contig_id_dic.keys() ) >= len( org_to_contig_id_dic_single_chr.keys() )

def test_read_gene_id_with_kos_labeled():
    sample_gene_filename = 'data/sce_kegg_genes.txt'
    gene_id_ko_list = proc.read_gene_id_with_kos_labeled(sample_gene_filename)
    for gene_id, ko_id in gene_id_ko_list:
        assert ko_id is not None

def test_read_gene_start_and_end_positions():
    sample_gene_id = 'sce:YAL062W'
    start_correct = 31567
    end_correct = 32940
    start, end = proc.read_gene_start_and_end_positions(sample_gene_id)
    assert start == start_correct
    assert end == end_correct

def test_read_gene_start_and_end_positions_2():
    sample_gene_id = 'sce:YAL068C'
    start_correct = 1807
    end_correct = 2169
    start, end = proc.read_gene_start_and_end_positions(sample_gene_id)
    assert start == start_correct
    assert end == end_correct

def test_read_gene_start_and_end_positions_3():
    sample_gene_id = 'eco:b0001'
    start_correct = 190
    end_correct = 255
    start, end = proc.read_gene_start_and_end_positions(sample_gene_id)
    assert start == start_correct
    assert end == end_correct

def test_make_kegg_gene_file_name():
    org_code = 'sce'
    gene_fname = 'sce_kegg_genes.txt'
    assert gene_fname == proc.make_kegg_gene_file_name(org_code)

def test_find_key_1():
    lst = ['string1.1', 'string2.3']
    key = 'string1'
    result = 'string1.1'
    assert result == proc.find_key(lst, key)
    key = 'string2'
    result = 'string2.3'
    assert result == proc.find_key(lst, key)
    key = 'string5'
    result = None
    assert result == proc.find_key(lst, key)
