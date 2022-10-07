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
