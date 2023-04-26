import pytest
from eval_somatic_control import read_vcf_data, read_known_tsv, check_variants_presence
import os
import sys
dirname = os.path.dirname(os.path.realpath(__file__))

def test_read_vcf_data():

    test_vcf = dirname + "/test.vcf"
    vcf_data = read_vcf_data(test_vcf)
    assert isinstance(vcf_data, list)
    assert len(vcf_data) > 0
    assert all(isinstance(variant, dict) for variant in vcf_data)

def test_read_known_tsv():

    tsv_data = read_known_tsv(dirname + "/known.tsv")
    assert isinstance(tsv_data, list)
    assert len(tsv_data) > 0
    assert all(isinstance(variant, dict) for variant in tsv_data)

def test_check_variants_presence(tmpdir):
    tsv_data = [
        {"chr": "chr1", "position": 1000, "ref": "A", "alt": "G", "gene": "GENE1", "hgvsp": "p.P1L", "vaf": 10},
        {"chr": "chr1", "position": 2000, "ref": "C", "alt": "T", "gene": "GENE2", "hgvsp": "p.P2L", "vaf": 15},
        {"chr": "chr1", "position": 3000, "ref": "C", "alt": "T", "gene": "GENE2", "hgvsp": "p.P2L", "vaf": 20},
        {"chr": "chr1", "position": 4000, "ref": "C", "alt": "T", "gene": "GENE2", "hgvsp": "p.P2L", "vaf": 25},
        {"chr": "chr1", "position": 5000, "ref": "C", "alt": "T", "gene": "GENE2", "hgvsp": "p.P2L", "vaf": 30},
        {"chr": "chr1", "position": 6000, "ref": "C", "alt": "T", "gene": "GENE2", "hgvsp": "p.P2L", "vaf": 75}
    ]
    
    vcf_data = [
        {"chr": "chr1", "pos": 1000, "ref": "A", "alt": "G", "vaf": 12},
        {"chr": "chr1", "pos": 2000, "ref": "C", "alt": "T", "vaf": 17},
        {"chr": "chr1", "pos": 3000, "ref": "C", "alt": "T", "vaf": 22},
        {"chr": "chr1", "pos": 4000, "ref": "C", "alt": "T", "vaf": 27},
        {"chr": "chr1", "pos": 5000, "ref": "C", "alt": "T", "vaf": 33},
        {"chr": "chr1", "pos": 6000, "ref": "C", "alt": "T", "vaf": 68}
    ]
    
    output_tsv = dirname + "/output.tsv"
    check_variants_presence(tsv_data, vcf_data, output_tsv)

    with open(output_tsv, "r") as f:
        content = f.read()

    expected_content = "chr\tposition\tref\talt\tgene\thgvsp\tvaf\tfound_vaf\tdetected\n" \
        "chr1\t1000\tA\tG\tGENE1\tp.P1L\t10\t12\tTrue\n" \
        "chr1\t2000\tC\tT\tGENE2\tp.P2L\t15\t17\tTrue\n" \
        "chr1\t3000\tC\tT\tGENE2\tp.P2L\t20\t22\tTrue\n" \
        "chr1\t4000\tC\tT\tGENE2\tp.P2L\t25\t27\tTrue\n" \
        "chr1\t5000\tC\tT\tGENE2\tp.P2L\t30\t33\tTrue\n" \
        "chr1\t6000\tC\tT\tGENE2\tp.P2L\t75\t68\tTrue\n" 
                       
    assert content == expected_content