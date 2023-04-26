# Variant Comparison Tool

This Python script compares variants from a VCF file against a list of known mutations provided in a TSV file. It generates an output TSV file containing the original TSV data with additional 'found_vaf' and 'detected' columns, and a plot comparing Expected VAF and Found VAF.


This tool assumes that you have a `known_variants.tsv` such as:
|chr  |position |ref             |alt|gene  |hgvsp               |vaf  |
|-----|---------|----------------|---|------|--------------------|-----|
|chr1 |115256530|G               |T  |NRAS  |p.Q61K              |12.5%|
|chr3 |41266101 |C               |A  |CTNNB1|p.S33Y              |32.5%|
|chr3 |41266133 |CCTT            |C  |CTNNB1|p.S45del            |10.0%|
|chr3 |178936091|G               |A  |PIK3CA|p.E545K             |9.0% |
|chr3 |178952085|A               |G  |PIK3CA|p.H1047R            |17.5%|
|chr4 |55599321 |A               |T  |KIT   |p.D816V             |10.0%|
|chr4 |55602765 |G               |C  |KIT   |p.L862L             |7.5% |
|chr5 |112175770|G               |A  |APC   |p.T1493T            |35.0%|
|chr7 |55241707 |G               |A  |EGFR  |p.G719S             |24.5%|

And an `input.vcf` to be tested whether if the variants have beed detected or not.

## Installation

1. Install the required Python packages:

```
pip3 install pandas seaborn matplotlib numpy
```

2. Clone the repository
```
git clone https://github.com/GENCARDIO/eval_somatic_control.git
cd eval_somatic_control
```

### Usage

To run the script, use the following command:

```
python eval_somatic_control.py --input_vcf sample.vcf --known_tsv known_variants.tsv --output_tsv output_file.tsv

Arguments

    --input_vcf: Input VCF file containing the query variants.
    --known_tsv: Input TSV file containing the list of known mutations.
    --output_tsv: Output TSV file with additional columns.
```

