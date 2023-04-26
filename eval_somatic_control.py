import os
import sys
import argparse
import csv
import gzip
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd

def plot_expected_vaf_vs_found_vaf(data_list, output_png):
    """
    Plots the Expected VAF against the Found VAF using Seaborn and saves the plot as a PNG file.

    :param data_list: List of dictionaries containing 'vaf' and 'found_vaf' keys.
    :param output_png: Output file name to save the plot as a PNG file.
    """
    vaf = [d['vaf'] for d in data_list]
    found_vaf = [d['found_vaf'] for d in data_list]

    # Convert data list to pandas DataFrame
    data_df = pd.DataFrame(data_list)

    # Calculate R-squared
    correlation_matrix = np.corrcoef(data_df['vaf'], data_df['found_vaf'])
    r_squared = correlation_matrix[0, 1] ** 2

    # Create the Seaborn plot with a regression line
    sns.set(style="darkgrid")
    sns.regplot(x="found_vaf", y="vaf", data=data_df)

    # Set axis labels
    plt.xlabel('Found VAF')
    plt.ylabel('Expected VAF')
    plt.title('Expected VAF vs Found VAF')
    # Add R-squared annotation
    plt.text(min(data_df['found_vaf']), max(data_df['vaf']), f'R-squared: {r_squared:.2f}', horizontalalignment='left', verticalalignment='top')

    # output the plot
    plt.savefig(output_png, dpi=130)


def check_variants_presence(tsv_data, vcf_data, output_csv):
    """
    Compares the variants from the TSV data and VCF data, and outputs a TSV file containing
    the original TSV data with additional 'found_vaf' and 'detected' columns. Also generates a
    plot comparing Expected VAF and Found VAF.

    :param tsv_data: List of dictionaries containing the TSV data.
    :param vcf_data: List of dictionaries containing the VCF data.
    :param output_csv: Output file name for the TSV file with additional columns.
    """
    results = []

    header_keys = list(tsv_data[0].keys())
    header_keys.append("found_vaf")
    header_keys.append("detected")

    for tsv_variant in tsv_data:
        output_variant = tsv_variant
        output_variant['detected'] = False
        output_variant["found_vaf"] = "."
        for vcf_variant in vcf_data:
            if (tsv_variant["chr"] == vcf_variant["chr"]
                    and tsv_variant["position"] == vcf_variant["pos"]
                    and tsv_variant["ref"] == vcf_variant["ref"]
                    and tsv_variant["alt"] == vcf_variant["alt"]):
                output_variant['detected'] = True
                output_variant["found_vaf"] = vcf_variant["vaf"]
                break

        results.append(output_variant)

    o = open(output_csv, "w")
    o.write('\t'.join(header_keys)+"\n")
    for item in results:
        out_str = '\t'.join(str(item[key]) for key in header_keys)
        o.write(out_str+"\n")
    o.close()

    output_dir = os.path.dirname(output_csv)
    output_name = os.path.basename(output_csv).replace(".tsv", ".png")

    output_png = os.path.join(output_dir, output_name)

    plot_expected_vaf_vs_found_vaf(results, output_png)


def read_vcf_data(input_vcf):
    """
    Reads the VCF file and returns a list of dictionaries containing the VCF data.

    :param input_vcf: Input VCF file.
    :return: List of dictionaries containing the VCF data.
    """
    data_list = []

    open_func = gzip.open if input_vcf.endswith(".gz") else open

    with open_func(input_vcf, "rt") as vcf_file:
        for line in vcf_file:
            if not line.startswith("#"):  # Skip header lines
                fields = line.strip().split("\t")

                # Find the index of the AD field in the FORMAT column
                format_fields = fields[8].split(":")
                ad_index = format_fields.index("AD") if "AD" in format_fields else None

                # Get the corresponding AD values from the sample column
                if ad_index is not None:
                    sample_fields = fields[9].split(":")
                    ref_depth, alt_depth = map(int, sample_fields[ad_index].split(","))
                    vaf = round(alt_depth / (ref_depth + alt_depth) * 100, 3)
                else:
                    vaf = None

                data_dict = {
                    "chr": fields[0],
                    "pos": int(fields[1]),
                    "ref": fields[3],
                    "alt": fields[4],
                    "vaf": vaf
                }
                data_list.append(data_dict)
    return data_list

def read_known_tsv(input_tsv: str):
    """
    Reads the TSV file containing known variants and returns a list of dictionaries
    with the TSV data.

    :param input_tsv: Input TSV file.
    :return: List of dictionaries containing the TSV data.
    """    
    data_list = []

    with open(input_tsv, "r") as tsv_file:
        tsv_reader = csv.reader(tsv_file, delimiter="\t")
        
        header = next(tsv_reader)  # Read header row
        required_fields = ["chr", "position", "ref", "alt", "gene", "hgvsp", "vaf"]

        # Check if all required fields are present in the header
        if not all(field in header for field in required_fields):
            raise ValueError("Invalid header. Missing required fields.")

        for row in tsv_reader:
            data_dict = {
                "chr": row[0].strip(),
                "position": int(row[1].strip()),
                "ref": row[2].strip(),
                "alt": row[3].strip(),
                "gene": row[4].strip(),
                "hgvsp": row[5].strip(),
                "vaf": float(row[6].strip().replace("%", ""))
            }
            data_list.append(data_dict)

    return data_list


def main():
    """
    Main function that processes command-line arguments and calls appropriate functions
    to read, process, and compare VCF and TSV data.
    """
    parser = argparse.ArgumentParser(description="Evaluate a VCF file agains a list of known mutations")
    parser.add_argument("--input_vcf", dest="input_vcf", required=True, help="Query variants in VCF format")
    parser.add_argument("--known_tsv", dest="known_tsv", required=True, help="Known variants in TSV format")
    parser.add_argument("--output_tsv", dest="output_tsv", required=True, help="Output results in TSV format")

    args = parser.parse_args()

    known_csv = args.known_tsv
    sample_vcf = args.input_vcf
    output_tsv = args.output_tsv
    
    sample_variants = read_vcf_data(sample_vcf)
    control_variants = read_known_tsv(known_csv)
    check_variants_presence(control_variants, sample_variants, output_tsv)


if __name__ == "__main__":
    main()
