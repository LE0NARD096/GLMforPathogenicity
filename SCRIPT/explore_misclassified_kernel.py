import re
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def read_reference_gene(filename):
    with open(filename, 'r') as f:
        header = f.readline().strip()
        sequence = ''.join(line.strip() for line in f)
    # Parse the header to get the gene coordinates
    match = re.match(r'^>(\S+):(\d+)-(\d+)', header)
    if match:
        seq_id = match.group(1)
        gene_start = int(match.group(2))
        gene_end = int(match.group(3))
    else:
        raise ValueError("Invalid reference gene header format.")
    return seq_id, gene_start, gene_end, sequence

def read_variant_info(filename):
    variant_list = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line == '':
                continue
            fields = line.split('\t')
            if len(fields) < 4:
                continue
            accession, spdi, pathogenicity, review_status = fields[:4]
            variant_list.append({
                'Accession': accession,
                'SPDI': spdi,
                'Pathogenicity': pathogenicity,
                'ReviewStatus': review_status
            })
    return variant_list

def read_accession_list(filename):
    accession_set = set()
    with open(filename, 'r') as f:
        for line in f:
            accession = line.strip()
            if accession == '' or accession.lower() == 'accession_code':
                continue
            accession_set.add(accession)
    return accession_set

def parse_spdi(spdi):
    fields = spdi.split(':')
    if len(fields) != 4:
        return None
    seq_id, position, deleted_seq, inserted_seq = fields
    try:
        position = int(position)
    except ValueError:
        return None
    mutation_type = determine_mutation_type(deleted_seq, inserted_seq)
    genomic_position = position + 1  # Convert zero-based to one-based
    return seq_id, genomic_position, deleted_seq, inserted_seq, mutation_type

def determine_mutation_type(deleted_seq, inserted_seq):
    if deleted_seq == inserted_seq:
        return 'No change'
    elif len(deleted_seq) == 1 and len(inserted_seq) == 1:
        return 'SNP'
    elif len(deleted_seq) > 0 and inserted_seq == '':
        return 'Deletion'
    elif deleted_seq == '' and len(inserted_seq) > 0:
        return 'Insertion'
    elif len(deleted_seq) > 0 and len(inserted_seq) > 0:
        if len(deleted_seq) == len(inserted_seq):
            return 'Substitution'
        else:
            return 'Indel'
    else:
        return 'Complex'

def main():
    reference_gene_file = '/NFSHOME/lmasci/DNABERT_2/DATA/BRCA1_WT.fasta'    # Update with your file path
    variant_info_file = '/NFSHOME/lmasci/DNABERT_2/DATA/filtered_variants_BRCA1.txt'          # Update with your file path
    accession_list_file = '/NFSHOME/lmasci/DNABERT_2/accession_code_uncertain.txt'      # Update with your file path

    # Read reference gene
    seq_id, gene_start, gene_end, gene_seq = read_reference_gene(reference_gene_file)
    gene_length = gene_end - gene_start + 1

    # Read variant information
    variants = read_variant_info(variant_info_file)

    # Read accession list
    accession_set = read_accession_list(accession_list_file)

    # Process variants
    data = []
    for variant in variants:
        accession = variant['Accession']

        # Filter variants by accession list
        if accession not in accession_set:
            continue

        spdi = variant['SPDI']
        pathogenicity = variant['Pathogenicity']
        spdi_info = parse_spdi(spdi)
        if spdi_info is None:
            continue
        spdi_seq_id, genomic_position, deleted_seq, inserted_seq, mutation_type = spdi_info

        # Ensure the variant is on the same chromosome
        if spdi_seq_id != seq_id:
            continue

        # Compute position on gene
        position_on_gene = genomic_position - gene_start + 1

        # Skip if the position is outside the gene
        if not (1 <= position_on_gene <= gene_length):
            continue

        data.append({
            'Accession': accession,
            'Pathogenicity': pathogenicity,
            'Position_on_gene': position_on_gene,
            'Mutation_type': mutation_type
        })

    # Create DataFrame
    df = pd.DataFrame(data)

    # Display the dataset
    #print(df)

    if df.empty:
        print("No variants matched the accession list.")
        return

    # Visualization - Kernel Density Plot
    plt.figure(figsize=(12, 6))

    # Plot the kernel density estimate
    sns.kdeplot(df['Position_on_gene'], fill=True, bw_adjust=0.1, color='blue')

    plt.xlim(1, gene_length)
    plt.xlabel('Position on Gene')
    plt.ylabel('Density')
    plt.title('Mutation Density along Reference Gene')

    # Save plot to file
    plt.savefig('mutation_density_on_gene_uncertain.png', format='png', dpi=300)  # Save as PNG with 300 DPI

    # Display plot
    plt.show()

if __name__ == '__main__':
    main()
