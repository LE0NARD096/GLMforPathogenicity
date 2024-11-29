import re
import pandas as pd
import matplotlib.pyplot as plt

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
    if len(deleted_seq) == 1 and len(inserted_seq) == 1:
        return 'SNP'
    elif len(deleted_seq) > 0 and inserted_seq == '':
        return 'Deletion'
    elif deleted_seq == '' and len(inserted_seq) > 0:
        return 'Insertion'
    else:
        return 'Indel'  # All other cases, including substitutions of length > 1

def main():
    variant_info_file = '/NFSHOME/lmasci/DNABERT_2/DATA/filtered_variants_BRCA2_USED.txt'  # Update with your file path
    accession_list_file = '/NFSHOME/lmasci/DNABERT_2/accession_code_misclassified_BRCA2.txt'         # Update with your file path

    # Read variant information
    variants = read_variant_info(variant_info_file)

    # Read accession list of misclassified variants
    misclassified_accessions = read_accession_list(accession_list_file)

    # Process variants to determine mutation types
    total_variants_data = []
    misclassified_variants_data = []

    for variant in variants:
        accession = variant['Accession']
        spdi = variant['SPDI']
        pathogenicity = variant['Pathogenicity']
        spdi_info = parse_spdi(spdi)
        if spdi_info is None:
            continue
        seq_id, genomic_position, deleted_seq, inserted_seq, mutation_type = spdi_info

        # Collect data for total variants
        total_variants_data.append({
            'Accession': accession,
            'Mutation_type': mutation_type
        })

        # Collect data for misclassified variants
        if accession in misclassified_accessions:
            misclassified_variants_data.append({
                'Accession': accession,
                'Mutation_type': mutation_type
            })

    # Create DataFrames
    df_total = pd.DataFrame(total_variants_data)
    df_misclassified = pd.DataFrame(misclassified_variants_data)

    if df_total.empty:
        print("No variants found in the filtered variants file.")
        return

    # Count total variants per mutation type
    total_counts = df_total['Mutation_type'].value_counts().rename('Total_Counts')

    # Count misclassified variants per mutation type
    misclassified_counts = df_misclassified['Mutation_type'].value_counts().rename('Misclassified_Counts')

    # Combine counts into a single DataFrame
    counts_df = pd.concat([total_counts, misclassified_counts], axis=1).fillna(0)

    # Calculate percentage of misclassified variants per mutation type
    counts_df['Misclassified_Percentage'] = (counts_df['Misclassified_Counts'] / counts_df['Total_Counts']) * 100

    print(counts_df)

    # Plotting
    counts_df = counts_df.reset_index()  # 'index' will be 'Mutation_type'
    counts_df = counts_df.rename(columns={'index': 'Mutation_type'})
    counts_df = counts_df.sort_values('Total_Counts', ascending=False)

    plt.figure(figsize=(10, 6))

    # Bar plot for total variants
    plt.bar(counts_df['Mutation_type'], counts_df['Total_Counts'], color='lightgray', label='Total Variants')

    # Bar plot for misclassified variants (overlayed)
    plt.bar(counts_df['Mutation_type'], counts_df['Misclassified_Counts'], color='red', label='Misclassified Variants')

    plt.xlabel('Mutation Type')
    plt.ylabel('Number of Variants')
    plt.title('Comparison of Total and Misclassified Variants by Mutation Type')
    plt.legend()

    plt.tight_layout()
    plt.savefig('mutation_type_comparison_BRCA2.png', format='png', dpi=300)
    plt.show()

    # Plot percentage of misclassified variants
    plt.figure(figsize=(10, 6))
    plt.bar(counts_df['Mutation_type'], counts_df['Misclassified_Percentage'], color='blue')

    plt.xlabel('Mutation Type')
    plt.ylabel('Percentage of Misclassified Variants')
    plt.title('Percentage of Misclassified Variants by Mutation Type')

    plt.tight_layout()
    plt.savefig('misclassified_percentage_BRCA2.png', format='png', dpi=300)
    plt.show()

if __name__ == '__main__':
    main()
