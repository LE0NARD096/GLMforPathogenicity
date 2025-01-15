import csv

# Read and parse the FASTA file
sequence_to_header = {}

with open('path_to/mutated_sequences_BRCA1.fasta', 'r') as fasta_file:
    header = ''
    sequence_lines = []
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            if header and sequence_lines:
                # Process previous sequence
                sequence = ''.join(sequence_lines)
                sequence_to_header[sequence] = {
                    'accession_code': accession_code,
                    'pathogenicity': pathogenicity,
                    'review_status': review_status
                }
                sequence_lines = []
            # Process new header
            header = line[1:]  # Remove '>'
            # Split header into three parts
            header_parts = header.split('_', 2)
            accession_code = header_parts[0]
            pathogenicity = header_parts[1] if len(header_parts) > 1 else ''
            review_status = header_parts[2] if len(header_parts) > 2 else ''
        else:
            # Sequence line
            sequence_lines.append(line)
    # Process last sequence
    if header and sequence_lines:
        sequence = ''.join(sequence_lines)
        sequence_to_header[sequence] = {
            'accession_code': accession_code,
            'pathogenicity': pathogenicity,
            'review_status': review_status
        }

# Read and parse the predictions CSV file (without 'true_label')
sequence_data_list = []

with open('path_to/labeled_uncertain_sequences.csv', 'r') as csv_file:
    reader = csv.DictReader(csv_file)
    for row in reader:
        # Extract fields
        sequence = row['sequence']
        predicted_class = row['predicted_class']
        probability_scores_str = row['probability_scores']
        # Parse probability_scores string to get two values
        probability_scores = probability_scores_str.strip('[]').split(',')
        if len(probability_scores) == 2:
            probability_0 = probability_scores[0].strip()
            probability_1 = probability_scores[1].strip()
        else:
            probability_0 = ''
            probability_1 = ''
        # Store the data
        sequence_data_list.append({
            'sequence': sequence,
            'predicted_class': predicted_class,
            'probability_0': probability_0,
            'probability_1': probability_1
        })

# Merge datasets by sequence
output_data = []

for seq_data in sequence_data_list:
    sequence = seq_data['sequence']
    # Retrieve header information from the FASTA data
    header_info = sequence_to_header.get(sequence, {
        'accession_code': '',
        'pathogenicity': '',
        'review_status': ''
    })
    # Prepare output row
    output_row = {
        'sequence': sequence,
        'predicted_class': seq_data['predicted_class'],
        'probability_0': seq_data['probability_0'],
        'probability_1': seq_data['probability_1'],
        'accession_code': header_info['accession_code'],
        'pathogenicity': header_info['pathogenicity'],
        'review_status': header_info['review_status']
    }
    output_data.append(output_row)

# Write the output to a CSV file
output_columns = ['sequence', 'predicted_class', 'probability_0', 'probability_1',
                  'accession_code', 'pathogenicity', 'review_status']

with open('path_to/output_uncertain_BRCA1.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=output_columns)
    writer.writeheader()
    for row in output_data:
        writer.writerow(row)
