import csv

# Input and output file paths
input_file = "path_to/mutated_sequences_BRCA2_500bp.fasta"
output_file = "path_to/labeled_sequences_BRCA2_500bp.csv"

# Initialize variables
sequences = []
current_sequence = ""
label = None

# Define acceptable pathogenicity labels
acceptable_pathogenicity = {
    "pathogenic": 1,
    "likely pathogenic": 1,
    "benign": 0,
    "likely benign": 0
}

# Define pathogenicity labels to exclude
exclude_pathogenicity = [
    "uncertain significance",
    "conflicting classifications of pathogenicity",
    # Add other labels to exclude as needed
]

# Define review statuses to exclude
exclude_review_status = [
    "no assertion criteria provided",
    "criteria provided, single submitter",
    # Add other review statuses to exclude as needed
]

# Read and parse the FASTA file
with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):  # Header line
            # Save the previous sequence if label is not None
            if current_sequence and label is not None:
                sequences.append((current_sequence, label))
            # Reset current sequence and label
            current_sequence = ""
            label = None

            # Parse the header
            header_line = line.lstrip(">")
            fields = header_line.split("_", 2)  # Split into at most 3 parts
            if len(fields) < 3:
                # Handle malformed rows
                print(f"Malformed header detected and skipped: {header_line}")
                continue  # Skip this sequence

            accession = fields[0]
            pathogenicity = fields[1].strip().lower()
            review_status = fields[2].strip().lower()

            # Exclude sequences based on review status
            exclude = False
            for exclude_status in exclude_review_status:
                if exclude_status.lower() == review_status:
                    exclude = True
                    break
            if exclude:
                continue  # Skip this sequence

            # Split pathogenicity on '/'
            pathogenicity_labels = [label.strip() for label in pathogenicity.split('/')]

            # Exclude sequences based on pathogenicity
            exclude = False
            for label_str in pathogenicity_labels:
                if label_str in exclude_pathogenicity:
                    exclude = True
                    break
            if exclude:
                continue  # Skip this sequence

            # Determine the label based on pathogenicity
            labels_set = set()
            for label_str in pathogenicity_labels:
                label_value = acceptable_pathogenicity.get(label_str)
                if label_value is not None:
                    labels_set.add(label_value)

            if not labels_set:
                # No acceptable labels found, skip this sequence
                continue
            elif len(labels_set) == 1:
                # All labels agree
                label = labels_set.pop()
            else:
                # Conflicting labels, skip the sequence
                continue

        else:
            current_sequence += line  # Build the sequence

    # Add the last sequence if label is not None
    if current_sequence and label is not None:
        sequences.append((current_sequence, label))

# Write the sequences and labels to a CSV file
with open(output_file, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["sequence", "label"])  # Header
    writer.writerows(sequences)  # Write each sequence and its label

print(f"CSV file saved as: {output_file}")
