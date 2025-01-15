import csv
import re

# Input and output file paths
input_file = "mutated_sequences.fasta"  # Replace with your FASTA filename
output_file = "uncertain_sequences.txt"

# Initialize variables
sequences = []
current_sequence = ""
matching_header = False

# Define target pathogenicity and review status
target_pathogenicity = "uncertain significance"
target_review_status = "criteria provided, multiple submitters, no conflicts"

# Read and parse the FASTA file
with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):  # Header line
            # Save the previous sequence if it matches the criteria
            if matching_header and current_sequence:
                sequences.append(current_sequence)
            # Reset current sequence and matching header flag
            current_sequence = ""
            matching_header = False

            # Parse the header
            header_line = line.lstrip(">")
            fields = header_line.split("_", 2)  # Split into at most 3 parts
            if len(fields) < 3:
                # Handle malformed rows
                print(f"Malformed header detected and skipped: {header_line}")
                continue  # Skip this sequence

            pathogenicity = fields[1].strip().lower()
            review_status = fields[2].strip().lower()

            # Check if the sequence matches the target pathogenicity and review status
            if pathogenicity == target_pathogenicity and review_status == target_review_status:
                matching_header = True

        else:
            if matching_header:
                current_sequence += line  # Build the sequence

    # Add the last sequence if it matches the criteria
    if matching_header and current_sequence:
        sequences.append(current_sequence)

# Write the sequences to the output file
with open(output_file, 'w') as outfile:
    for sequence in sequences:
        outfile.write(sequence + "\n")

print(f"Sequences saved to: {output_file}")
