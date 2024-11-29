#!/bin/bash

# Input and output files
INPUT_FILE="filtered_variants_BRCA2.txt"
SEQUENCE_OUTPUT="mutated_sequences_BRCA2.fasta"

# Remove any existing output file
rm -f "$SEQUENCE_OUTPUT"

# Loop through each line in the input file, skipping the header
tail -n +2 "$INPUT_FILE" | while read -r line; do
    # Read the line into an array using tab as delimiter
    IFS=$'\t' read -r -a fields <<< "$line"

    # Check that we have at least 4 fields
    if [[ ${#fields[@]} -lt 4 ]]; then
        echo "Warning: Not enough fields in line: $line. Skipping this entry." >&2
        continue
    fi

    # Assign variables
    accession="${fields[0]}"
    spdi="${fields[1]}"
    pathogenicity="${fields[2]}"
    review_status="${fields[3]}"

    # Check if SPDI is present and valid
    if [[ -z "$spdi" ]]; then
        echo "Warning: SPDI is missing for $accession. Skipping this entry." >&2
        continue
    fi

    # Parse SPDI components
    chrom=$(echo "$spdi" | cut -d':' -f1)       # Chromosome/Accession
    pos=$(echo "$spdi" | cut -d':' -f2)         # Position (zero-based)
    ref=$(echo "$spdi" | cut -d':' -f3)         # Reference allele
    alt=$(echo "$spdi" | cut -d':' -f4)         # Alternate allele

    # Check for incomplete SPDI
    if [[ -z "$chrom" || -z "$pos" ]]; then
        echo "Warning: Incomplete SPDI for $accession. Skipping this entry." >&2
        continue
    fi

    # Check if position is a number
    if ! [[ "$pos" =~ ^[0-9]+$ ]]; then
        echo "Warning: Invalid position in SPDI for $accession: $pos. Skipping this entry." >&2
        continue
    fi

    # Convert SPDI position to zero-based index (it already is)
    zero_based_pos=$pos

    # Calculate ref_length and alt_length
    ref_length=${#ref}
    alt_length=${#alt}

    # Calculate the change in length
    delta_length=$((alt_length - ref_length))

    # Calculate the desired fetched sequence length
    fetched_seq_length=$((300 - delta_length))

    # Ensure fetched_seq_length is at least 1
    if [[ "$fetched_seq_length" -lt 1 ]]; then
        echo "Error: Fetched sequence length is less than 1 for $accession" >&2
        continue
    fi

    # Calculate left and right context lengths
    left_context_length=$(( (fetched_seq_length - ref_length) / 2 ))
    right_context_length=$(( fetched_seq_length - ref_length - left_context_length ))

    # Calculate start and end positions
    start=$(( zero_based_pos - left_context_length ))
    end=$(( zero_based_pos + ref_length + right_context_length - 1 ))

    # Adjust start position to be at least 0
    if [[ "$start" -lt 0 ]]; then
        start=0
    fi

    # Fetch the reference sequence
    ref_seq=$(efetch -db nucleotide -id "$chrom" -seq_start $((start + 1)) -seq_stop $((end + 1)) -format fasta | sed '1d' | tr -d '\n')

    # Ensure we fetched a sequence
    if [[ -z "$ref_seq" ]]; then
        echo "Error: Unable to fetch sequence for $chrom from positions $((start + 1)) to $((end + 1)) for $accession" >&2
        continue
    fi

    # Ensure the fetched sequence length matches expected
    expected_seq_length=$((end - start + 1))
    seq_length=${#ref_seq}
    if [[ "$seq_length" -ne "$expected_seq_length" ]]; then
        echo "Error: Fetched sequence length ($seq_length) does not match expected length ($expected_seq_length) for $accession" >&2
        continue
    fi

    # Calculate mutation position relative to the fetched sequence
    mutation_pos=$((zero_based_pos - start))

    # Ensure mutation position is valid
    if [[ "$mutation_pos" -lt 0 || "$((mutation_pos + ref_length - 1))" -ge "$seq_length" ]]; then
        echo "Error: Invalid mutation position for $accession" >&2
        continue
    fi

    # Extract the reference allele from the sequence
    ref_base=${ref_seq:$mutation_pos:$ref_length}

    # Compare the extracted reference allele with the expected one
    if [[ "$ref_base" != "$ref" ]]; then
        echo "Error: Reference base(s) at position do not match expected allele for $accession" >&2
        continue
    fi

    # Apply mutation
    mutated_seq="${ref_seq:0:$mutation_pos}${alt}${ref_seq:$((mutation_pos + ref_length))}"

    # Ensure the mutated sequence is 300 bp long
    mutated_seq_length=${#mutated_seq}
    if [[ "$mutated_seq_length" -ne 300 ]]; then
        echo "Error: Mutated sequence length is not 300 bp for $accession" >&2
        continue
    fi

    # Save sequence to output with header containing accession, pathogenicity, and review status
    echo -e ">${accession}_${pathogenicity}_${review_status}\n$mutated_seq" >> "$SEQUENCE_OUTPUT"
done
