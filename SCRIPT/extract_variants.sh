#!/bin/bash

# Input and output files
ID_FILE="clinvar_ids_BRCA2.txt"
LOG_FILE="process_log_BRCA.txt"
OUTPUT_FILE="filtered_variants_BRCA2.txt"

# Initialize log and output files
echo "Starting extraction process on $(date)" > "$LOG_FILE"
echo -e "Accession\tSPDI\tPathogenicity\tReviewStatus" > "$OUTPUT_FILE"

# Loop through each variant ID in the input file
while read -r ID; do
    echo "Processing ID: $ID" | tee -a "$LOG_FILE"

    # Extract SPDI, pathogenicity, and review status
    esummary -db clinvar -id "$ID" -mode xml | \
    xtract -pattern DocumentSummary \
           -element accession \
           -group variation_set \
             -block variation \
               -sep ";" -element canonical_spdi \
           -group germline_classification \
             -element description \
             -element review_status \
           -group clinical_significance \
             -element description \
             -element review_status | \
    while IFS=$'\t' read -r accession spdi germ_desc germ_review clin_desc clin_review; do
        # Prefer germline_classification if available, else use clinical_significance
        pathogenicity="${germ_desc:-$clin_desc}"
        review_status="${germ_review:-$clin_review}"

        # Check if pathogenicity annotation is present
        if [[ -n "$pathogenicity" ]]; then
            echo -e "$accession\t$spdi\t$pathogenicity\t$review_status" >> "$OUTPUT_FILE"
            echo "ID $ID has pathogenicity annotation." | tee -a "$LOG_FILE"
        else
            echo "ID $ID does not have a pathogenicity annotation." | tee -a "$LOG_FILE"
        fi
    done

    if [ $? -eq 0 ]; then
        echo "Successfully processed ID: $ID" | tee -a "$LOG_FILE"
    else
        echo "Error processing ID: $ID" | tee -a "$LOG_FILE"
    fi

    # Optional: Add a delay to avoid rate limits
    sleep 0.5
done < "$ID_FILE"

echo "Extraction completed on $(date)" >> "$LOG_FILE"