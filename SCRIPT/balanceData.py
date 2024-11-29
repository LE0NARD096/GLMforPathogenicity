import pandas as pd

# Define a function to split the sequence into 300bp chunks with 150bp overlaps
def split_sequence(sequence, window_size=500, overlap=400):
    chunks = []
    for i in range(0, len(sequence) - window_size + 1, window_size - overlap):
        chunk = sequence[i:i + window_size]
        chunks.append(chunk)
    return chunks

# Main code to split the BRCA1 sequence and append the label
def main():
    # Read the BRCA1 sequence from a file (after fetching with efetch)
    with open('/NFSHOME/lmasci/DNABERT_2/DATA/BRCA1_WT.fasta', 'r') as f:
        brca1_sequence = ''.join([line.strip() for line in f.readlines() if not line.startswith('>')])

    # Split the sequence into chunks of 300bp with 150bp overlap
    split_sequences = split_sequence(brca1_sequence)

    # Create a DataFrame to store sequences and labels
    data = {
        "sequence": split_sequences,
        "label": [0] * len(split_sequences)  # Label 0 for all sequences
    }

    df = pd.DataFrame(data)

    # Save to CSV
    df.to_csv('/NFSHOME/lmasci/DNABERT_2/DATA/brca1_wild_sequences.csv', index=False)

    print("Sequences saved to _ .csv")

# Run the main function
if __name__ == "__main__":
    main()
