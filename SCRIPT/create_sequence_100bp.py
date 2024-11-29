import pandas as pd

# Input and output file paths
input_file = '/NFSHOME/lmasci/DNABERT_2/DATA/labeled_sequences_BRCA1.csv'
output_file = '/NFSHOME/lmasci/DNABERT_2/DATA/labeled_100bp_BRCA1.csv'

# Set the desired chunk size (e.g., 100 base pairs)
chunk_size = 100

# Load the input CSV file into a DataFrame
df = pd.read_csv(input_file)

# Prepare an empty list to collect output data
output_data = []

# Process each row in the DataFrame
for index, row in df.iterrows():
    sequence = row['sequence']
    label = int(row['label'])
    
    # Split the sequence into chunks of the specified size
    chunks = [sequence[i:i+chunk_size] for i in range(0, len(sequence), chunk_size)]
    
    # Initialize chunk labels based on the original label
    chunk_labels = [0] * len(chunks)
    if label == 0:
        # All chunks remain labeled as 0
        pass
    elif label == 1:
        # First and last chunks are labeled 0, middle chunks may be labeled 1
        for idx, chunk in enumerate(chunks):
            start_pos = idx * chunk_size
            end_pos = start_pos + len(chunk) - 1
            # Label middle chunks (from position 100 to 199) as 1
            if start_pos >= 100 and end_pos <= 199:
                chunk_labels[idx] = 1
    else:
        print(f"Invalid label (must be 0 or 1): {label}", file=sys.stderr)
        continue
    
    # Append the chunks with their corresponding labels to the output data
    for chunk, chunk_label in zip(chunks, chunk_labels):
        output_data.append({'sequence': chunk, 'label': chunk_label})

# Convert output data to a DataFrame
output_df = pd.DataFrame(output_data)

# Write the output DataFrame to a new CSV file
output_df.to_csv(output_file, index=False)
