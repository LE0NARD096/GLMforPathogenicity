import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
from torch.nn import functional as F
import csv

# Load the fine-tuned model and tokenizer
tokenizer = AutoTokenizer.from_pretrained("path_to_model_folder", trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained("path_to_model_folder", trust_remote_code=True)

# Load the dataset (assuming a CSV file with a 'sequence' column)
input_file = 'path_to_input'  # Update this path
output_file = 'path_to_output'         # Update this path

# Read sequences from the input CSV file
dataset = []
with open(input_file, 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        dataset.append(row)

# Extract sequences
dna_sequences = [row['sequence'] for row in dataset]

# Tokenize the input sequences
inputs = tokenizer(dna_sequences, return_tensors='pt', padding=True, truncation=True, max_length=512)

# Move model and inputs to GPU if available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
model.to(device)
inputs = {key: val.to(device) for key, val in inputs.items()}

# Set model to evaluation mode
model.eval()

# Forward pass through the model to get the logits
with torch.no_grad():
    outputs = model(**inputs)
    logits = outputs.logits

# Apply softmax to get the probabilities
probabilities = F.softmax(logits, dim=1)

# Prepare list to store predictions
predictions = []

# Process each sequence and save prediction
for i, dna_sequence in enumerate(dna_sequences):
    # Get the predicted class by taking the index of the highest probability
    predicted_class_idx = torch.argmax(probabilities[i]).item()
    predicted_class_prob = probabilities[i, predicted_class_idx].item()

    # Assuming binary classification (0: non-functional, 1: functional)
    class_labels = {0: "non-functional", 1: "functional"}
    predicted_class = class_labels[predicted_class_idx]

    # Prepare the output dictionary
    output_entry = {
        'index': i,
        'sequence': dna_sequence,
        'predicted_class': predicted_class_idx,
        'predicted_label': predicted_class,
        'probability_scores': probabilities[i].tolist(),
    }

    predictions.append(output_entry)

# Write predictions to a CSV file
with open(output_file, 'w', newline='') as csvfile:
    fieldnames = ['index', 'sequence', 'predicted_class', 'predicted_label', 'probability_scores']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in predictions:
        writer.writerow(row)
