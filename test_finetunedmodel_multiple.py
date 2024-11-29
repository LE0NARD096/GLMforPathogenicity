import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
from torch.nn import functional as F
import csv

# Load the fine-tuned model and tokenizer
tokenizer = AutoTokenizer.from_pretrained("/NFSHOME/lmasci/DNABERT_2/finetune/output/BRCATOT_TEST", trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained("/NFSHOME/lmasci/DNABERT_2/finetune/output/BRCATOT_TEST", trust_remote_code=True)

# Load the dataset manually (assuming a CSV file with 'sequence' and 'label' columns)
dataset = []
with open('/NFSHOME/lmasci/DNABERT_2/DATA/syntheticData_BRCA2.csv', 'r') as csvfile:
    reader = csv.DictReader(csvfile)
    for i, row in enumerate(reader):
        if i >= 15000:  # Limit to 100 rows
            break
        dataset.append(row)

# Extract sequences and labels
dna_sequences = [row['sequence'] for row in dataset]
true_labels = [int(row['label']) for row in dataset]

# Tokenize the input sequences
inputs = tokenizer(dna_sequences, return_tensors='pt', padding=True, truncation=True, max_length=512)

# Forward pass through the model to get the logits
outputs = model(**inputs)
logits = outputs.logits

# Apply softmax to get the probabilities
probabilities = F.softmax(logits, dim=1)

# Prepare lists to store correctly and incorrectly labeled sequences
correctly_labeled = []
incorrectly_labeled = []

# Process each sequence, determine correctness, and save information
for i, dna_sequence in enumerate(dna_sequences):
    # Get the predicted class by taking the index of the highest probability
    predicted_class_idx = torch.argmax(probabilities[i]).item()
    predicted_class_prob = probabilities[i, predicted_class_idx].item()
    true_label = true_labels[i]

    # Assuming binary classification (0: non-functional, 1: functional), define the class labels
    class_labels = {0: "non-functional", 1: "functional"}
    predicted_class = class_labels[predicted_class_idx]

    # Prepare the output dictionary
    output_entry = {
        'index': i,
        'sequence': dna_sequence,
        'true_label': true_label,
        'predicted_class': predicted_class_idx,
        'probability_scores': probabilities[i].tolist(),
    }

    # Classify as correct or incorrect
    if true_label == predicted_class_idx:
        correctly_labeled.append(output_entry)
    else:
        incorrectly_labeled.append(output_entry)

# Write correctly labeled sequences to a CSV file
with open('/NFSHOME/lmasci/DNABERT_2/DATA/correct/TOT_BRCA2.csv', 'w', newline='') as csvfile:
    fieldnames = ['index', 'sequence', 'true_label', 'predicted_class', 'probability_scores']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in correctly_labeled:
        writer.writerow(row)

# Write incorrectly labeled sequences to a CSV file
with open('/NFSHOME/lmasci/DNABERT_2/DATA/incorrect/TOT_BRCA2.csv', 'w', newline='') as csvfile:
    fieldnames = ['index', 'sequence', 'true_label', 'predicted_class', 'probability_scores']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for row in incorrectly_labeled:
        writer.writerow(row)
