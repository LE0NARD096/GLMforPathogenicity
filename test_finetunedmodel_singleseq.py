#import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
#from transformers import AutoModel, AutoTokenizer

from torch.nn import functional as F

# Load the fine-tuned model and tokenizer
tokenizer = AutoTokenizer.from_pretrained("path_to_model_folder", trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained("path_to_model_folder", trust_remote_code=True)
#model = AutoModel.from_pretrained("output/BRCA_model_sample_balance_third",trust_remote_code=True)

#print(model)

# Input DNA sequence
dna_sequence = "GCAAAACCCCTAATCTAAGCATAGCATTCAATTTTGGCCCTCTGTTTCTACCTAGTTCTGCTTGAATGTTTTCATCACTGGAACCTATTTCATTAATACTGGAGCCCACTTCATTAGTACTGGAACCTACTTCATTAATATTGCTTGAGTTGGCTTCTTTAAAAACATTTTCTCTAATGTTATTACGGCTAATTGTGCTCACTGTACTTGGAATGTTCTCATTTCCCATTTCTCTTTCAGGTGACATTGAATGTTCCTCAAAGTTTTCCTCTAGCAGATTTTTCTTACATTTAGTTTTAA"

# Tokenize the input sequence
inputs = tokenizer(dna_sequence, return_tensors='pt', padding=True, truncation=True, max_length=512)

# Forward pass through the model to get the logits
outputs = model(**inputs)
logits = outputs.logits

# Print the raw logits
print(f"Logits: {logits}")

# Apply softmax to get the probabilities
probabilities = F.softmax(logits, dim=1)

# Print the probabilities
print(f"Probabilities: {probabilities}")
