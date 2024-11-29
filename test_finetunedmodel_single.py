#import torch
from transformers import AutoModelForSequenceClassification, AutoTokenizer
#from transformers import AutoModel, AutoTokenizer

from torch.nn import functional as F

# Load the fine-tuned model and tokenizer
tokenizer = AutoTokenizer.from_pretrained("/NFSHOME/lmasci/DNABERT_2/finetune/output/BRCATOT_TEST", trust_remote_code=True)
model = AutoModelForSequenceClassification.from_pretrained("/NFSHOME/lmasci/DNABERT_2/finetune/output/BRCATOT_TEST", trust_remote_code=True)
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

# Get the predicted class by taking the index of the highest probability
#predicted_class_idx = torch.argmax(probabilities, dim=1).item()
#
# Assuming binary classification (0: non-functional, 1: functional), define the class labels
#class_labels = {0: "non-functional", 1: "functional"}
#
## Print the predicted class and its probability
#predicted_class = class_labels[predicted_class_idx]
#predicted_class_prob = probabilities[0, predicted_class_idx].item()
#
#print(f"Predicted class: {predicted_class} (Probability: {predicted_class_prob:.4f})")