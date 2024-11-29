import pandas as pd
from sklearn.model_selection import train_test_split

# Load your dataset from a CSV file
data = pd.read_csv('/NFSHOME/lmasci/DNABERT_2/DATA/syntheticData_BRCATOTAL_500bp.csv')

# Define the size of each set (example: 70% train, 15% dev, 15% test)
train_size = 0.7
dev_size = 0.15
test_size = 0.15

# First, split the data into train and temp (test + dev)
train_data, temp_data = train_test_split(data, test_size=(1 - train_size), random_state=42)

# Then split the temp data into dev and test
dev_data, test_data = train_test_split(temp_data, test_size=(test_size / (dev_size + test_size)), random_state=42)

# Save the split datasets to new CSV files
train_data.to_csv('/NFSHOME/lmasci/DNABERT_2/sample_BRCATOT_500bp/train.csv', index=False)
dev_data.to_csv('/NFSHOME/lmasci/DNABERT_2/sample_BRCATOT_500bp/dev.csv', index=False)
test_data.to_csv('/NFSHOME/lmasci/DNABERT_2/sample_BRCATOT_500bp/test.csv', index=False)

print(f"Train set: {len(train_data)} samples")
print(f"Dev set: {len(dev_data)} samples")
print(f"Test set: {len(test_data)} samples")
