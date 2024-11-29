import os
os.environ["CUDA_VISIBLE_DEVICES"] = "0"
import csv
import json
import logging
from dataclasses import dataclass, field
from typing import Union, Any, Dict, Sequence, Tuple, List, Optional
import torch
import transformers
import sklearn
import numpy as np
from torch.utils.data import Dataset, Subset
from sklearn.model_selection import KFold, StratifiedKFold
from transformers.models.bert.configuration_bert import BertConfig
from peft import (
    LoraConfig,
    get_peft_model,
    get_peft_model_state_dict,
)

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@dataclass
class ModelArguments:
    model_name_or_path: Optional[str] = field(default="zhihan1996/DNABERT-2-117M")
    use_lora: bool = field(default=False, metadata={"help": "Whether to use LoRA"})
    lora_r: int = field(default=8, metadata={"help": "Hidden dimension for LoRA"})
    lora_alpha: int = field(default=32, metadata={"help": "Alpha for LoRA"})
    lora_dropout: float = field(default=0.05, metadata={"help": "Dropout rate for LoRA"})
    lora_target_modules: str = field(default="query,value", metadata={"help": "Target modules for LoRA"})

@dataclass
class DataArguments:
    data_path: str = field(default=None, metadata={"help": "Path to the data directory."})
    kmer: int = field(default=-1, metadata={"help": "k-mer for input sequence. -1 means not using k-mer."})

@dataclass
class TrainingArguments(transformers.TrainingArguments):
    cache_dir: Optional[str] = field(default=None)
    run_name: str = field(default="run")
    optim: str = field(default="adamw_torch")
    model_max_length: int = field(default=512, metadata={"help": "Maximum sequence length."})
    gradient_accumulation_steps: int = field(default=1)
    per_device_train_batch_size: int = field(default=8)
    per_device_eval_batch_size: int = field(default=8)
    num_train_epochs: int = field(default=3)
    fp16: bool = field(default=False)
    logging_steps: int = field(default=100)
    save_steps: int = field(default=100)
    eval_steps: int = field(default=100)
    evaluation_strategy: str = field(default="steps")
    warmup_steps: int = field(default=50)
    weight_decay: float = field(default=0.01)
    learning_rate: float = field(default=1e-4)
    save_total_limit: int = field(default=3)
    load_best_model_at_end: bool = field(default=True)
    output_dir: str = field(default="output")
    find_unused_parameters: bool = field(default=False)
    checkpointing: bool = field(default=False)
    dataloader_pin_memory: bool = field(default=False)
    eval_and_save_results: bool = field(default=True)
    save_model: bool = field(default=False)
    seed: int = field(default=42)

def safe_save_model_for_hf_trainer(trainer: transformers.Trainer, output_dir: str):
    """Collects the state dict and dumps to disk."""
    state_dict = trainer.model.state_dict()
    if trainer.args.should_save:
        cpu_state_dict = {key: value.cpu() for key, value in state_dict.items()}
        del state_dict
        trainer._save(output_dir, state_dict=cpu_state_dict)

def generate_kmer_str(sequence: str, k: int) -> str:
    """Generate k-mer string from DNA sequence."""
    return " ".join([sequence[i:i+k] for i in range(len(sequence) - k + 1)])

def load_or_generate_kmer(data_path: str, texts: List[str], k: int) -> List[str]:
    """Load or generate k-mer string for each DNA sequence."""
    kmer_path = data_path.replace(".csv", f"_{k}mer.json")
    if os.path.exists(kmer_path):
        logger.info(f"Loading k-mer from {kmer_path}...")
        with open(kmer_path, "r") as f:
            kmer = json.load(f)
    else:
        logger.info(f"Generating k-mer...")
        kmer = [generate_kmer_str(text, k) for text in texts]
        with open(kmer_path, "w") as f:
            logger.info(f"Saving k-mer to {kmer_path}...")
            json.dump(kmer, f)
    return kmer

class SupervisedDataset(Dataset):
    """Dataset for supervised fine-tuning."""

    def __init__(self,
                 tokenizer: transformers.PreTrainedTokenizer,
                 data_path: str,
                 kmer: int = -1):

        super(SupervisedDataset, self).__init__()

        # Load data from the disk
        with open(data_path, "r") as f:
            reader = csv.reader(f)
            header = next(reader)  # Skip header
            data = list(reader)

        texts = [d[0] for d in data]
        labels = [int(d[1]) for d in data]

        if kmer != -1:
            logger.info(f"Using {kmer}-mer as input...")
            texts = load_or_generate_kmer(data_path, texts, kmer)

        output = tokenizer(
            texts,
            return_tensors="pt",
            padding="longest",
            max_length=tokenizer.model_max_length,
            truncation=True,
        )

        self.input_ids = output["input_ids"]
        self.attention_mask = output["attention_mask"]
        self.labels = torch.tensor(labels)
        self.num_labels = len(set(labels))

    def __len__(self):
        return len(self.input_ids)

    def __getitem__(self, idx) -> Dict[str, torch.Tensor]:
        return {
            'input_ids': self.input_ids[idx],
            'attention_mask': self.attention_mask[idx],
            'labels': self.labels[idx]
        }

@dataclass
class DataCollatorForSupervisedDataset(object):
    """Collate examples for supervised fine-tuning."""

    tokenizer: transformers.PreTrainedTokenizer

    def __call__(self, instances: Sequence[Dict]) -> Dict[str, torch.Tensor]:
        input_ids = torch.stack([instance['input_ids'] for instance in instances])
        attention_mask = torch.stack([instance['attention_mask'] for instance in instances])
        labels = torch.stack([instance['labels'] for instance in instances])
        return {
            'input_ids': input_ids,
            'attention_mask': attention_mask,
            'labels': labels
        }

def calculate_metrics(predictions: np.ndarray, labels: np.ndarray):
    return {
        "accuracy": sklearn.metrics.accuracy_score(labels, predictions),
        "f1": sklearn.metrics.f1_score(labels, predictions, average="macro", zero_division=0),
        "precision": sklearn.metrics.precision_score(labels, predictions, average="macro", zero_division=0),
        "recall": sklearn.metrics.recall_score(labels, predictions, average="macro", zero_division=0),
        "matthews_correlation": sklearn.metrics.matthews_corrcoef(labels, predictions),
    }

def compute_metrics(eval_pred):
    predictions, labels = eval_pred
    return calculate_metrics(predictions, labels)

def preprocess_logits_for_metrics(logits, labels):
    # If logits is a tuple, extract the first element
    if isinstance(logits, tuple):
        logits = logits[0]
    return torch.argmax(logits, dim=-1)

def train():
    parser = transformers.HfArgumentParser((ModelArguments, DataArguments, TrainingArguments))
    model_args, data_args, training_args = parser.parse_args_into_dataclasses()

    # Load tokenizer
    tokenizer = transformers.AutoTokenizer.from_pretrained(
        model_args.model_name_or_path,
        cache_dir=training_args.cache_dir,
        model_max_length=training_args.model_max_length,
        padding_side="right",
        use_fast=True,
        trust_remote_code=True
    )

    # Load the full dataset
    full_dataset = SupervisedDataset(
        tokenizer=tokenizer,
        data_path=os.path.join(data_args.data_path, "syntheticData_BRCATOTAL.csv"),  # Make sure this path is correct
        kmer=data_args.kmer
    )

    data_collator = DataCollatorForSupervisedDataset(tokenizer=tokenizer)

    # Set up K-Fold Cross Validation 
    #kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=training_args.seed)
    kf = KFold(n_splits=5, shuffle=True, random_state=training_args.seed)
    fold_results = []

    for fold, (train_indices, val_indices) in enumerate(kf.split(full_dataset)):
        logger.info(f"Starting fold {fold + 1}/{kf.get_n_splits()}")

        # Create data loaders for the current fold
        train_subset = Subset(full_dataset, train_indices)
        val_subset = Subset(full_dataset, val_indices)

        # Re-initialize the model for each fold
        config = BertConfig.from_pretrained(
            model_args.model_name_or_path,
            trust_remote_code=True,
            output_attentions=False,
            output_hidden_states=False,
            num_labels=full_dataset.num_labels,
            return_dict=True
        )

        model = transformers.AutoModelForSequenceClassification.from_pretrained(
            model_args.model_name_or_path,
            config=config,
            trust_remote_code=True
        )

        if model_args.use_lora:
            lora_config = LoraConfig(
                r=model_args.lora_r,
                lora_alpha=model_args.lora_alpha,
                target_modules=model_args.lora_target_modules.split(","),
                lora_dropout=model_args.lora_dropout,
                bias="none",
                task_type="SEQ_CLS",
                inference_mode=False,
            )
            model = get_peft_model(model, lora_config)
            model.print_trainable_parameters()

        # Define Trainer
        trainer = transformers.Trainer(
            model=model,
            tokenizer=tokenizer,
            args=training_args,
            data_collator=data_collator,
            train_dataset=train_subset,
            eval_dataset=val_subset,
            compute_metrics=compute_metrics,
            preprocess_logits_for_metrics=preprocess_logits_for_metrics,
        )

        # Train the model
        trainer.train()

        # Evaluate the model
        eval_results = trainer.evaluate(eval_dataset=val_subset)
        logger.info(f"Fold {fold + 1} evaluation results: {eval_results}")
        fold_results.append(eval_results)

        # Optionally save the model for the current fold
        if training_args.save_model:
            fold_output_dir = os.path.join(training_args.output_dir, f"fold_{fold + 1}")
            trainer.save_state()
            safe_save_model_for_hf_trainer(trainer=trainer, output_dir=fold_output_dir)

    # Compute average metrics across all folds
    avg_results = {}
    metric_keys = fold_results[0].keys()
    for key in metric_keys:
        avg_results[key] = np.mean([result[key] for result in fold_results])

    logger.info("Average cross-validation results:")
    for key, value in avg_results.items():
        logger.info(f"{key}: {value}")

    # Optionally save average results to a file
    if training_args.eval_and_save_results:
        results_path = os.path.join(training_args.output_dir, "cross_validation_results.json")
        os.makedirs(training_args.output_dir, exist_ok=True)
        with open(results_path, "w") as f:
            json.dump(avg_results, f, indent=4)

if __name__ == "__main__":
    train()
