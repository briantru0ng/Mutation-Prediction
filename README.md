(this readme.md was created out of a lingering need to always have a read me. but i'm busy and i just used chatgpt to make the read me for me since i have other things due on Cinco de Mayo)
# Mutation-Prediction

## Overview

This project evaluates whether custom transformer-based models, fine-tuned from pretrained protein language models such as **ESM2** and **ProtBERT**, can outperform classical alignment-based tools like **SIFT** in predicting the pathogenicity of single amino acid variants (SAVs). We curated a balanced dataset from **ClinVar**, focusing on key cancer-related genes (e.g., TP53, BRCA1/2, KRAS, PIK3CA and many more), and built a comprehensive pipeline involving both frozen and fine-tuned transformer models.


## Models Compared

| Method            | Feature Type         | Trainable? | Comments                          |
|------------------|----------------------|------------|-----------------------------------|
| **SIFT**         | Hand-engineered      | ❌         | High pathogenic recall, low precision |
| **Pretrained EMS2** | Transformer embeddings | ❌         | Mean-pooled, logistic regression   |
| **Pretrained ProtBERT** | Transformer embeddings | ❌         | 768-dim, pooled last hidden state |
| **Custom EMS2**   | Fine-tuned Transformer | ✅         | Additional encoder layers, partial unfreezing |

---
## Results

| Model              | Accuracy | F1 (Patho/Benign) | Notes                          |
|-------------------|----------|-------------------|--------------------------------|
| **SIFT**          | ~70%     | 0.77 / 0.63       | High recall for pathogenic SAVs |
| **Pretrained EMS2** | ~65%   | 0.65 / 0.66       | Balanced but shallow performance |
| **Pretrained ProtBERT** | ~66% | 0.66 / 0.66    | Similar to EMS2                |
| **Custom EMS2**    | **~82%** | **0.82 / 0.81**   | Best overall performance       |

---

##  How to Run

### 1. Set up your environment

```bash
# Create a new virtual environment
python -m venv .venv
source .venv/bin/activate    # On Windows: .venv\Scripts\activate

# Install required packages
pip install -r requirements.txt 
```

### 2. Run experiments

Depending on your experiment of interest:

- **Baseline SIFT Analysis**:
  - Open and run `sift.ipynb`
  - View results and plots in `sift_result.ipynb`

- **Pretrained Transformer Models (no fine-tuning)**:
  - Run `pretrained-ems.ipynb` for ESM2-based embedding and classification
  - Run `pretrained-protbert.ipynb` for ProtBERT-based embedding and classification

- **Custom Transformer Model (fine-tuned)**:
  - Run `custom-ems.ipynb` to fine-tune ESM2 with added Transformer encoder layers
