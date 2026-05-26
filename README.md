# Signature Scoring for Drug–Disease Transcriptomic Reversal Analysis

## Overview

This package provides a Python script (`signature_scoring.py`) for calculating similarity or reversal scores between drug-induced and disease-associated transcriptomic signatures.

It supports four signature-matching methods:

- **KS**: Kolmogorov–Smirnov-based score
- **XSUM**: extreme-sum score
- **XCOS**: cosine similarity score
- **SCS**: signature concordance score

The script is designed for transcriptome-based drug repositioning analysis, where:

- rows represent genes
- columns represent signatures
- values represent expression changes (e.g., log2 fold change)

---

## Requirements

- Python 3.9 or later is recommended

Install dependencies with:

```bash
pip install -r requirements.txt
```

Required packages:

- numpy
- pandas
- tqdm

---

## Input files

### 1. Drug expression matrix

A tab-separated file in which:

- rows = genes
- columns = drug signatures
- values = log2FC or other signed differential expression values

Example:

```tsv
gene  Drug_A  Drug_B
Gene1 1.25  -0.68
Gene2 -0.44 0.92
Gene3 0.77  -1.10
```

### 2. Disease expression matrix

A tab-separated file in which:

- rows = genes
- columns = disease signatures
- values = log2FC or other signed differential expression values

Example:

```tsv
gene  Disease_1 Disease_2
Gene1 -1.05 0.62
Gene2 0.83  -0.74
Gene3 -0.56 1.18
```

### 3. Optional pair file

If provided, only specified drug–disease pairs will be scored.

Required columns:

- `drug`
- `disease`

Example:

```tsv
drug  disease
Drug_A  Disease_1
Drug_B  Disease_7
Drug_C  Disease_7
```

---

## Output

The script writes one output file per disease signature into the specified output directory.

Each output file is a tab-separated table containing:

- `disease`
- `drug`
- `method`
- `score`

Example:

```tsv
disease drug  method  score
Disease_1 Drug_A  KS  -0.42
Disease_1 Drug_A  XSUM  -13.57
Disease_1 Drug_A  XCOS  -0.61
Disease_1 Drug_A  SCS 84
```

Output filenames are automatically converted to Windows-safe names.

---

## Usage

### 1. Full scoring mode

This mode calculates scores for all drug × disease combinations.

```bash
python signature_scoring.py drug_log2fc.tsv disease_log2fc.tsv score_output/
```

### 2. Pair mode

This mode calculates scores only for specified drug–disease pairs.

```bash
python signature_scoring.py drug_log2fc.tsv disease_log2fc.tsv score_output/ pair_file.tsv
```

### 3. Import as a Python module

```python
from signature_scoring import run_scoring

run_scoring(
    drug_expr_path="drug_log2fc.tsv",
    disease_expr_path="disease_log2fc.tsv",
    output_dir="score_output",
    pair_file="pair_file.tsv"  # set to None if not needed
)
```

---

## Implemented methods

### KS
A Kolmogorov–Smirnov-like enrichment score based on the positions of disease up/down genes in the ranked drug signature.

### XSUM
Calculates the sum of drug expression values over disease up-regulated genes minus the sum over disease down-regulated genes.

### XCOS
Computes cosine similarity between the selected extreme genes from the drug signature and the disease signature.

### SCS
Calculates a concordance-style score based on overlap between up/down gene sets of drug and disease signatures.

---

## Default parameters

- Default methods:
  - `KS`
  - `XSUM`
  - `XCOS`
  - `SCS`
- Default `topN` for each method: `200`

---

## Notes

1. Drug and disease matrices are aligned using the intersection of gene identifiers.
2. Missing values are removed within each signature before ranking.
3. Non-numeric values in the input matrices are automatically converted to `NaN`.
4. If no overlapping genes exist between the drug and disease matrices, the script will raise an error.
5. If an output file for a disease already exists, that disease will be skipped.

---

## Suggested package contents

```text
signature_scoring.py
requirements.txt
README.md
```

Optional:

```text
example_input/
example_output/
pair_file.tsv
```

---

## Data availability

The datasets used in this project are available on Figshare:

[Figshare link]


## Author

ZHANG SULIN

## Date

2026-01-26
