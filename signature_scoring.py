"""
@filename:signature_scoring.py
@author:Wintie
@time:2026-01-26
"""

import os
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
import re


# ============================================================
# I/O & Preprocess
# ============================================================

def safe_filename(name):
    """
    Make Windows-safe filename from signature name.
    """
    return re.sub(r'[\\/:*?"<>|]', '_', str(name))


def read_expression_matrix(path):
    """
    Read expression matrix.
    Assumption:
      - rows: genes
      - columns: signatures
      - values: log2FC
    """
    df = pd.read_csv(path, sep="\t", engine="python", index_col=0)
    # 强制转数值，避免 object 导致 XSUM / XCOS 失效
    df = df.apply(pd.to_numeric, errors="coerce")
    if df.empty:
        raise ValueError(f"Empty expression matrix: {path}")
    return df


def read_pair_file(pair_file):
    """
    Read pair file with columns: drug, disease
    """
    pair_df = pd.read_csv(pair_file, sep=None, engine="python")
    if not {"drug", "disease"}.issubset(pair_df.columns):
        raise ValueError("pair_file must contain columns: drug, disease")
    return pair_df


def align_matrices(drug_expr, disease_expr):
    """
    Align gene universe between drug and disease matrices.
    """
    common_genes = drug_expr.index.intersection(disease_expr.index)
    if len(common_genes) == 0:
        raise ValueError("No overlapping genes between drug and disease matrices.")

    drug_expr = drug_expr.loc[common_genes]
    disease_expr = disease_expr.loc[common_genes]
    return drug_expr, disease_expr


def preprocess_signature(vec):
    """
    Only do ranking once.
    """
    vec = vec.dropna()
    if vec.empty:
        return None

    ranked = vec.sort_values(ascending=False)
    return {
        "ranked": ranked,
        "vector": ranked
    }


def get_up_down(sig, topN):
    ranked = sig["ranked"]
    if topN * 2 > len(ranked):
        return set(), set()

    up = set(ranked.iloc[:topN].index)
    down = set(ranked.iloc[-topN:].index)
    return up, down


# ============================================================
# Scoring Methods (LOGIC UNCHANGED)
# ============================================================

def ks_enrichment(ranked_genes, query_genes):
    if len(query_genes) == 0:
        return 0.0

    gene_to_rank = {g: i + 1 for i, g in enumerate(ranked_genes)}
    ranks = sorted([gene_to_rank[g] for g in query_genes if g in gene_to_rank])

    if len(ranks) == 0:
        return 0.0

    n = len(ranked_genes)
    m = len(ranks)

    d = np.arange(1, m + 1) / m - np.array(ranks) / n
    a = d.max()
    b = -d.min() + 1 / m

    return a if a > b else -b


def ks_score(drug_sig, disease_sig, topN):
    ranked = drug_sig["ranked"].index.tolist()
    d_up, d_down = get_up_down(disease_sig, topN)

    score_up = ks_enrichment(ranked, d_up)
    score_down = ks_enrichment(ranked, d_down)

    if score_up * score_down <= 0:
        return score_up - score_down
    else:
        return 0.0


def xsum_score(drug_sig, disease_sig, topN):
    drug_vec = drug_sig["vector"]
    d_up, d_down = get_up_down(disease_sig, topN)

    score_up = drug_vec.loc[drug_vec.index.intersection(d_up)].sum()
    score_down = drug_vec.loc[drug_vec.index.intersection(d_down)].sum()

    return score_up - score_down


def xcos_score(drug_sig, disease_sig, topN):
    drug_ranked = drug_sig["ranked"]
    disease_vec = disease_sig["vector"]

    if topN * 2 > len(drug_ranked):
        return np.nan

    selected_genes = drug_ranked.iloc[:topN].index.union(
        drug_ranked.iloc[-topN:].index
    )

    common = selected_genes.intersection(disease_vec.index)
    if len(common) == 0:
        return np.nan

    v1 = drug_ranked.loc[common].values
    v2 = disease_vec.loc[common].values

    denom = np.linalg.norm(v1) * np.linalg.norm(v2)
    if denom == 0:
        return np.nan

    return float(np.dot(v1, v2) / denom)


def scs_score(drug_sig, disease_sig, topN):
    d_up, d_down = get_up_down(disease_sig, topN)
    drug_up, drug_down = get_up_down(drug_sig, topN)

    concordant = len(d_up & drug_down) + len(d_down & drug_up)
    discordant = len(d_up & drug_up) + len(d_down & drug_down)

    return concordant - discordant


# ============================================================
# Main Runner (Disease-level checkpoint + pair mode)
# ============================================================

def run_scoring(
    drug_expr_path,
    disease_expr_path,
    output_dir,
    methods=("KS", "XSUM", "XCOS", "SCS"),
    topN_params=None,
    pair_file=None,
    verbose=True
):
    os.makedirs(output_dir, exist_ok=True)

    drug_expr = read_expression_matrix(drug_expr_path)
    disease_expr = read_expression_matrix(disease_expr_path)

    drug_expr, disease_expr = align_matrices(drug_expr, disease_expr)

    if topN_params is None:
        topN_params = {m: 200 for m in methods}

    # preprocess drug signatures
    drug_sigs = {}
    for drug in drug_expr.columns:
        sig = preprocess_signature(drug_expr[drug])
        if sig is not None:
            drug_sigs[drug] = sig

    # pair mode
    if pair_file is not None:
        pair_df = read_pair_file(pair_file)
        disease_list = pair_df["disease"].unique()
    else:
        disease_list = disease_expr.columns

    if verbose:
        disease_list = tqdm(disease_list, desc="Processing diseases")

    for disease in disease_list:
        if disease not in disease_expr.columns:
            continue

        safe_name = safe_filename(disease)
        out_file = os.path.join(output_dir, f"{safe_name}.tsv")
        if os.path.exists(out_file):
            continue

        disease_sig = preprocess_signature(disease_expr[disease])
        if disease_sig is None:
            continue

        records = []

        # select drugs
        if pair_file is not None:
            drugs_to_run = pair_df.loc[
                pair_df["disease"] == disease, "drug"
            ].unique()
        else:
            drugs_to_run = drug_sigs.keys()

        for drug in drugs_to_run:
            if drug not in drug_sigs:
                continue
            drug_sig = drug_sigs[drug]

            for method in methods:
                topN = topN_params.get(method, 200)

                try:
                    if method == "KS":
                        score = ks_score(drug_sig, disease_sig, topN)
                    elif method == "XSUM":
                        score = xsum_score(drug_sig, disease_sig, topN)
                    elif method == "XCOS":
                        score = xcos_score(drug_sig, disease_sig, topN)
                    elif method == "SCS":
                        score = scs_score(drug_sig, disease_sig, topN)
                    else:
                        continue
                except Exception:
                    score = np.nan

                records.append({
                    "disease": disease,
                    "drug": drug,
                    "method": method,
                    "score": score
                })

        if len(records) == 0:
            continue

        df_out = pd.DataFrame.from_records(records)
        df_out.to_csv(out_file, sep="\t", index=False)


# ============================================================
# CLI
# ============================================================

if __name__ == "__main__":
    if len(sys.argv) not in (4, 5):
        print(
            "Usage:\n"
            "  python signature_scoring.py <drug_expr.tsv> <disease_expr.tsv> <output_dir>\n"
            "  python signature_scoring.py <drug_expr.tsv> <disease_expr.tsv> <output_dir> <pair_file.tsv>\n"
        )
        sys.exit(1)

    run_scoring(
        drug_expr_path=sys.argv[1],
        disease_expr_path=sys.argv[2],
        output_dir=sys.argv[3],
        pair_file=sys.argv[4] if len(sys.argv) == 5 else None,
        verbose=True
    )


"""
========================
USAGE EXAMPLES
========================

1) 全量计算（默认，drug × disease 全组合）

python signature_scoring.py \
    drug_log2fc.tsv \
    disease_log2fc.tsv \
    score_output/


2) 指定组合计算（pair 模式）

pair_file.tsv 格式：
--------------------
drug        disease
Drug_A      Disease_1
Drug_B      Disease_7
Drug_C      Disease_7
--------------------

运行：
python signature_scoring.py \
    drug_log2fc.tsv \
    disease_log2fc.tsv \
    score_output/ \
    pair_file.tsv


3) 在 Python 中作为模块调用

from signature_scoring import run_scoring

run_scoring(
    drug_expr_path="drug_log2fc.tsv",
    disease_expr_path="disease_log2fc.tsv",
    output_dir="score_output",
    pair_file="pair_file.tsv",   # 不需要则设为 None
)
"""
