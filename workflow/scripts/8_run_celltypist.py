import argparse
import sys

import numpy as np
import pandas as pd
import anndata as ad
import celltypist


def pick_final_labels(pred, use_majority_voting):
    labs = pred.predicted_labels
    if isinstance(labs, pd.DataFrame):
        if use_majority_voting and "majority_voting" in labs.columns:
            return labs["majority_voting"].astype(str)
        if "predicted_labels" in labs.columns:
            return labs["predicted_labels"].astype(str)
        return labs.iloc[:, 0].astype(str)
    return labs.astype(str)


def compute_confidence(prob):
    top1 = prob.max(axis=1)
    top2 = prob.apply(lambda r: np.partition(r.values, -2)[-2], axis=1)
    return top1, top1 - top2


def simplify_label(label):
    s = str(label)
    if s.startswith(("L2", "L3", "L4", "L5", "L6")):
        return "Excitatory Neurons"
    if s.startswith("InN") or s.startswith("IN"):
        return "Inhibitory Neurons"
    if s.startswith("Astro"):
        return "Astrocyte"
    if s.startswith("Micro"):
        return "Microglia"
    if s.startswith("Oligo"):
        return "Oligodendrocytes"
    if s.startswith("OPC"):
        return "OPCs"
    if s.startswith("Endo"):
        return "Endothelial"
    if s.startswith(("Macro", "Myeloid", "DC", "B ", "T ", "NK")):
        return "Immune"
    return "Other/Unknown"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--majority-voting", required=True)
    parser.add_argument("--min-prob", required=True, type=float)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    majority_voting = args.majority_voting.lower() in ("true", "1", "yes")

    print("PYTHON EXEC:", sys.executable)
    print("Loading h5ad:", args.h5ad)
    adata = ad.read_h5ad(args.h5ad)
    print("Cells:", adata.n_obs, "Genes:", adata.n_vars)
    print("majority_voting =", majority_voting, "| min_prob =", args.min_prob)

    pred = celltypist.annotate(adata, model=args.model, majority_voting=majority_voting)
    prob = pred.probability_matrix
    top1, margin = compute_confidence(prob)

    labels = pick_final_labels(pred, use_majority_voting=majority_voting).reindex(prob.index).astype(str)
    simple = labels.map(simplify_label)

    if args.min_prob > 0:
        low = top1 < args.min_prob
        labels.loc[low] = "Unassigned"
        simple.loc[low] = "Unassigned"

    out = pd.DataFrame({
        "cell": prob.index.astype(str),
        "celltypist_label": labels.values,
        "celltype_simple": simple.values,
        "max_prob": top1.values,
        "margin_top1_top2": margin.values,
    })
    out.to_csv(args.out, index=False)
    print("Saved:", args.out)
    print(out["celltype_simple"].value_counts().head(20).to_string())


if __name__ == "__main__":
    main()