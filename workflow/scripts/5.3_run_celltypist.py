import argparse
import sys

import numpy as np
import pandas as pd
import anndata as ad
import celltypist

#this script runs CellTypist on a given h5ad file and saves the results to a CSV file.
# It also simplifies the cell type labels into broader categories.


def pick_final_labels(pred, use_majority_voting):
    """
    Return final 1D Series of labels for each cell.
    """
    labs = pred.predicted_labels
    if isinstance(labs, pd.DataFrame):
        if use_majority_voting and "majority_voting" in labs.columns:
            return labs["majority_voting"].astype(str)
        if "predicted_labels" in labs.columns:
            return labs["predicted_labels"].astype(str)
        return labs.iloc[:, 0].astype(str)
    return labs.astype(str)


def compute_confidence(prob):
    """
    prob : pred.probability_matrix , where there is the predicted probability for each cell (rows) and each cell type (columns).
    Compute confidence metrics from the probability matrix. Returns (max_prob, margin_top1_top2).
    -- max_prob: the highest predicted probability for each cell.
    -- margin_top1_top2: the difference between the top predicted probability and the second highest,
    """
    top1 = prob.max(axis=1)
    top2 = prob.apply(lambda r: np.partition(r.values, -2)[-2], axis=1)
    return top1, top1 - top2


def simplify_label(label):
    """
    Map detailed CellTypist labels to broader categories for easier interpretation with my previous annotations. 
    """
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
    """
    Main function to run CellTypist on an h5ad file and save the results.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--majority-voting", required=True)
    ## If the model supports majority voting,
    #this will use the majority_voting column for final labels instead of the default predicted_labels. 
    # This can provide more robust predictions by considering multiple models or iterations, but it requires that the model was trained with majority voting enabled.
    # If you want to use majority voting : False  is the default , 
    parser.add_argument("--min-prob", required=True, type=float)
     # This allows the user to set a confidence threshold for the predictions. 
    # If the maximum predicted probability for a cell is below this threshold, the cell will be labeled as "Unassigned"
    # instead of being assigned to a specific cell type. This can help improve the reliability of the annotations by filtering out low-confidence predictions.
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    majority_voting = args.majority_voting.lower() in ("true", "1", "yes")

    print("PYTHON EXEC:", sys.executable)
    print("Loading h5ad:", args.h5ad)
    adata = ad.read_h5ad(args.h5ad)
    print("Cells:", adata.n_obs, "Genes:", adata.n_vars)
    print("majority_voting =", majority_voting, "| min_prob =", args.min_prob)

    # RUN CellTypist
    print("Running CellTypist with model:", args.model)
    pred = celltypist.annotate(adata, model=args.model, majority_voting=majority_voting)
    prob = pred.probability_matrix
    top1, margin = compute_confidence(prob)

    labels = pick_final_labels(pred, use_majority_voting=majority_voting).reindex(prob.index).astype(str)
    simple = labels.map(simplify_label)

    if args.min_prob > 0:
        low = top1 < args.min_prob # the probability we assigned in the argumnets, 
        #if the max probability for a cell is below this threshold, we consider it low confidence
        #.loc is used to selecct those cells with low confidence in a pandas series object
        labels.loc[low] = "Unassigned" # For low-confidence cells, we set the label to "Unassigned" instead of the predicted cell type
        simple.loc[low] = "Unassigned" # # For low-confidence cells, we set the label to "Unassigned"

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