import sys

import numpy as np
import pandas as pd
import anndata as ad
import celltypist

# Send all output and errors to the Snakemake log file
log_file = open(snakemake.log[0], "w")
sys.stdout = log_file
sys.stderr = log_file

h5ad            = snakemake.input["h5ad"]
model           = snakemake.params["model"]
majority_voting = snakemake.params["majority_voting"]
min_prob        = float(snakemake.params["min_prob"])
out_path        = snakemake.output["predictions"]

# The config can give majority_voting as a boolean (YAML true/false) or as a string
majority_voting = str(majority_voting).lower() in ("true", "1", "yes")

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
    All settings come from the Snakemake rule (input, params, output) defined at the top.
    - majority_voting: if the model supports it, the majority_voting column is used for the
      final labels instead of the default predicted_labels (more robust, but slower).
    - min_prob: confidence threshold. If the maximum predicted probability of a cell is below it,
      the cell is labeled "Unassigned" instead of a specific cell type.
    """
    print("PYTHON EXEC:", sys.executable)
    print("Loading h5ad:", h5ad)
    adata = ad.read_h5ad(h5ad)
    print("Cells:", adata.n_obs, "Genes:", adata.n_vars)
    print("majority_voting =", majority_voting, "| min_prob =", min_prob)

    # RUN CellTypist
    print("Running CellTypist with model:", model)
    pred = celltypist.annotate(adata, model=model, majority_voting=majority_voting)
    prob = pred.probability_matrix
    top1, margin = compute_confidence(prob)

    labels = pick_final_labels(pred, use_majority_voting=majority_voting).reindex(prob.index).astype(str)
    simple = labels.map(simplify_label)

    if min_prob > 0:
        low = top1 < min_prob # the probability we assigned in the config,
        #if the max probability for a cell is below this threshold, we consider it low confidence
        #.loc is used to select those cells with low confidence in a pandas series object
        labels.loc[low] = "Unassigned" # For low-confidence cells, we set the label to "Unassigned" instead of the predicted cell type
        simple.loc[low] = "Unassigned"

    out_df = pd.DataFrame({
        "cell": prob.index.astype(str),
        "celltypist_label": labels.values,
        "celltype_simple": simple.values,
        "max_prob": top1.values,
        "margin_top1_top2": margin.values,
    })
    out_df.to_csv(out_path, index=False)
    print("Saved:", out_path)
    print(out_df["celltype_simple"].value_counts().head(20).to_string())


if __name__ == "__main__":
    main()
    log_file.flush()