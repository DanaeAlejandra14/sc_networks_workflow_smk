# CLAUDE.md

Snakemake (>= 9) pipeline that builds **gene co-expression networks from snRNA-seq of the
prefrontal cortex (ROSMAP, Alzheimer's disease)**, one per phenotype x cell type, with
ARACNe-AP on SuperCell metacells, and then analyzes those networks.

Talk to the user in Spanish. Code, code comments and commit messages are in English;
comments in `config/config.yaml` are in Spanish.

## Pipeline

1. Demultiplex 10x batches by individual, attach phenotype (`1*`)
2. QC, doublets (scDblFinder), LogNormalize (`2*`-`4*`)
3. Cell type annotation: reference atlas + CellTypist for unassigned cells (`5*`)
4. Split by phenotype and cell type, merge, Harmony integration (`6*`-`8*`)
5. Subsample every phenotype to the smallest one (`8.2`; Astrocyte: 31,330 cells each), SuperCell metacells per gamma (`9*`)
6. Gene set for ARACNe (`9.5_select_aracne_genes`): one set per cell type, shared by every phenotype and gamma. On single cells of all phenotypes pooled: detected in >= 5% of cells (`aracne.genes.min_frac_cells`), MT-/ribosomal removed, then the 4000 most variable (vst, `n_variable`). Output `results/aracne_genes/{cell_type}/genes.txt`
7. ARACNe-AP: matrix export, MI threshold, 100 bootstraps, consolidation (`10*`, shell rules). Every selected gene is a regulator (`centroids.txt`)
8. Network analysis (`11*`-`19*`, see below)

Phenotypes: `control`, `MCI`, `AD-NC_ASYM` (AD pathology, no symptoms), `AD-NC_SYM`.
What gets built is set by `aracne.targets` and `supercell.gamma_selected` in the config (now: 4 phenotypes x Astrocyte x gamma 60/70/80 = 12 networks, ~522/448/392 metacells).

## Network analysis (scripts 11-19)

| # | Rule | Output |
|---|---|---|
| 11 | `network_metrics` | `node_metrics.tsv`: PageRank, degree, strength, betweenness of every gene |
| 12 | `plot_networks` | network figures (top 100 genes by MI), per phenotype + panel |
| 13 | `plot_top_hubs` | top 10 PageRank hubs per phenotype |
| 14 | `ego_network` | ego network tables of each gene in `ego_networks.genes` |
| 15 | `ego_enrichment` | GO enrichment of the ego network (env `enrichment.yaml`) |
| 16 | `plot_ego_network` | ego network panel, neighbors dot matrix, bootstrap support |
| 17 | `plot_ego_enrichment` | GO dot plot |
| 18 | `network_proximity` | proximity of the gene to `proximity.gene_sets`, degree-matched null, final + bootstrap networks |
| 19 | `plot_proximity` | proximity z-scores per phenotype + Wilcoxon vs `proximity.reference` |

Layout:
- `results/aracne/{phenotype}/{cell_type}/gamma{gamma}/`: ARACNe output (`network.txt`, `bootstrap_networks/`, `centroids.txt`, `scmatrix.tsv`)
- `results/net_analysis/{phenotype}/{cell_type}/gamma{gamma}/`: tables, one job per network
- `results/net_analysis/comparisons/{cell_type}/gamma{gamma}/`: tables comparing phenotypes
- `results/plots/{net_vis,hubs,ego,proximity}/{cell_type}/gamma{gamma}/`: figures. `net_vis` is only for network figures

## Conventions

- **One rule per script**, same number and name. Calculations run per network (`{phenotype}` wildcard); figures read the per-phenotype tables and compare phenotypes (one job per cell type and gamma).
- **Parameters live in the config**, passed through `params`; genes of interest too (`ego_networks.genes`, `hubs.highlight_genes`, `proximity.gene_sets`). Adding a gene to `ego_networks.genes` creates its ego and proximity analyses with no code changes.
- **Shared R code** in `workflow/scripts/R/`, loaded with `source(snakemake@input[[...]])` and declared as a rule input so edits re-run the rules:
  - `theme_pub.R`: phenotype order (`phenotype_levels`) and colors, cell type colors, `pal_highlight`, `theme_pub()`, `theme_pub_void()`, color scales, `save_pub()` (PNG + PDF). Every plot uses it.
  - `network_utils.R`: `read_aracne_network()` (final or bootstrap files; collapses A->B / B->A into one edge, max MI) and `aracne_graph()` (undirected igraph, MI as weight).
- **Phenotypes in clinical order** in every figure: control -> MCI -> AD-NC_ASYM -> AD-NC_SYM.
- **Scripts write to the log** with `sink()` like the existing ones.
- **Plot style the user prefers:** no captions; subtitles only when asked; large gene labels in boxes; one shared legend in multi-panel figures (patchwork `guides = "collect"` with shared scale limits).
- **Conda envs:** `network_analysis.yaml` for 11-14, 16, 18, 19; `enrichment.yaml` (clusterProfiler, slow to install) for 15 and 17; `supercell.yaml` for 8.3, 9*, 10. Changing an env file re-runs every rule that uses it.
- **SuperCell is not on conda:** `workflow/envs/supercell.post-deploy.sh` installs it from GitHub (pinned commit) when Snakemake creates the env. Never install packages inside scripts: parallel jobs race on the install (00LOCK) and fail.

## Pitfalls

- **igraph `distances()` uses the `weight` attribute by default.** ARACNe graphs carry MI as weight, so always pass `weights = NA` to count edges.
- **ggplot2 >= 3.5 only draws legend keys for values present in a layer.** For legends that must match across patchwork panels, use `show.legend = c(<aes> = TRUE, <other aes> = FALSE)` and set scale `limits`.
- **Final ARACNe networks are very sparse and fragmented** (~4k genes, ~1.5 edges per gene; a gene can be missing, e.g. SLC1A2 in control). Bootstrap networks are dense and connected (~15k genes, ~450k edges). Distance-based metrics are only meaningful on the bootstraps.
- **Bootstrap networks are large:** reading one takes ~5 s; 1000 BFS on one take ~3 min. Keep random sampling small (proximity uses 100 sources x 10 sets).
- **`10_export_aracne_matrices.R` changed (a comment) after the matrices were made**, so Snakemake's code trigger wants to rebuild `scmatrix.tsv` / `centroids.txt` and all of ARACNe when a rule asks for `centroids.txt` (rule 15). If the matrices are fine: `snakemake --cleanup-metadata results/aracne/*/Astrocyte/gamma10/{scmatrix.tsv,centroids.txt}`.
- **ARACNe cost grows with the square of the number of genes** (all genes are regulators); gamma barely matters. Measured with ~11-13k genes and 627-1566 metacells: ~10.5 h and ~21 GB per bootstrap. That is why the gene set is capped (`aracne.genes.n_variable`).
- **ARACNe hangs after a `java.lang.OutOfMemoryError`:** the worker threads die but the JVM keeps running (1 core at 100%, no output), so Snakemake waits forever. Check `grep OutOfMemory logs/aracne/*.log` when bootstraps take too long; if so, kill the java and raise `aracne.memory.bootstrap_gb`. Java gets `-Xmx` = memory - 2 GB.
- **ARACNe does not respect `--threads`:** a bootstrap with 4 threads uses ~5 cores (up to ~20 at times), so the real load is higher than `--cores`.
- **A gene filter by fraction of metacells depends on gamma:** bigger metacells sum more cells, so more genes pass. That is why genes are selected on single cells, once per cell type.
- **ARACNe runs on the server:** `--cores 36 --resources mem_gb=300` (9 bootstraps at a time with 4 threads and 12 GB each).
- **This machine has only `results/aracne/`** (no `data/`, `resources/` or intermediate results); the upstream pipeline runs on the server. Always dry-run (`-n`) first and check the job list does not include upstream rules.

## Known caveats of the current networks

The networks in `results/aracne/*/Astrocyte/gamma10` are preliminary: 10 bootstraps, different metacell counts per phenotype (106-250), MT and ribosomal genes not filtered (they appear as hubs). The new run (gamma 60/70/80, 100 bootstraps, same number of cells per phenotype, shared 4000-gene set without MT/ribosomal) fixes those. Still open: phenotype confounded with sex (AD-NC_ASYM mostly male, MCI mostly female) and sex-chromosome genes (XIST, chrY) are kept. Interpret differences between phenotypes as hypotheses.

## Running

```bash
# Snakemake env: workflow/envs/snakemake.yaml (here: /home/tiamat/miniforge3/envs/snakemake)
snakemake -n --use-conda <targets>          # always check the DAG first
snakemake -c4 --use-conda <targets>
```

`results/`, `logs/`, `.snakemake/`, `data/` and `resources/` are git-ignored. The repo is shared
(origin: DanaeAlejandra14/sc_networks_workflow_smk): work on a branch and open a PR to `main`.
