#!/usr/bin/env bash
# Snakemake runs this once, inside the activated env, right after creating supercell.yaml.
# SuperCell is not on conda, so it is installed from GitHub, pinned to a commit (v1.0)
# for reproducibility. Installing it here (and not inside the scripts) avoids parallel
# jobs racing to install it at the same time.
set -euo pipefail
Rscript -e 'remotes::install_github("GfellerLab/SuperCell@f0a041e92a611dba650f0ab1f37d7075ff25f6e9", upgrade = "never")'
Rscript -e 'library(SuperCell); message("SuperCell ", packageVersion("SuperCell"), " installed")'
