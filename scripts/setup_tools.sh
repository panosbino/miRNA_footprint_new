#!/bin/bash
# Run this ONCE on Dardel (login node -- needs internet access) before
# running mireact_comparison.R. Sets up the shared tools/ infrastructure
# under the NEW miRNA_footprint_new structure.

set -e  # stop on first error, don't silently continue past a failed step

BASE_DIR="/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR="${BASE_DIR}/tools"

mkdir -p "${TOOLS_DIR}"
cd "${TOOLS_DIR}"

echo "=== Cloning miReact ==="
git clone https://github.com/muhligs/miReact.git

echo "=== Setting up shared R library ==="
mkdir -p "${TOOLS_DIR}/R_library"

echo "=== Verifying tarbase.rds came along with the clone ==="
ls -la "${TOOLS_DIR}/miReact/data/tarbase.rds"

echo "=== Installing Regmex into the shared R library ==="
Rscript -e '
  .libPaths(c("'"${TOOLS_DIR}"'/R_library", .libPaths()))
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", lib = "'"${TOOLS_DIR}"'/R_library")
  devtools::install_github("muhligs/Regmex", dependencies = FALSE, lib = "'"${TOOLS_DIR}"'/R_library")
'

echo "=== Done. Add this to your .Renviron for the library to load automatically: ==="
echo "R_LIBS_USER=${TOOLS_DIR}/R_library"
