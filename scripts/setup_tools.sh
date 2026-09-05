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
if [ -d "${TOOLS_DIR}/miReact/.git" ]; then
  echo "miReact already cloned at ${TOOLS_DIR}/miReact -- checking it's valid rather than re-cloning..."
  if [ -f "${TOOLS_DIR}/miReact/data/tarbase.rds" ]; then
    echo "Looks like a complete, valid clone (tarbase.rds present). Skipping clone step."
  else
    echo "ERROR: ${TOOLS_DIR}/miReact exists but looks INCOMPLETE (no data/tarbase.rds)."
    echo "This is likely a partial clone from an interrupted earlier run. Remove it and re-run:"
    echo "  rm -rf ${TOOLS_DIR}/miReact"
    exit 1
  fi
else
  git clone https://github.com/muhligs/miReact.git "${TOOLS_DIR}/miReact"
fi

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

echo "=== Installing bayesReact into the SAME shared R library (not a separate location) ==="
Rscript -e '
  .libPaths(c("'"${TOOLS_DIR}"'/R_library", .libPaths()))
  if (requireNamespace("bayesReact", quietly = TRUE, lib.loc = "'"${TOOLS_DIR}"'/R_library")) {
    cat("bayesReact already installed -- skipping (re-run devtools::install_github manually if you need to update it)\n")
  } else {
    devtools::install_github("JakobSkouPedersenLab/bayesReact", dependencies = TRUE, lib = "'"${TOOLS_DIR}"'/R_library")
  }
'

echo "=== Done. Add this to your .Renviron for the library to load automatically: ==="
echo "R_LIBS_USER=${TOOLS_DIR}/R_library"
