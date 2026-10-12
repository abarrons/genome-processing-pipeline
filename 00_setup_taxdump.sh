#!/bin/bash
## ONE-TIME SETUP: download and unpack the NCBI taxdump for BlobTools2.
## Run this ONCE (not per batch, not in an array) before the blob+filter step.
## Not a SLURM array job -- just run on a login node or as a simple sbatch.
##
## Usage:  bash 00_setup_taxdump.sh
##
## The resulting directory path is what you pass as TAXDUMP in the pipeline.

set -euo pipefail

## Where to put it (shared location, used by all batches)
TAXDUMP=/dfs10/bio/abarrons/databases/taxdump

mkdir -p "$TAXDUMP"
cd "$TAXDUMP"

echo "Downloading NCBI new_taxdump into $TAXDUMP ..."
## new_taxdump includes the taxid lineage files BlobTools2 needs
curl -L -o new_taxdump.tar.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz

echo "Unpacking ..."
tar -xzf new_taxdump.tar.gz

echo "Done. TAXDUMP=$TAXDUMP"
echo "Contents:"
ls -1 "$TAXDUMP" | head
