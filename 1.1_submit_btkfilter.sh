#!/bin/bash
## Launcher for the BlobToolKit (v2) combined script.
## Auto-counts assemblies and submits the array.
##
## Usage:  ./submit_btkfilter.sh <batch>     e.g. ./submit_btkfilter.sh 3rdbatch
##
## Run 00_setup_taxdump.sh ONCE before using this.

set -euo pipefail

BATCH="${1:?Usage: ./submit_btkfilter.sh <batch>   e.g. ./submit_btkfilter.sh 3rdbatch}"
PROCESS_EXEC=/dfs10/bio/abarrons/CA_curto/isolates/process_executables
ASSEMBLIES=/dfs10/bio/abarrons/CA_curto/isolates/${BATCH}/spades_out/assemblies

N=$(ls "$ASSEMBLIES"/*_scaffolds.fasta 2>/dev/null | wc -l)
if [ "$N" -eq 0 ]; then
  echo "ERROR: no *_scaffolds.fasta found in $ASSEMBLIES" >&2
  echo "       (check that '$BATCH' is spelled correctly and assemblies exist)" >&2
  exit 1
fi

echo "Batch:    $BATCH"
echo "Isolates: $N  ->  --array=1-$N"

sbatch --export=ALL,BATCH="$BATCH" \
       --array=1-"$N" \
       ${PROCESS_EXEC}/2.2_blobtools2_filter.sh
