#!/bin/bash
## Launcher: counts R1 files in the batch and submits the SPAdes array
## with the correct --array range. No need to edit array size by hand.
##
## Usage: ./1.0_submit_batch.sh <batch>   e.g. ./1.0_submit_batch.sh 3rdbatch
## (Submits 2.1_SPAdes2bwtblst.sh as an array job)

set -euo pipefail

## BATCH argument passed to 2.1_SPAdes2bwtblst.sh via --export
BATCH="${1:?Usage: [CLEAN=1] ./1.0_submit_batch.sh <batch>   e.g. ./1.0_submit_batch.sh 3rdbatch}"
PROCESS_EXEC=/dfs10/bio/abarrons/CA_curto/isolates/process_executables
## ==========================================================================

CLEAN="${CLEAN:-0}"

RAW=/dfs10/bio/abarrons/CA_curto/isolates/${BATCH}/raw_seqs/seqs

## Count R1 files = number of samples = array size
N=$(ls "$RAW"/*READ1-Sequences.txt.gz 2>/dev/null | wc -l)

if [ "$N" -eq 0 ]; then
  echo "ERROR: no *READ1-Sequences.txt.gz files found in $RAW" >&2
  exit 1
fi

echo "Batch: $BATCH"
echo "Samples found: $N  ->  submitting --array=1-$N"
echo "Cleanup: CLEAN=$CLEAN$( [ "$CLEAN" = "1" ] && echo '  (prior outputs will be wiped; raw reads kept)' )"

sbatch  --export=ALL,BATCH="$BATCH",CLEAN="$CLEAN"  \
        --array=1-"$N" \
        ${PROCESS_EXEC}/2.1_SPAdes2bwtblst.sh
