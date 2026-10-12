#!/bin/bash

#SBATCH --job-name=db_CAN      ## Name of the job.
#SBATCH -A JMARTINY_LAB             ## Account to charge; personal or lab
#SBATCH -p standard                 ## Partition/queue name
#SBATCH --time=24:00:00
#SBATCH --nodes=1                   ## (-N) number of nodes the job will use
#SBATCH --cpus-per-task=8          ## number of CPUs to be used per task
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid
#SBATCH --array=1-394               ## This value (1-64) depends on the number of samples i.e. 1-num. samples

export BASE=/dfs10/bio/abarrons/CA_curto
export WD=$BASE/isolates/all_good_genomes
PROGRAM_PATH=/data/homezvol0/abarrons/.conda/envs/run_dbcan_env/bin
export db_can_db_PATH=/data/homezvol0/abarrons/db
OUTDIR=$WD/db_CAN_out
THREADS=8

cd $WD

## load conda modules

module load anaconda/2024.06

conda init bash
conda activate run_dbcan_env



## get sample names and index each file for array task
mapfile -t FASTAS < <( find . -maxdepth 1 -type f -name '*.fna' ! -name '*Frigo*' | sort )
SAMPLE_NAME="${FASTAS[$SLURM_ARRAY_TASK_ID-1]}"
ISOLATE_ID=$(basename "$SAMPLE_NAME" .fna)

mkdir -p "$OUTDIR/$ISOLATE_ID"


# Run CAZyme annotation
$PROGRAM_PATH/run_dbcan CAZyme_annotation \
  --input_raw_data $SAMPLE_NAME \
  --mode prok \
  --output_dir $OUTDIR/$ISOLATE_ID \
  --db_dir $db_can_db_PATH \
  --threads $THREADS


  # ─── 4. CGC info generation (gff_process) ───────────────────────────────────────
  # Assumes that after your initial run_dbcan you have:
  #   $WORKDIR/uniInput.gff   (the GFF produced by the CAZyme annotation step)
  $PROGRAM_PATH/run_dbcan \
    gff_process \
      --input_gff "$OUTDIR/$ISOLATE_ID/uniInput.gff" \
      --gff_type prodigal \
      --db_dir "$db_can_db_PATH" \
      --output_dir "$OUTDIR/$ISOLATE_ID" \
      --threads $THREADS

  #____5. CGC finder____

  $PROGRAM_PATH/run_dbcan \
    cgc_finder \
      --output_dir "$OUTDIR/$ISOLATE_ID" \
      --use_null_genes \
      --num_null_gene 5 \
      --use_distance \
      --base_pair_distance 15000 \
      --additional_genes TC,TF,STP
