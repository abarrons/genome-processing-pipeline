#!/bin/bash

#SBATCH --job-name=prodigal2hmmscan  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=12:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --cpus-per-task=4               ## number of CPUs to be used per task
#SBATCH --mem-per-cpu=4GB
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid
#SBATCH --array=1-395%5                   ## This value (1-64) depends on the number of samples (1-no. of samples)

## Input: Expects high-quality passing assemblies and any reference/outgroup
## genomes in $WD ($BASE/all_good_genomes) as *.fna files.

## This sets environmnetal variables to be used by slurm in order to use all possible RAM
## Set environmental variables
BASE=/dfs10/bio/abarrons/CA_curto/isolates
REFDB=/dfs10/bio/abarrons/core_genes_db/hmmprofiles
WD=$BASE/all_good_genomes
## Load programs: prodigal
module load prodigal/2.6.3

## Go to working directory and make folders for output files
cd $WD

mkdir -p "$WD/translated_seqs" "$WD/hmmscan_out/filtered_seqs" "$WD/hmmscan_out/seqs_ids"

## Set environmental variables specific for the job
mapfile -t FASTAS < <( find . -maxdepth 1 -type f -name '*.fna' | sort )
SAMPLE_NAME="${FASTAS[$SLURM_ARRAY_TASK_ID-1]}"
ISOLATE_ID=$(basename "$SAMPLE_NAME" .fna)


## Prodigal will translate nucleotide sequences to protein sequences for subsequent analysis
prodigal \
  -a $WD/translated_seqs/${ISOLATE_ID}.faa \
  -d $WD/translated_seqs/${ISOLATE_ID}_nucleotid_seqs.fna \
  -i $SAMPLE_NAME


## Load programs: HMMER
module load hmmer/3.3

#get the list of proteins
cd $REFDB
proteins=$(for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done)

## Go to working directory and make folder for output files
cd $WD/translated_seqs
mkdir $WD/hmmscan_out


## hmmscan will identify the protien sequences that mache to the core gene sequences in the .hmm profiles (downloaded from Alex's github: elevation_community)
## make a while read line loop to go through each gene

for protein in $proteins
  do
        hmmsearch --tblout $WD/hmmscan_out/${ISOLATE_ID}.${protein}.hmm.txt -E 1e-10 --cpu 4 \
              $REFDB/${protein}.hmm ${ISOLATE_ID}.faa
  done


#load BBmap suit we eill need it to filter the fasta files
module load bbmap/38.96


## Go to working directory and make folder for output files
cd $WD/hmmscan_out
mkdir filtered_seqs
mkdir seqs_ids

## make a loop to go throughg all the combination of genomes and proteins
## to extract the sequence identifiers for each protein in each genome and
## get that sequence from the protein fasta file of that genome

#first go to the hmmscan outp folder

for protein in $proteins
  do

    # first get the highest socore match for each protein and store the protein name and sequenceID
    # to filter that sequence from the protein fasta file
      hmmscanout=$(cat ${ISOLATE_ID}.${protein}.hmm.txt | awk '$1 !~ /^#/ && ($6 > max_score || NR == 1) { max_score = $6; best_entry = $0 } END { split(best_entry, fields); prot_name = fields[3]; seqID = fields[1]; print prot_name, seqID }')

    # Extract the variables from the output
      prot_name=$(echo "$hmmscanout" | awk '{print $1}')
      seqID=$(echo "$hmmscanout" | awk '{print $2}')

      echo "$seqID" > $WD/hmmscan_out/seqs_ids/${ISOLATE_ID}.${protein}.tmp.txt

    # now let's use those variables to get the mathcing sequence from each genome
    # we will use filterbyname fromm the BBmap tool suite
      filterbyname.sh \
        in=$WD/translated_seqs/${ISOLATE_ID}_nucleotid_seqs.fna \
        out=$WD/hmmscan_out/filtered_seqs/${ISOLATE_ID}.${prot_name}.fna \
        names=$WD/hmmscan_out/seqs_ids/${ISOLATE_ID}.${protein}.tmp.txt \
        ow=t \
        include=t \
        fixjunk


    # Replace header in the FASTA file
      fasta_file="filtered_seqs/${ISOLATE_ID}.${prot_name}.fna"
      tmp_file="filtered_seqs/${ISOLATE_ID}.${prot_name}.temp.fna"

   # Use awk to replace specific parts of the header
   awk -v var1="$ISOLATE_ID" '/^>/{ sub(">[^ ]*", ">" var1); } 1' "$fasta_file" > "$tmp_file"

   # Rename the temporary file to overwrite the original FASTA file
   mv "$tmp_file" "$fasta_file"

 done

#rm $WD/tmp.txt
