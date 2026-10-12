#!/bin/bash

#SBATCH --job-name=clustalo2raxml  ## Name of the job.
#SBATCH -A JMARTINY_LAB                ## Account to charge; personal or lab
#SBATCH -p standard                    ## Partition/queue name
#SBATCH --time=36:00:00
#SBATCH --nodes=1                      ## (-N) number of nodes the job will use
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=32              ## number of CPUs to be used per task
#SBATCH -o %x.%A.%a.out             ## File to which STDOUT will be written, %j inserts jobid
#SBATCH -e %x.%A.%a.err             ## File to which STDERR will be written, %j inserts jobid

# Set folder paths
BASE=/dfs10/bio/abarrons/CA_curto/isolates/all_good_genomes
WD=$BASE/hmmscan_out
REFDB=/dfs10/bio/abarrons/core_genes_db/hmmprofiles
PROGRAM_PATH=/dfs10/bio/abarrons/programs
CONDA_PATH=/data/homezvol0/abarrons/.conda/envs/catfasta2phyml/bin
THREADS=32

#get the list of proteins
cd $REFDB
proteins=$(for f in *.hmm; do printf '%s\n' "${f%.hmm}"; done)


cd $WD/filtered_seqs

# make output directory for alignmets
rm -rf $BASE/clustalo_out
mkdir $BASE/clustalo_out

  # Iterate over the protein names to make an alignment for each protein
  for protein in $proteins; do
    # Find files matching the string
    files=$(find . -name "*${protein}*" | cut -c 3-)

    # Create a temporary file for the concatenated sequences
    tmp_file=$(mktemp)

    # Concatenate sequences from different files into a single multi-sequence FASTA file
    cat $files > $tmp_file

    # Perform alignment with clustalo if files are found
    if [[ -s $tmp_file ]]; then
      $PROGRAM_PATH/clustalo -i "$tmp_file" -o "$BASE/clustalo_out/${protein}.aln" --threads=$THREADS -v
    else
      echo "No files found for string: $protein"
    fi

    rm $tmp_file
    #rm $files
  done


cd $BASE/clustalo_out

#do a clean up of the folder

rm *.fix.*



# Iterateq  over the genome names
for protein in $proteins; do

  # go over each alingment file to fix the header with seqkit replace

    echo "Fixing headers for ${protein}.aln"
  # Extract the sequences belonging to the genome and save them to the output file
  # using seqkit toolsuite, we will use subcomand grep
    $PROGRAM_PATH/seqkit replace -p "(.+)${protein}(.+)" -r "$1" ${protein}.aln | awk '/^>/ { sub(/ .*/, ""); } 1' > ${protein}.fix.aln

done

# now run catfasta2phyml to concatenate the aligned sequences for the phylogenetic tree

rm concat_alignment.fna

$CONDA_PATH/catfasta2phyml --fasta --concatenate --verbose *.fix.aln > concat_alignment.fna

module load openmpi/4.0.3/gcc.8.4.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

$PROGRAM_PATH/standard-RAxML/raxmlHPC-PTHREADS \
  -f a \
  -m GTRGAMMA \
  -p 4321 \
  -s concat_alignment.fna \
  -n final21mrkrgn -T $THREADS -x 100 -N 100
