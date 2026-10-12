# Curtobacterium Isolates Genomics Pipeline

This repository contains the SLURM HPC pipeline used to process, quality-filter, assemble, and analyze *Curtobacterium* bacterial isolate genomes from raw Illumina sequencing reads through to phylogenetics, comparative genomics, and carbohydrate-active enzyme (CAZyme) annotation.

---

## Pipeline Overview

The pipeline is organized in sequential stages indicated by file numbering:

```text
00_setup_taxdump.sh
1.0_submit_batch.sh  ──>  2.1_SPAdes2bwtblst.sh
1.1_submit_btkfilter.sh  ──>  2.2_blobtools2_filter.sh
3.0_consolidate_assemblies.sh
4.0_checkM2.sh
         │
         ▼ (Quality curation into all_good_genomes/)
5.1a_prodigal2hmmscan_nucl.sh  ──>  5.2a_clustalo2RAxML_nucl.sh   (Nucleotide tree)
5.1b_prodigal2hmmscan_prot.sh  ──>  5.2b_clustalo2RAxML_prot.sh   (Protein tree)
6.1_fastANI.sh
6.2_fastAAI.sh
7.0_run_dbCAN.sh
```

---

## Workflow Stages

### Stage 0: Database Setup
- **`00_setup_taxdump.sh`**  
  One-time setup script to download and unpack the NCBI `new_taxdump` required for BlobToolKit v2 taxonomy assignments.

### Stage 1 & 2: Assembly & Taxonomic Filtering (Per Batch)
- **`1.0_submit_batch.sh`**  
  Launcher that auto-counts read pairs for a given sequencing batch and submits `2.1_SPAdes2bwtblst.sh` as a SLURM job array.
- **`2.1_SPAdes2bwtblst.sh`**  
  Assembles raw reads using SPAdes, indexes and maps reads back to scaffolds using Bowtie2 and Samtools, and runs BLAST against the NCBI nt database.
- **`1.1_submit_btkfilter.sh`**  
  Launcher that auto-counts assemblies in a batch and submits `2.2_blobtools2_filter.sh` as a SLURM job array.
- **`2.2_blobtools2_filter.sh`**  
  BlobToolKit (BlobTools2) QC and taxonomic decontamination. Assesses dominant genus per scaffold set and applies length ($\ge 5\,\text{kb}$), coverage ($\ge 3\times$), and conditional GC filters to isolate target *Curtobacterium* contigs.

### Stage 3 & 4: Consolidation & Quality Assessment
- **`3.0_consolidate_assemblies.sh`**  
  Consolidates filtered assemblies across all sequencing batches into a central directory. Detects cross-batch sample ID duplicates, routes duplicates safely to a separate folder, and generates `consolidation_report.tsv`.
- **`4.0_checkM2.sh`**  
  Runs CheckM2 across assemblies to estimate genome completeness and contamination, outputting `quality_report.tsv`.

---

### Transition from Stage 4 to Stage 5: Curating `all_good_genomes/`

> **Note on Handover:** Stages 5, 6, and 7 read from `$BASE/all_good_genomes/`.  
> Assemblies passing quality criteria (e.g., completeness $\ge 90\%$, contamination $< 5\%$) are curated into this working folder as `*.fna` files. Any reference strains or outgroups (such as *Frigobacterium*) used for phylogenetic rooting or comparative genomics are also placed here.

---

### Stage 5: Phylogenomics (Core-Gene Concatenation Trees)
Core genes are identified from HMM profiles via Prodigal and HMMER (`hmmsearch`), aligned with Clustal Omega, concatenated with `catfasta2phyml`, and reconstructed with RAxML:
- **Nucleotide Tree:**
  - `5.1a_prodigal2hmmscan_nucl.sh`: Extracts nucleotide core gene coding sequences.
  - `5.2a_clustalo2RAxML_nucl.sh`: Aligns and concatenates nucleotide markers, then infers a maximum-likelihood phylogeny using RAxML (`GTRGAMMA`).
- **Protein Tree:**
  - `5.1b_prodigal2hmmscan_prot.sh`: Extracts translated amino acid core gene sequences.
  - `5.2b_clustalo2RAxML_prot.sh`: Aligns and concatenates protein markers, then infers a maximum-likelihood phylogeny using RAxML (`PROTGAMMAWAG`).

### Stage 6: Pairwise Genomic Distances
- **`6.1_fastANI.sh`**  
  Calculates whole-genome Average Nucleotide Identity (ANI) across isolate genomes.
- **`6.2_fastAAI.sh`**  
  Calculates whole-genome Average Amino Acid Identity (AAI) across isolate genomes using FastAAI.

### Stage 7: Functional Annotation
- **`7.0_run_dbCAN.sh`**  
  Annotates Carbohydrate-Active Enzymes (CAZymes) using `run_dbcan`, performs GFF preprocessing, and identifies CAZyme Gene Clusters (CGCs).

---

## Subdirectories

- **`NCBI_sub_execs/`**  
  Scripts for NCBI submission packaging, metadata mapping, and NCBI Foreign Contamination Screen (FCS) corrections (e.g., `4.4_fcs_fix.py`).
- **`legacy/`**  
  Earlier iterations, alternative tree-building methods (FastTree, bcgTree), and exploration scripts preserved for reproducibility and provenance.

---

## Software Dependencies & Environment

The pipeline relies on SLURM on an HPC cluster and the following tools:
- **Assembly & Mapping:** SPAdes ($\ge 3.15$), Bowtie2, Samtools, BBMap
- **Taxonomic Screening:** BlobToolKit (v2), BLAST+
- **Quality Assessment:** CheckM2
- **Phylogenetics:** Prodigal ($\ge 2.6.3$), HMMER ($\ge 3.3$), Clustal Omega, `catfasta2phyml`, RAxML ($\ge 8.2$)
- **Comparative Genomics:** fastANI, FastAAI
- **Functional Annotation:** `run_dbcan` ($\ge 4.x$)

*Note: Update base directory paths (`BASE`, `PROGRAM_PATH`, database paths) in script headers to match your cluster environment prior to job submission.*
