#############################################################
####1 Reads quality control with bbduk (v38.90)##############
#############################################################

#!/bin/bash
# ==============================================================================
# Script for FastQ Quality Control and Adapter Trimming using BBDuk
# ==============================================================================

# Exit immediately if a command exits with a non-zero status,
# if an uninitialized variable is used, or if any command in a pipeline fails.
set -euo pipefail

# ------------------------------------------------------------------------------
# 1. DEFINE PATH VARIABLES
# ------------------------------------------------------------------------------
# Define input/output file paths and QC statistics file
R1_IN="path/to/R1.fastq.gz"
R2_IN="path/to/R2.fastq.gz"
R1_OUT="path/to/R1_clean.fastq.gz"
R2_OUT="path/to/R2_clean.fastq.gz"
STATS="path/to/stats.txt"

bbduk.sh \
  in="${R1_IN}" \
  in2="${R2_IN}" \
  out="${R1_OUT}" \
  out2="${R2_OUT}" \
  stats="${STATS}" \
  ftm=5 \
  tpe \
  tbo \
  qtrim=rl \
  trimq=25 \
  minlen=50 \
  ref=adapters,phix \
  -Xmx100g

#############################################################
####2 remove host contamination using BBSplit (v38.90)#######
#############################################################
#!/bin/bash
# ==============================================================================
# Script for Reference-Based Read Splitting (Mapping) using BBSplit
# ==============================================================================

# Exit immediately if a command exits with a non-zero status,
# if an uninitialized variable is used, or if any command in a pipeline fails.
set -euo pipefail

# ------------------------------------------------------------------------------
# 1. DEFINE PATH AND PARAMETER VARIABLES
# ------------------------------------------------------------------------------
# Define input/output file paths, reference genomes, and statistics file
IN1="path/to/input_R1.fastq.gz"
IN2="path/to/input_R2.fastq.gz"
OUTU1="path/to/unmapped_R1.fastq.gz"
OUTU2="path/to/unmapped_R2.fastq.gz"
REF="path/to/reference_folder_or_fastas" # here we use human reference genome (CHM13) as host reference genome
BASENAME="path/to/mapped_output_%.fastq.gz"
REFSTATS="path/to/mapping_stats.txt"

bbsplit.sh \
  in1="${IN1}" \
  in2="${IN2}" \
  outu1="${OUTU1}" \
  outu2="${OUTU2}" \
  ref="${REF}" \
  basename="${BASENAME}" \
  refstats="${REFSTATS}" \
  t=45 \
  -Xmx100g

###################################################################################
####3 taxonomic profiling use MetaPhlAn4.0 (v4.0.6, database: vJun23_202403)#######
###################################################################################

#!/bin/bash
# ==============================================================================
# Script for Taxonomic Profiling using MetaPhlAn
# ==============================================================================

# Exit immediately if a command exits with a non-zero status,
# if an uninitialized variable is used, or if any command in a pipeline fails.
set -euo pipefail

# ------------------------------------------------------------------------------
# 1. DEFINE PATH AND PARAMETER VARIABLES
# ------------------------------------------------------------------------------
# Define input alignment, output abundance file, and database directory
BOWTIE2OUT="path/to/sample_bowtie2.bz2"   # Input BowTie2 alignment output
ABUT2="path/to/profiled_metagenome.txt"   # Output taxonomic abundance profile
REF="path/to/metaphlan_database/"         # Path to the MetaPhlAn database


metaphlan \
  "${BOWTIE2OUT}" \
  --input_type bowtie2out \
  -t rel_ab \
  --nproc 40 \
  --unclassified_estimation \
  -o "${ABUT2}" \
  --bowtie2db "${REF}"
