#!/usr/bin/env bash

# ==============================================================================
# LSDV Nanopore Assembly Pipeline
# Developer : Vinay Rajput (srm.vinay0005@gmail.com)
#
# A production-grade workflow for processing Oxford Nanopore reads:
# NanoPlot QC -> Host Depletion (Human/Vero/RSV) -> Flye Assembly -> 
# Medaka Polishing -> RagTag Scaffolding -> Mapping & Coverage Analysis -> QUAST
# ==============================================================================

# Strict error handling
set -euo pipefail
IFS=$'\n\t'

# ------------------------------------------------------------------------------
# GLOBAL CONSTANTS & DEFAULTS
# ------------------------------------------------------------------------------
SCRIPT_NAME=$(basename "$0")
THREADS=4
INPUT_DIR=""
REF_DIR=""
WORKDIR="LSDV_analysis"
GENOME_SIZE="151k"
MEDAKA_MODEL="r1041_e82_400bps_sup_v5.0.0"

# Required Dependencies
readonly DEPENDENCIES=(
    "NanoPlot" "minimap2" "samtools" "gzip" "zcat" 
    "flye" "medaka_consensus" "ragtag.py" "awk" "quast.py" "seqkit"
)

# ------------------------------------------------------------------------------
# LOGGING & ERROR HANDLING FUNCTIONS
# ------------------------------------------------------------------------------
log_info()    { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] \033[0;34m[INFO]\033[0m $*"; }
log_success() { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] \033[0;32m[SUCCESS]\033[0m $*"; }
log_warn()    { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] \033[0;33m[WARN]\033[0m $*" >&2; }
log_error()   { echo -e "[$(date +'%Y-%m-%d %H:%M:%S')] \033[0;31m[ERROR]\033[0m $*" >&2; }

# Production Cleanup Trap
cleanup() {
    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        log_error "Pipeline aborted prematurely due to an unexpected failure."
    fi
    exit $exit_code
}
trap cleanup EXIT

# ------------------------------------------------------------------------------
# HELP & USAGE MESSAGE
# ------------------------------------------------------------------------------
usage() {
    cat << EOF

================================================================================
LSDV NANOPORE ASSEMBLY PIPELINE — USAGE MANUAL
================================================================================

This automated pipeline runs structural QC, serial host depletion, de novo assembly,
polishing, and scaffolding on Oxford Nanopore data.

USAGE:
  $SCRIPT_NAME -i <input_dir> -r <ref_dir> [options]

REQUIRED ARGUMENTS:
  -i DIR    Path to directory containing input raw, gzipped Nanopore files (*.fastq.gz).
  -r DIR    Path to directory containing reference FASTA structural assets.
            The pipeline strictly expects the following file names inside this directory:
              ├── human_GRCh38.fa
              ├── Chlorocebus_sabaeus.fa
              ├── RSV_reference.fasta
              └── LSDV_reference.fasta

OPTIONAL OPTIONS:
  -w DIR    Output working directory path. (Default: "${WORKDIR}")
  -t INT    Number of CPU threads to allocate across execution steps. (Default: ${THREADS})
  -g STR    Estimated genome size flag passed directly into Flye. (Default: "${GENOME_SIZE}")
  -m STR    Polishing model parameter passed directly into Medaka. (Default: "${MEDAKA_MODEL}")
  -h        Display this detailed help configuration screen and exit.

EXECUTION EXAMPLES:
  # Basic execution using defaults
  $SCRIPT_NAME -i ./Fastq_files -r ./Host_references

  # Production execution tuning high-resource compute allocations
  $SCRIPT_NAME -i ./Fastq_files -r ./Host_references -w ./run_2026_output -t 32 -g 155k

DEVELOPER & ENVIRONMENTAL REQUIREMENTS:
  Before launching, ensure all core bioinformatics dependencies are accessible in
  your active \$PATH. If utilizing an isolated Conda deployment environment:
  
    $ conda activate lsdv_pipeline_env
    $ chmod +x $SCRIPT_NAME
    $ ./$SCRIPT_NAME -i <input_dir> -r <ref_dir>

================================================================================
EOF
    exit 1
}

# ------------------------------------------------------------------------------
# ARGUMENT PARSING
# ------------------------------------------------------------------------------
while getopts "i:r:w:t:g:m:h" opt; do
    case "${opt}" in
        i) INPUT_DIR="${OPTARG}" ;;
        r) REF_DIR="${OPTARG}" ;;
        w) WORKDIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        g) GENOME_SIZE="${OPTARG}" ;;
        m) MEDAKA_MODEL="${OPTARG}" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validate Required Inputs
if [ -z "${INPUT_DIR}" ] || [ -z "${REF_DIR}" ]; then
    log_error "Missing required parameters. Option flags '-i' and '-r' are mandatory configurations."
    usage
fi

# Define reference paths explicitly after target directory confirmation
HUMAN_REF="${REF_DIR}/human_GRCh38.fa"
VERO_REF="${REF_DIR}/Chlorocebus_sabaeus.fa"
RSV_REF="${REF_DIR}/RSV_reference.fasta"
LSDV_REF="${REF_DIR}/LSDV_reference.fasta"

# ------------------------------------------------------------------------------
# VALIDATION CHECKS (FAIL-EARLY DESIGN)
# ------------------------------------------------------------------------------
log_info "Initiating system pre-flight environment checks..."

# 1. Dependency Environment Verification
for tool in "${DEPENDENCIES[@]}"; do
    if ! command -v "$tool" &> /dev/null; then
        log_error "Missing critical environmental dependency: '${tool}' is not installed or accessible in PATH."
        exit 1
    fi
done

# 2. File and Folder Assertions
if [ ! -d "${INPUT_DIR}" ]; then
    log_error "Input directory does not exist: ${INPUT_DIR}"
    exit 1
fi

readonly REQUIRED_REFS=("${HUMAN_REF}" "${VERO_REF}" "${RSV_REF}" "${LSDV_REF}")
for ref in "${REQUIRED_REFS[@]}"; do
    if [ ! -f "$ref" ]; then
        log_error "Critical reference file asset missing: ${ref}"
        exit 1
    fi
done

# Find target sample fastq files safely
shopt -s nullglob
FASTQ_FILES=("${INPUT_DIR}"/*.fastq.gz)
shopt -u nullglob

if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    log_error "No valid long-read input archive targets found inside: ${INPUT_DIR}/*.fastq.gz"
    exit 1
fi

# ------------------------------------------------------------------------------
# PIPELINE SETUP & INDEXING
# ------------------------------------------------------------------------------
log_info "Initializing production workspace environments..."
mkdir -p "${WORKDIR}"

# Initialize global tracking output
STATS_FILE="${WORKDIR}/Pipeline_Statistics.tsv"
if [ ! -f "${STATS_FILE}" ]; then
    echo -e "Sample\tRawReads\tHumanRemoved\tVeroRemoved\tRSVRemoved\tRemainingReads\tRemainingPct\tAssemblySize(bp)\tContigs\tN50(bp)\tMeanCoverage" > "${STATS_FILE}"
fi

# Helper reader function
count_reads() {
    zcat "$1" | awk 'END{print (NR==0?0:NR/4)}'
}

# Building Minimap2 Indexes
log_info "Verifying and updating target minimap2 indices..."
for ref in "${REQUIRED_REFS[@]}"; do
    base_name=$(basename "$ref" | cut -d. -f1)
    mmi_idx="${WORKDIR}/${base_name}.mmi"
    if [ ! -f "$mmi_idx" ] || [ "$ref" -nt "$mmi_idx" ]; then
        log_info "Generating structural index mapping file: ${mmi_idx}"
        minimap2 -d "$mmi_idx" "$ref" 2>/dev/null
    fi
done

# ------------------------------------------------------------------------------
# CORE PIPELINE EXECUTION LOOP
# ------------------------------------------------------------------------------
for FASTQ in "${FASTQ_FILES[@]}"; do
    SAMPLE=$(basename "$FASTQ" .fastq.gz)
    SAMPLE_DIR="${WORKDIR}/${SAMPLE}"
    
    log_info "================================================================="
    log_info "Processing Target Sample: ${SAMPLE}"
    log_info "================================================================="
    
    mkdir -p "${SAMPLE_DIR}"
    
    # 1. Base Analysis Metrics
    log_info "Computing total raw molecule counts..."
    RAW_READS=$(count_reads "$FASTQ")
    log_info "Identified Raw Sequences: ${RAW_READS}"
    
    if [ "$RAW_READS" -eq 0 ]; then
        log_warn "Target payload file '${FASTQ}' is entirely empty. Skipping sample pipeline track."
        continue
    fi

    # 2. Quality Evaluation
    log_info "Executing NanoPlot Diagnostic Metrics..."
    NanoPlot --fastq "$FASTQ" -t "${THREADS}" -o "${SAMPLE_DIR}/nanoplot" -p "${SAMPLE}_" &>/dev/null

    # 3. Serial Biosphere Contamination Depletion
    CURRENT_FASTQ="$FASTQ"
    
    # Step 3a: Human Removal
    log_info "Depleting Human Host contamination..."
    minimap2 -ax map-ont -t "${THREADS}" "${WORKDIR}/human_GRCh38.mmi" "$CURRENT_FASTQ" 2>/dev/null \
        | samtools view -@ "${THREADS}" -b -o "${SAMPLE_DIR}/human.bam"
    HUMAN_REMOVED=$(samtools view -c -F 4 "${SAMPLE_DIR}/human.bam")
    samtools fastq -f 4 "${SAMPLE_DIR}/human.bam" 2>/dev/null | gzip > "${SAMPLE_DIR}/no_human.fastq.gz"
    CURRENT_FASTQ="${SAMPLE_DIR}/no_human.fastq.gz"

    # Step 3b: Vero Removal
    log_info "Depleting Vero Host contamination..."
    minimap2 -ax map-ont -t "${THREADS}" "${WORKDIR}/Chlorocebus_sabaeus.mmi" "$CURRENT_FASTQ" 2>/dev/null \
        | samtools view -@ "${THREADS}" -b -o "${SAMPLE_DIR}/vero.bam"
    VERO_REMOVED=$(samtools view -c -F 4 "${SAMPLE_DIR}/vero.bam")
    samtools fastq -f 4 "${SAMPLE_DIR}/vero.bam" 2>/dev/null | gzip > "${SAMPLE_DIR}/no_human_vero.fastq.gz"
    CURRENT_FASTQ="${SAMPLE_DIR}/no_human_vero.fastq.gz"

    # Step 3c: RSV Removal
    log_info "Depleting cross-contaminating RSV sequences..."
    minimap2 -ax map-ont -t "${THREADS}" "${WORKDIR}/RSV_reference.mmi" "$CURRENT_FASTQ" 2>/dev/null \
        | samtools view -@ "${THREADS}" -b -o "${SAMPLE_DIR}/rsv.bam"
    RSV_REMOVED=$(samtools view -c -F 4 "${SAMPLE_DIR}/rsv.bam")
    samtools fastq -f 4 "${SAMPLE_DIR}/rsv.bam" 2>/dev/null | gzip > "${SAMPLE_DIR}/${SAMPLE}.non_host.fastq.gz"
    
    # Filtered metrics output
    REMAINING_READS=$(count_reads "${SAMPLE_DIR}/${SAMPLE}.non_host.fastq.gz")
    REMAINING_PCT=$(awk -v r="$REMAINING_READS" -v t="$RAW_READS" 'BEGIN{printf "%.2f",(r/t)*100}')
    log_info "Depletion Summary -> Retained: ${REMAINING_READS} (${REMAINING_PCT}%)"

    if [ "$REMAINING_READS" -lt 10 ]; then
        log_warn "Insufficient target viral reads remain after depletion for standard processing. Skipping Assembly."
        continue
    fi

    # 4. De Novo Core Assembly
    log_info "Executing Flye De Novo Assembler..."
    flye --nano-hq "${SAMPLE_DIR}/${SAMPLE}.non_host.fastq.gz" --threads "${THREADS}" --genome-size "${GENOME_SIZE}" --out-dir "${SAMPLE_DIR}/flye" &>/dev/null

    if [ ! -f "${SAMPLE_DIR}/flye/assembly.fasta" ]; then
        log_error "Flye failed to generate a continuous sequence graph topology for ${SAMPLE}."
        continue
    fi

    # 5. Medaka Polishing
    log_info "Polishing structural draft alignments via Medaka Consensus Neural Network..."
    medaka_consensus -i "${SAMPLE_DIR}/${SAMPLE}.non_host.fastq.gz" -d "${SAMPLE_DIR}/flye/assembly.fasta" -o "${SAMPLE_DIR}/medaka" -t "${THREADS}" -m "${MEDAKA_MODEL}" &>/dev/null

    # 6. Scaffolding Refinement
    log_info "Scaffolding consensus drafts using RagTag reference tracks..."
    ragtag.py scaffold "${LSDV_REF}" "${SAMPLE_DIR}/medaka/consensus.fasta" -t "${THREADS}" -o "${SAMPLE_DIR}/ragtag" &>/dev/null
    FINAL_ASSEMBLY="${SAMPLE_DIR}/ragtag/ragtag.scaffold.fasta"

    # 7. Post-Production Validation Alignments
    log_info "Calculating processing maps & depth vectors against final consensus build..."
    minimap2 -ax map-ont -t "${THREADS}" "${FINAL_ASSEMBLY}" "${SAMPLE_DIR}/${SAMPLE}.non_host.fastq.gz" 2>/dev/null \
        | samtools sort -@ "${THREADS}" -o "${SAMPLE_DIR}/final.bam"
    samtools index "${SAMPLE_DIR}/final.bam"

    # 8. Statistical Calculations & Metric Extractions
    samtools depth -aa "${SAMPLE_DIR}/final.bam" > "${SAMPLE_DIR}/per_base_coverage.depth.txt"
    MEAN_COVERAGE=$(awk '{sum+=$3} END{print (NR==0?0:sprintf("%.2f",sum/NR))}' "${SAMPLE_DIR}/per_base_coverage.depth.txt")

    log_info "Evaluating structural output via QUAST Framework..."
    quast.py "${FINAL_ASSEMBLY}" -o "${SAMPLE_DIR}/quast" --threads "${THREADS}" &>/dev/null

    CONTIGS=$(grep -c "^>" "${FINAL_ASSEMBLY}" || true)
    ASSEMBLY_SIZE=$(seqkit stats -T "${FINAL_ASSEMBLY}" 2>/dev/null | awk 'NR==2{print $5}')
    N50=$(grep "^N50" "${SAMPLE_DIR}/quast/report.tsv" | cut -f2 || echo "N/A")

    # 9. Record Database Storage Outputs
    echo -e "${SAMPLE}\t${RAW_READS}\t${HUMAN_REMOVED}\t${VERO_REMOVED}\t${RSV_REMOVED}\t${REMAINING_READS}\t${REMAINING_PCT}\t${ASSEMBLY_SIZE}\t${CONTIGS}\t${N50}\t${MEAN_COVERAGE}" >> "${STATS_FILE}"

    cat > "${SAMPLE_DIR}/summary.txt" << EOF
Sample: ${SAMPLE}
Raw Reads: ${RAW_READS}
Human Removed: ${HUMAN_REMOVED}
Vero Removed: ${VERO_REMOVED}
RSV Removed: ${RSV_REMOVED}
Remaining Reads: ${REMAINING_READS}
Remaining Percentage: ${REMAINING_PCT}%
Assembly Size: ${ASSEMBLY_SIZE}
Contigs: ${CONTIGS}
N50: ${N50}
Mean Coverage: ${MEAN_COVERAGE}
Final Assembly Path: ${FINAL_ASSEMBLY}
Coverage Profile Path: ${SAMPLE_DIR}/per_base_coverage.depth.txt
EOF

    log_success "Processing cycle complete for sequence footprint: ${SAMPLE}"
done

log_success "All pipeline matrix tasks have concluded successfully."
column -t "${STATS_FILE}"
