# ==========================================================
# WORKFLOW FOR THE IDENTIFICATION OF ANTIMICROBIAL TARGETS
#
# 1. Genome retrieval
# 2. Plasmid filtering
# 3. Topology prediction with DeepTMHMM
# 4. Transmembrane protein filtering
# 5. Human homology identification
# 6. Human non-homology filtering
# 7. DEG database retrieval
# 8. Essentiality  identification
# 9. Essentiality filtering
# 10. CD-HIT clustering
# 11. Cluster filtering
# 12.Functional annotation
# 13. Candidate proteins prioritization
# 14. Plots generation
# ==========================================================
set -euo pipefail
source config.sh

log() { echo "[$(date '+%H:%M:%S')] $1"
}
check_dependencies() {
    log "[CHECK] Verifying required tools..."
    for tool in "$@"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "[ERROR] Required tool not found: $tool"
            exit 1
        fi
    done
    log "Available dependency $*"    
}

EMAIL="${1:-}"
ACCESSIONS="${2:-}"
if [[ -z "$EMAIL" ]]; then
    read -p "Enter your email required for NCBI policy: " EMAIL
fi

if [[ -z "$ACCESSIONS" ]]; then
    read -p "Enter the desired Accession numbers (comma-separated): " ACCESSIONS
fi

# Initial inputs validation
if [[ "$EMAIL" != *"@"* ]]; then
    echo "[ERROR] Invalid email"
    exit 1
fi

if [[ -z "$ACCESSIONS" ]]; then
    echo "[ERROR] No accessions provided"
    exit 1
fi


download_ncbi() {
    log "[STEP 1] Downloading genomes"
    "$PY_BASE" scripts/download_ncbi.py \
    "$EMAIL" "$ACCESSIONS" \
     --bacteria_dir $BACTERIA_DIR \
     --human_dir $HUMAN_DIR
}

decompress_files() {
    log "Decompressing files"
    INPUT_DIR="$1"
    OUTPUT_DIR="$2"
    "$PY_BASE" scripts/decompress.py \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR"  
}

plasmid_filter() {
    log "[STEP 2] Filtering plasmid sequences..."
    "$PY_BASE" scripts/plasmid_filter.py \
    --input_dir "$DECOMPRESSED_BACTERIA_DIR" \
    --output_dir "$PLASMID_FILTERED_OUT" \
    --plasmid_filtered "$PLASMID_STATS"
}

run_deeptmhmm() {
    log "[STEP 3] Topology prediction..."
    "$PY_DEEPTMHMM" scripts/run_deeptmhmm.py \
      --input_dir "$PLASMID_FILTERED_OUT" \
      --chunks_dir "$CHUNK_DIR" \
      --output_dir "$PROTEIN_PREDICTION_OUT" \
      --deeptmhmm "$DEEPTMHMM"
}

transmembrane_filter() {
    log "[STEP 4] Transmembrane filtering..."
    "$PY_BASE" scripts/transmembrane_filter.py \
      --plasmid_dir "$PLASMID_FILTERED_OUT" \
      --topology_dir "$PROTEIN_PREDICTION_OUT" \
      --output_dir "$TRANSMEMBRANE_DIR" \
      --transmembrane_filtered "$TRANSMEMBRANE_STATS"
}

run_blastp() {
    QUERY_DIR="$1"
    DB_FASTA="$2"
    OUTPUT_DIR="$3"
    DB_DIR="$4"
    STEP_NAME="$5"
    check_dependencies blastp
    log "[STEP $STEP_NAME - Running BLASTp"
    "$PY_BASE" scripts/blastp.py \
        --query_dir "$QUERY_DIR" \
        --db_fasta "$DB_FASTA" \
        --output_dir "$OUTPUT_DIR" \
        --db_dir "$DB_DIR"
}

blastp_filter() {
    TSV_DIR="$1"
    INPUT_DIR="$2"
    OUTPUT_DIR="$3"
    STATS_FILE="$4"
    MODE="$5"
    STEP_NAME="$6"
    log "[STEP $STEP_NAME - BLASTp filtering ($MODE)"
    "$PY_BASE" scripts/blastp_filter.py \
        --tsv_dir "$TSV_DIR" \
        --input_dir "$INPUT_DIR" \
        --output_dir "$OUTPUT_DIR" \
        --stats_file "$STATS_FILE" \
        --mode "$MODE"
}

download_deg() {
    log "[STEP 7] Downloading DEG database..."
    "$PY_BASE" scripts/download_deg.py \
       --url "$ESSENTIAL_URL" \
       --output_dir "$DEG_DIR"
}

merge_fastas() {
    log "Merging FASTA files"
    "$PY_BASE" scripts/merged.py \
        --input_dir "$ESSENTIALITY_DIR" \
        --output_dir "$CD_HIT_IN"
}

run_cdhit() {
    log "[STEP 10] CD-HIT clustering"
    check_dependencies cd-hit
    "$PY_BASE" scripts/cd_hit.py \
        --input_dir "$CD_HIT_IN/merged.faa" \
        --output_dir "$CD_HIT_OUT" \
        --cluster_threshold 0.9 \
        --cluster_size 5
}

filter_clusters() {
    log "[STEP 11] Filtering clusters"
    read -p "Minimum number of proteins required per cluster [default: 5]: " MIN_PROTEINS
    MIN_PROTEINS=${MIN_PROTEINS:-5}
    read -p "Minimum number of organisms per cluster to ensure conservation [default: 5]: " MIN_ORGANISMS
    MIN_ORGANISMS=${MIN_ORGANISMS:-5}
    "$PY_BASE" scripts/cluster_filter.py \
        --input_dir "$CLUSTERS_RAW" \
        --output_dir "$CLUSTERS_FILTERED" \
        --stats_file "$CLUSTERS_STATS" \
        --min_proteins "$MIN_PROTEINS" \
        --min_organisms "$MIN_ORGANISMS" \
        --finalcluster_dir "$CD_HIT_OUT" 
}

run_eggnog() {
    log "[STEP 12] Running eggNOG-mapper for functional annotation"
    check_dependencies "$(dirname "$PY_EGGNOG")/emapper.py"
    "$PY_EGGNOG" scripts/eggNOG_mapper.py \
        --input_fasta "$CD_HIT_OUT/representatives.faa" \
        --output_dir "$EGGNOG_OUT" \
        --db_dir "$DECOMPRESSED_EGGNOG_DIR" \
        --db_type "2" \
        --emapper "envs/eggnog_env/bin/emapper.py"
}

run_prioritization() {
    log "[STEP 13] Prioritizing antimicrobial candidates..."
    read -p "Enter GO pre-filtering threshold [default: 5]: " GO_THRESHOLD
    GO_THRESHOLD=${GO_THRESHOLD:-5}
    "$PY_BASE" scripts/prioritization.py \
        --input_file "$EGGNOG_OUT/eggnog_annotation.emapper.annotations" \
        --output_file "$PRIORITIZATION_OUT/candidates.tsv" \
        --go_threshold "$GO_THRESHOLD"
}

run_plots() {
    log "[STEP 14] Generating plots..."
    "$PY_BASE" analysis/plot.py \
        --stats_dir "$STATS_DIR" \
        --plots_dir "$PLOTS_DIR"
}


main() {
    log "Pipeline started"
    download_ncbi
    decompress_files "$BACTERIA_DIR" "$DECOMPRESSED_BACTERIA_DIR"
    decompress_files "$HUMAN_DIR" "$DECOMPRESSED_HUMAN_DIR"
    plasmid_filter
    run_deeptmhmm
    transmembrane_filter
    run_blastp "$TRANSMEMBRANE_DIR" "$DECOMPRESSED_HUMAN_DIR" "$BLAST_NONHOMOLOGOUS_OUT" "$BLAST_NONHOMOLOGOUS_DB" "5] Non-homology"
    blastp_filter "$BLAST_NONHOMOLOGOUS_OUT" "$TRANSMEMBRANE_DIR" "$NONHOMOLOGOUS_DIR" "$NONHOMOLOGOUS_STATS" "exclude" "6] Non-homology filtering"
    download_deg 
    decompress_files "$DEG_DIR" "$DECOMPRESSED_DEG_DIR"
    run_blastp "$NONHOMOLOGOUS_DIR" "$DECOMPRESSED_DEG_DIR" "$BLAST_ESSENTIALITY_OUT" "$BLAST_ESSENTIALITY_DB" "8] Essentiality"
    blastp_filter "$BLAST_ESSENTIALITY_OUT" "$NONHOMOLOGOUS_DIR" "$ESSENTIALITY_DIR" "$ESSENTIALITY_STATS" "include" "9] Essentiality filtering"
    merge_fastas
    run_cdhit
    filter_clusters
    decompress_files "$EGGNOG_DIR" "$DECOMPRESSED_EGGNOG_DIR"
    run_eggnog
    run_prioritization
    run_plots
    log "Pipeline finished"
}


main "$@"