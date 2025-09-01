#!/bin/bash
# set -e pipefail

t_start="$(TZ=America/New_York date +'%Y%m%d_%H%M%S')"
JOB_NAME="infercnv_azure_${t_start}"
N_PARALLEL=4
N_THREADS=4

HMM=true
DENOISE=true
CUTOFF=0.1
REF_GROUP_NAMES="normal"

usage() {
  echo "Usage: $0 -I DATA_INPUT_DIR -O OUTPUT_DIR [-N JOB_NAME] [-P N_PARALLEL] [-T N_THREADS]"
  echo "  -I   (required) Path to input data directory (h5ad/mtx)"
  echo "  -O   (required) Path to output directory (created if missing)"
  echo "  -N   Optional job name (default: infercnv_azure_<timestamp>)"
  echo "  -P   Optional process parallelism (default: 4)"
  echo "  -T   Optional threads per process (default: 4)"
  echo "  -R   Optional reference group names, space-separated (default: normal)"
  echo "  -h   Show this help"
}

# Parse flags
DATA_FOLDER=""
OUT_FOLDER=""
while getopts ":I:O:N:P:R:T:M:h" opt; do
  case "$opt" in
    I) DATA_FOLDER="$OPTARG" ;;
    O) OUT_FOLDER="$OPTARG" ;;
    N) JOB_NAME="$OPTARG" ;;
    P) N_PARALLEL="$OPTARG" ;;
    T) N_THREADS="$OPTARG" ;;
    R) REF_GROUP_NAMES="$OPTARG" ;;
    M) MALIG_NAME="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Error: Unknown option -$OPTARG" >&2; usage; exit 1 ;;
    :)  echo "Error: Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done
shift $((OPTIND-1))

# check required arguments
[[ -n "$DATA_FOLDER" && -n "$OUT_FOLDER" ]] || { echo "Error: -I and -O are required."; usage; exit 1; }

[[ -d "$DATA_FOLDER" ]] || { echo "Error: DATA_FOLDER does not exist: $DATA_FOLDER" >&2; exit 1; }
mkdir -p "$OUT_FOLDER"

# environment variables
export JOB_NAME DATA_FOLDER OUT_FOLDER N_PARALLEL N_THREADS HMM DENOISE CUTOFF REF_GROUP_NAMES MALIG_NAME

# begin and make yml file
echo "[Local] Start time: $t_start"
envsubst < job_infercnv.yml.tmpl > "$OUT_FOLDER"/job_infercnv.resolved.yml
echo "[Local] Generated AML job YAML at $OUT_FOLDER/job_infercnv.resolved.yml"

# upload data to Azure Blob Storage
azcopy sync "$DATA_FOLDER" \
  "https://exvivocoldeastus.blob.core.windows.net/data/broad_infercnv_data/${JOB_NAME}?<SAS_TOKEN>" \
  --recursive --delete-destination=false
echo "[Local] Uploaded data from $DATA_FOLDER to Azure Blob Storage."

# sh file that loads necessary paths into yml
amlt run "$OUT_FOLDER"/job_infercnv.resolved.yml --replace "$JOB_NAME" -d "$JOB_NAME" 
echo "[Local] Submitted AML job: $JOB_NAME"
echo "[Local] This program will check job status every 2 minutes..."

# auto-check status of the job
while true; do
    T="$(TZ=America/New_York date +'%Y%m%d_%H%M%S')"
    S="$(amlt status "$JOB_NAME")"
    if echo "$S" | grep -qE "queued|scheduling|preparing|starting|running"; then
        echo "[$T] status check"
        # echo "-------most recent log-------"
        # { amlt logs "$JOB_NAME" 2>/dev/null || true; } | tail -n 20
        # echo "-----------------------------"
        sleep 120
    elif echo "$S" | grep -qE "failed|canceled"; then
        echo "[$T] job failed, exiting"
        CONT=0
        break
    elif echo "$S" | grep -qE "completed|pass"; then
        echo "[$T] job completed, proceeding to download results"
        CONT=1
        break
    else
        echo "[$T] status not recognized, exiting"
        CONT=0
        break
    fi
done

if [[ $CONT -eq 0 ]]; then
    echo "[Local] Job did not complete successfully. Please check the AML portal for details."
else
  # download results from Azure Blob Storage
  azcopy sync \
    "https://exvivocoldeastus.blob.core.windows.net/projects/Projects/broad_infercnv_out/${JOB_NAME}?<SAS_TOKEN>" \
    "$OUT_FOLDER" \
    --recursive --delete-destination=false
fi