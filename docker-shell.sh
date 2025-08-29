#!/bin/bash
set -Eeuo pipefail

t_start="$(TZ=America/New_York date +'%Y%m%d_%H%M%S')"
JOB_NAME="infercnv_azure_${t_start}"
N_PARALLEL=4
N_THREADS=4

usage() {
  echo "Usage: $0 -I DATA_INPUT_DIR -O OUTPUT_DIR [-N JOB_NAME] [-P N_PARALLEL] [-T N_THREADS]"
  echo "  -I   (required) Path to input data directory (h5ad/mtx)"
  echo "  -O   (required) Path to output directory (created if missing)"
  echo "  -N   Optional job name (default: infercnv_azure_<timestamp>)"
  echo "  -P   Optional process parallelism (default: 4)"
  echo "  -T   Optional threads per process (default: 4)"
  echo "  -h   Show this help"
}

# Parse flags
DATA_FOLDER=""
OUT_FOLDER=""
while getopts ":I:O:N:P:T:h" opt; do
  case "$opt" in
    I) DATA_FOLDER="$OPTARG" ;;
    O) OUT_FOLDER="$OPTARG" ;;
    N) JOB_NAME="$OPTARG" ;;
    P) N_PARALLEL="$OPTARG" ;;
    T) N_THREADS="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "Error: Unknown option -$OPTARG" >&2; usage; exit 1 ;;
    :)  echo "Error: Option -$OPTARG requires an argument." >&2; usage; exit 1 ;;
  esac
done
shift $((OPTIND-1))

# check required arguments
[[ -n "$DATA_FOLDER" && -n "$OUT_FOLDER" ]] || { echo "Error: -I and -O are required."; usage; exit 2; }

[[ -d "$DATA_FOLDER" ]] || { echo "Error: DATA_FOLDER does not exist: $DATA_FOLDER" >&2; exit 1; }
mkdir -p "$OUT_FOLDER"

# environment variables
export JOB_NAME DATA_FOLDER OUT_FOLDER N_PARALLEL N_THREADS

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
amlt run "$OUT_FOLDER"/job_infercnv.resolved.yml "$JOB_NAME" -d "$JOB_NAME" 
echo "[Local] Submitted AML job: $JOB_NAME"
echo "[Local] This program will check job status every 2 minutes..."

# auto-check status of the job
while true; do
    T="$(TZ=America/New_York date +'%Y%m%d_%H%M%S')"
    S="$(amlt status "$JOB_NAME")"
    if echo "$S" | grep -qE "queued|scheduling|preparing|starting|running"; then
        echo "[$T] status check"
        echo "-------most recent log-------"
        echo "$(amlt logs "$JOB_NAME" --tail 10)"
        echo "-----------------------------"
        sleep 120
    elif echo "$S" | grep -qE "failed|canceled"; then
        echo "[$T] job failed, exiting"
        exit 1
    elif echo "$S" | grep -qE "completed"; then
        echo "[$T] job completed, proceeding to download results"
        break
    else
        echo "[$T] status not recognized, exiting"
        exit 1
    fi
done

# download results from Azure Blob Storage
azcopy sync \
  "https://exvivocoldeastus.blob.core.windows.net/projects/Projects/broad_infercnv_out/${JOB_NAME}?<SAS_TOKEN>" \
  "$OUT_FOLDER" \
  --recursive --delete-destination=false