#!/bin/bash
#set -e

mkdir -p user_runs

# environment variables
export SECRETS_DIR=$(pwd)/../secrets
export JOB_NAME="infercnv_test_$(TZ=America/New_York date +'%Y%m%d_%H%M%S')"
export DATA_FOLDER="/home/b-harryli/workspace/tmp_infercnv_adata"
export N_PARALLEL=4
export OUT_FOLDER="user_runs/${JOB_NAME}"

mkdir -p "$OUT_FOLDER"
envsubst < job_infercnv.yml.tmpl > "$OUT_FOLDER"/job_infercnv.resolved.yml

# sh file that loads necessary paths into yml
amlt run "$OUT_FOLDER"/job_infercnv.resolved.yml -d "$JOB_NAME"