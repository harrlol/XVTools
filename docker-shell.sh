#!/bin/bash
#set -e

# environment variables
export BASE_DIR=$(pwd)
export SECRETS_DIR=$(pwd)/../secrets
export JOB_NAME="infercnv_test_$(TZ=America/New_York date +'%Y%m%d_%H%M%S')"
export DATA_FOLDER="/home/b-harryli/workspace/tmp_infercnv_adata"
export N_PARALLEL=4
export OUT_FOLDER="tmp_out_8_27"


envsubst < job_infercnv.yml.tmpl > job_infercnv.resolved.yml

# sh file that loads necessary paths into yml
amlt run job_infercnv.resolved.yml -d "$JOB_NAME"