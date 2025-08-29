if (!requireNamespace("infercnv", quietly = TRUE)) {
  stop("The 'infercnv' package is not installed. Please run:\n\n  BiocManager::install('infercnv')\n")
}

library(infercnv)
library(argparse)
options(scipen = 100)

# note in this version we assume all non-malignant cells are "normal"
parser <- ArgumentParser()
parser$add_argument("--matrix", required=TRUE, help="Path to raw count matrix")
parser$add_argument("--annotations", required=TRUE, help="Path to cell annotations")
parser$add_argument("--gene_order", required=TRUE, help="Path to gene order file")
parser$add_argument("--out_dir", required=TRUE, help="Output directory")
parser$add_argument("--cutoff", type="numeric", default=0.1, help="Cutoff for CNV detection (default: 0.1)")
parser$add_argument("--denoise", type="logical", default=TRUE, help="Denoise the data (default: TRUE)")
parser$add_argument("--HMM", type="logical", default=TRUE, help="Use HMM (default: TRUE)")
parser$add_argument("--cluster_by_groups", type="logical", default=TRUE, help="Cluster by groups (default: TRUE)")
parser$add_argument("--num_threads", type="integer", default=4, help="Number of threads to use (default: 4)")
parser$add_argument("--ref_group_names", nargs="+", default=c("normal"), help="Reference group names (default: normal)")

args <- parser$parse_args()

# create object
# note currently all normal cells are not distinguished
infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = args$matrix,
  annotations_file  = args$annotations,
  delim             = "\t",
  gene_order_file   = args$gene_order,
  ref_group_names   = args$ref_group_names
)

# run infercnv
infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff = args$cutoff,
  out_dir = args$out_dir,
  cluster_by_groups = args$cluster_by_groups, 
  denoise = args$denoise,
  HMM = args$HMM,
  num_threads = args$num_threads
)