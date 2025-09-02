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
parser$add_argument("--ref_group_names", required=TRUE, nargs="+", help="Reference group names (default: normal)")

# optional parameters
parser$add_argument("--cutoff", type="numeric", default=0.1, help="Cutoff for CNV detection (default: 0.1)")
parser$add_argument("--analysis_mode", type="character", default="subclusters", help="Analysis mode (default: subclusters)")
parser$add_argument("--cluster_by_groups", type="logical", default=FALSE, help="Cluster by groups (default: FALSE)")
parser$add_argument("--denoise", type="logical", default=TRUE, help="Denoise the data (default: TRUE)")
parser$add_argument("--tumor_subcluster_partition_method", type="character", default="random_trees", help="Tumor subcluster partition method (default: random_trees)")
parser$add_argument("--tumor_subcluster_pval", type="numeric", default=0.05, help="Tumor subcluster p-value (default: 0.05)")
parser$add_argument("--num_threads", type="integer", default=2, help="Number of threads to use (default: 4)")
parser$add_argument("--sd_amplifier", type="numeric", default=1.5, help="Standard deviation amplifier (default: 1.5)")
parser$add_argument("--noise_logistic", type="logical", default=TRUE, help="Use noise logistic (default: TRUE)")
parser$add_argument("--HMM", type="logical", default=FALSE, help="Use HMM (default: FALSE)")
parser$add_argument("--BayesMaxPNormal", type="numeric", default=0.2, help="Bayes maximum P(Normal) (default: 0.2)")

args <- parser$parse_args()

# create object
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
  analysis_mode= args$analysis_mode,
  out_dir = args$out_dir,
  cluster_by_groups = args$cluster_by_groups, 
  denoise = args$denoise,
  tumor_subcluster_partition_method = args$tumor_subcluster_partition_method,
  tumor_subcluster_pval = args$tumor_subcluster_pval,
  num_threads = args$num_threads,
  sd_amplifier = args$sd_amplifier,
  noise_logistic = args$noise_logistic,
  HMM = args$HMM,
  BayesMaxPNormal = args$BayesMaxPNormal
)
