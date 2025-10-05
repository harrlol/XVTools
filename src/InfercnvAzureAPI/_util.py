import pandas as pd
import matplotlib.pyplot as plt
import datetime
import numpy as np
import os.path as osp
import subprocess
import fnmatch
import os
import urllib.request
import glob
from tqdm import tqdm
from google.oauth2 import service_account
from google.cloud import storage


def auto_expand(path):
    """
    Detects and expands the path if it contains '~'.
    """
    if '~' in path:
        return osp.expanduser(path)
    return path

## Google Cloud Storage helpers ##
def parse_gs_path_helper(gs_path):
    assert gs_path.startswith("gs://"), "Not a valid GCS path"
    path_parts = gs_path[5:].split("/", 1)
    bucket = path_parts[0]
    prefix = path_parts[1] if len(path_parts) > 1 else ""
    return bucket, prefix


def download_terra_data(SECRET, gcs_dir, local_dir, 
                        SCOPES=['https://www.googleapis.com/auth/devstorage.read_only'], 
                        pattern=None, excluded_suffixes=None):
    """ 
    Downloads data from a Terra bucket to a local directory.
    Args:
        secret (str): Path to the service account JSON key file.
        gcs_dir (str): The GCS bucket directory to download from.
        local_dir (str): The local directory to save the downloaded files.
    """
    # prep gcs
    credentials = service_account.Credentials.from_service_account_file(
        SECRET,
        scopes=SCOPES
    )
    bucket, prefix = parse_gs_path_helper(gcs_dir)
    client = storage.Client(credentials=credentials)
    bucket = client.get_bucket(bucket)
    blobs = bucket.list_blobs()

    # prep local
    if not osp.exists(local_dir):
        os.makedirs(local_dir, exist_ok=True)
    if pattern is None:
        pattern = "*"
    if excluded_suffixes is None:
        excluded_suffixes = ()
    else:
        excluded_suffixes = tuple(excluded_suffixes)
    matched_blobs = [
        blob for blob in blobs
        if fnmatch.fnmatch(blob.name[len(prefix):], pattern)
        and not blob.name.endswith(excluded_suffixes)
    ]

    # scan and download
    for blob in tqdm(matched_blobs, desc="Downloading files", unit="file"):
        tqdm.write(f"→ {blob.name}")
        relative_path = blob.name[len(prefix):]
        local_file_path = os.path.join(local_dir, relative_path)
        os.makedirs(os.path.dirname(local_file_path), exist_ok=True)
        blob.download_to_filename(local_file_path)

    return

## anndata/infercnv helpers ##
def h5ad_to_matrix(adata):
    expr = adata.raw.X if adata.raw is not None else adata.X
    if not isinstance(expr, np.ndarray):
        expr = expr.toarray() 
    df = pd.DataFrame(expr.T,
                      index=adata.var_names,
                      columns=adata.obs_names)
    return df


def extract_infercnv_sampleann(adata, cell_type_col="cell_type", auto_infer=True):
    """ 
    Extracts sample annotations for InferCNV from an AnnData object.
    Will not distinguish between different normal cell types
    Assumes that we malignant cells are labeled with "malignant" in the specified cell_type_col.
    """

    labels = adata.obs[cell_type_col].astype(str)
    if auto_infer:
        is_malignant = labels.str.lower().str.contains("malignant")
        annotated_labels = labels.copy()
        annotated_labels[is_malignant] = "malignant"
        annotated_labels[~is_malignant] = "normal"
    else:
        # if auto_infer is False, we assume the labels are already in the correct format
        annotated_labels = labels

    sample_annotations = pd.DataFrame({
        'cell_id': adata.obs_names,
        'label': annotated_labels
    })

    return sample_annotations


def cnv_state_cnt(cluster_geno, save_path=None, show=False):
    """
    Count how many times each state appears for each grouping in a given infercnv result
    """
    level_counts = cluster_geno.apply(pd.Series.value_counts, axis=1).fillna(0).astype(int)
    level_counts.plot(kind="bar", stacked=True, figsize=(12, 6))
    plt.xlabel("Row")
    plt.ylabel("Count of Levels")
    plt.title("Count of Float Levels per Row")
    plt.legend(title="Level")
    plt.xticks(rotation=45)
    plt.tight_layout()
    if save_path:
        plt.savefig(auto_expand(save_path))
    if show:
        plt.show()
    plt.close()

    return


def store_infercnv_in_adata(adata, infercnv_out_path):
    """ 
    Store inferCNV results in the AnnData object.
    """

    # parse infercnv output
    pat = "HMM_CNV_predictions.*_genes.dat"
    infercnv_files = glob.glob(os.path.join(infercnv_out_path, pat))
    if len(infercnv_files) == 0:
        raise ValueError("No inferCNV files found in the specified directory. Please check the path and pattern.")
    
    # grouping to cnv
    infercnv_df = pd.read_csv(infercnv_files[0], sep="\t")
    sub_infercnv_df = infercnv_df[['cell_group_name', 'state', 'gene']]
    collapsed = sub_infercnv_df.groupby('cell_group_name').apply(lambda x: dict(zip(x['gene'], x['state']))).to_dict()
    cluster_geno = pd.DataFrame(collapsed).T
    cluster_geno.index = cluster_geno.index.str.split(".", n=1).str[1]

    # matching/expanding to adata genes, qc
    cluster_geno_expanded = cluster_geno.reindex(columns=adata.var_names, fill_value=np.nan)
    dropped_genes = cluster_geno.columns.difference(adata.var_names)
    if len(dropped_genes) > 0:
        print(f"[inferCNV] Dropped {len(dropped_genes)} genes not found in adata.var_names:\n{', '.join(dropped_genes[:10])}...")

    cnv_state_cnt(cluster_geno_expanded, save_path=osp.join(infercnv_out_path, "infer_cnv_state_freq_by_group.png"))

    # cell to grouping
    pat = "infercnv.observation_groupings.txt"
    infercnv_groupings_files = glob.glob(os.path.join(infercnv_out_path, pat))
    if len(infercnv_groupings_files) == 0:
        raise ValueError("No inferCNV groupings files found in the specified directory. Please check the path and pattern.")
    infercnv_groupings_df = pd.read_csv(infercnv_groupings_files[0], sep=" ")
    sub_infercnv_groupings_df = infercnv_groupings_df[['Dendrogram Group']]

    # to prevent "Dendrogram Group_x" from being created
    if 'Dendrogram Group' in adata.obs.columns:
        adata.obs = adata.obs.drop(columns=['Dendrogram Group'])

    # merge the infercnv groupings into the adata obs
    adata.obs = adata.obs.merge(sub_infercnv_groupings_df, left_index=True, right_index=True, how='left')
    valid_mask = adata.obs['Dendrogram Group'].notna()
    valid_groups = adata.obs.loc[valid_mask, 'Dendrogram Group']

    # add to obs
    # store the CNV states by cell in the adata obsm
    cnv_matrix = np.full((adata.n_obs, cluster_geno_expanded.shape[1]), np.nan)  # initialize with NaN
    cnv_matrix[valid_mask.values] = cluster_geno_expanded.loc[valid_groups].values
    adata.obsm["cnv_states"] = cnv_matrix

    # also store grouping_to_cnv_state in uns
    adata.uns["cnv_group_to_states"] = cluster_geno_expanded

    return adata


def prep_h5ad_for_infercnv(adata, infercnv_out_path, cell_type_col="cell_type", auto_infer=True, gene_ordering_file=None):
    # adata to matrix + sample annotation + gene ordering
    matrix_adata = h5ad_to_matrix(adata)
    size_gb = matrix_adata.memory_usage(deep=True).sum() / (1024**3)
    matrix_path = osp.abspath(osp.join(infercnv_out_path, 'singleCell.counts.matrix'))

    # dynamically adjust buffer size based on size of the matrix
    ## this does nothing but it at least looks cool
    if size_gb > 1:
        buffer_size = 1024 * 1024 * 64 
    else:
        buffer_size = -1
    
    with open(matrix_path, 'w', buffering=buffer_size) as f:
        matrix_adata.to_csv(f, sep="\t")

    sample_annotations = extract_infercnv_sampleann(adata, cell_type_col=cell_type_col, auto_infer=auto_infer)
    sample_annotations_path = osp.abspath(osp.join(infercnv_out_path, 'cellAnnotations.txt'))
    sample_annotations.to_csv(sample_annotations_path, sep="\t", index=False, header=False)

    # detect if gene ordering file is provided, if not, download the default one
    if gene_ordering_file is not None:
        if not osp.exists(gene_ordering_file):
            raise FileNotFoundError(f"Gene ordering file {gene_ordering_file} does not exist.")
        gene_order_path = auto_expand(gene_ordering_file)
    else:
        gene_order_path = osp.abspath(osp.join(infercnv_out_path, "gene_ordering_file.txt"))
        urllib.request.urlretrieve("https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt", gene_order_path)

    return matrix_path, sample_annotations_path, gene_order_path

def _r_bool(x):  # helper
    return "TRUE" if bool(x) else "FALSE"

## main function
def call_infercnv(matrix_path, sample_annotations_path, gene_order_path, infercnv_out_path, ref_group_names, 
                  adata=None, **kwargs):
    """
    Run inferCNV with the provided mtx/annotations/gene_ordering files.
    """

    cutoff = kwargs.get("cutoff")
    denoise = kwargs.get("denoise")
    HMM = kwargs.get("HMM")
    cluster_by_groups = kwargs.get("cluster_by_groups")
    num_threads = kwargs.get("num_threads")
    analysis_mode = kwargs.get("analysis_mode")
    tumor_subcluster_partition_method = kwargs.get("tumor_subcluster_partition_method")
    tumor_subcluster_pval = kwargs.get("tumor_subcluster_pval")
    sd_amplifier = kwargs.get("sd_amplifier")
    noise_logistic = kwargs.get("noise_logistic")
    BayesMaxPNormal = kwargs.get("BayesMaxPNormal")

    if not ref_group_names:
        raise ValueError("ref_group_names must be a non-empty list (auto-infer in CLI if needed).")

    r_script_path="/app/R/src/infercnv_helper.R"

    rscript_cmd = [
        "Rscript", r_script_path,
        "--matrix", matrix_path,
        "--annotations", sample_annotations_path,
        "--gene_order", gene_order_path,
        "--out_dir", infercnv_out_path,
        "--ref_group_names", *list(ref_group_names),
    ]    

    # add additional parameters if specified
    if cutoff is not None:
        rscript_cmd += ["--cutoff", str(cutoff)]
    if analysis_mode:
        rscript_cmd += ["--analysis_mode", str(analysis_mode)]
    if cluster_by_groups is not None:
        rscript_cmd += ["--cluster_by_groups", _r_bool(cluster_by_groups)]
    if denoise is not None:
        rscript_cmd += ["--denoise", _r_bool(denoise)]
    if tumor_subcluster_partition_method:
        rscript_cmd += ["--tumor_subcluster_partition_method", str(tumor_subcluster_partition_method)]
    if tumor_subcluster_pval is not None:
        rscript_cmd += ["--tumor_subcluster_pval", str(tumor_subcluster_pval)]
    if num_threads is not None:
        rscript_cmd += ["--num_threads", str(int(num_threads))]
    if sd_amplifier is not None:
        rscript_cmd += ["--sd_amplifier", str(sd_amplifier)]
    if noise_logistic is not None:
        rscript_cmd += ["--noise_logistic", _r_bool(noise_logistic)]
    if HMM is not None:
        rscript_cmd += ["--HMM", _r_bool(HMM)]
    if BayesMaxPNormal is not None:
        rscript_cmd += ["--BayesMaxPNormal", str(BayesMaxPNormal)]

    print("Running inferCNV with command:")
    print(" ".join(rscript_cmd))

    # do the call
    result = subprocess.run(
        rscript_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(f"inferCNV failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")

    if adata is not None:
        adata = store_infercnv_in_adata(adata, infercnv_out_path)
        return adata

    return None