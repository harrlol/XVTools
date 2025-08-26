import os
import os.path as osp
import subprocess
import fnmatch
import urllib.request
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


def infercnv(adata, infercnv_out_path=None, cell_type_col="cell_type", auto_infer=True, gene_ordering_file=None, 
             cutoff=0.1, denoise=True, HMM=True, cluster_by_groups=True, num_threads=4):
    """
    Run inferCNV on the given AnnData object, effectively adds CNV as a metadata to each cell.
    """
    # check specification of tempdir
    if infercnv_out_path is None:
        dt = datetime.datetime.now()
        infercnv_out_path = osp.join(osp.expanduser("~"), 'infercnv_temp', dt.strftime('%Y-%m-%d_%H-%M-%S'))
    else:
        infercnv_out_path = auto_expand(infercnv_out_path)

    # make tempdir
    if not isinstance(infercnv_out_path, str):
        raise TypeError("infercnv_out_path must be a string path.")
    if not osp.exists(infercnv_out_path):
        os.makedirs(infercnv_out_path)

    # adata to matrix + sample annotation + gene ordering
    matrix_adata = h5ad_to_matrix(adata)
    size_gb = matrix_adata.memory_usage(deep=True).sum() / (1024**3)
    matrix_path = osp.abspath(osp.join(infercnv_out_path, 'infercnv_matrix.mtx'))

    # dynamically adjust buffer size based on size of the matrix
    ## this does nothing but it at least looks cool
    if size_gb > 1:
        buffer_size = 1024 * 1024 * 64 
    else:
        buffer_size = -1
    
    with open(matrix_path, 'w', buffering=buffer_size) as f:
        matrix_adata.to_csv(f, sep="\t")

    sample_annotations = extract_infercnv_sampleann(adata, cell_type_col=cell_type_col, auto_infer=auto_infer)
    sample_annotations_path = osp.abspath(osp.join(infercnv_out_path, 'cell_annotations.txt'))
    sample_annotations.to_csv(sample_annotations_path, sep="\t", index=False, header=False)

    # detect if gene ordering file is provided, if not, download the default one
    if gene_ordering_file is not None:
        if not osp.exists(gene_ordering_file):
            raise FileNotFoundError(f"Gene ordering file {gene_ordering_file} does not exist.")
        gene_order_path = auto_expand(gene_ordering_file)
    else:
        gene_order_path = osp.abspath(osp.join(infercnv_out_path, "gene_order.txt"))
        urllib.request.urlretrieve("https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt", gene_order_path)

    # note the script location is hard coded!
    r_script_path = osp.abspath(osp.join(osp.dirname(__file__), "../../R/src/infercnv_helper.R"))
    # subprocess call of R script
    rscript_cmd = [
        "Rscript",
        r_script_path,
        "--matrix", matrix_path,
        "--annotations", sample_annotations_path,
        "--gene_order", gene_order_path,
        "--out_dir", infercnv_out_path,
        "--cutoff", str(cutoff),
    ]

    # add additional parameters if specified
    if denoise:
        rscript_cmd += ["--denoise", "TRUE"]
    if HMM:
        rscript_cmd += ["--HMM", "TRUE"]
    if cluster_by_groups:
        rscript_cmd += ["--cluster_by_groups", "TRUE"]
    if num_threads:
        rscript_cmd += ["--num_threads", str(num_threads)]

    # do the call
    result = subprocess.run(
        rscript_cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(f"inferCNV failed:\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}")
    
    adata = store_infercnv_in_adata(adata, infercnv_out_path)

    return adata