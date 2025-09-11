import os
import argparse
import scanpy as sc
from pathlib import Path
import urllib.request
from ._util import *
from multiprocessing import Pool

def _worker(f_path, out_folder, cell_type_col, worker_msg, **kwargs):
    f_path = Path(f_path)
    sample_name = f_path.stem
    infercnv_out_path = Path(out_folder) / sample_name
    infercnv_out_path.mkdir(parents=True, exist_ok=True)
    print(f"[Azure] Created directory: {infercnv_out_path}")

    print(f"[Azure] Reading {sample_name} at {f_path} ...")
    if worker_msg == "mtx":
        adata = None
        matrix_path = str(f_path / "singleCell.counts.matrix")
        sample_annotations_path = str(f_path / "cellAnnotations.txt")
        gene_order_path = str(f_path / "gene_ordering_file.txt")

        # patch 8/29, allow gene_ordering_file.txt to be downloaded on the spot under mtx
        if not Path(gene_order_path).exists():
            print("[Azure] Fetching gene ordering file...")
            urllib.request.urlretrieve("https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt", gene_order_path)
            if Path(gene_order_path).exists():
                print(f"[Azure] Gene ordering file successfully saved at {gene_order_path}")
            else:
                raise FileNotFoundError(
                    f"[Azure][ERROR] Download attempted but gene ordering file not found at {gene_order_path}"
                )

        missing = [p for p in (matrix_path, sample_annotations_path, gene_order_path) if not Path(p).exists()]
        if missing:
            raise FileNotFoundError(f"MTX mode expects files next to {f_path} but missing: {missing}")
    elif worker_msg == "h5ad":
        adata = sc.read_h5ad(f_path)
        matrix_path, sample_annotations_path, gene_order_path = prep_h5ad_for_infercnv(
            adata, infercnv_out_path=str(infercnv_out_path), auto_infer=False, cell_type_col=cell_type_col
        )
    else:
        raise ValueError("worker_msg must be either 'mtx' or 'h5ad'.")
    
    # 9/11 patch: require ref_group_names
    if kwargs.get("ref_group_names", None) is None:
        raise ValueError("[Azure][ERROR] 'ref_group_name' must be provided.")
    
    # 9/11 patch: check ref_group_names is surjective
    ref_names = kwargs.get("ref_group_names", [])
    ann = pd.read_csv(sample_annotations_path, sep="\t", header=0, index_col=0)
    ann_groups = set(ann[cell_type_col].unique())
    missing_refs = [r for r in ref_names if r not in ann_groups]
    candidate_malig = [g for g in ann_groups if g not in ref_names]
    if missing_refs:
        raise ValueError(f"[Azure][ERROR] ref_group_names {missing_refs} not found in sample annotations groups {ann_groups}.")
    if len(candidate_malig) == 0:
        raise ValueError(f"[Azure][ERROR] No malignant group found. All groups {ann_groups} are listed as references.")
    else:
        print(f"[Azure] The following cell type label will be considered malignant: {candidate_malig} (all others are references)")
    
    # pass paths to infercnv
    print(f"[Azure] Starting infercnv for {sample_name} ...")
    adata = call_infercnv(matrix_path, sample_annotations_path, gene_order_path, str(infercnv_out_path), 
                          adata=adata, **kwargs)

    # adata is not none only in h5ad mode
    if adata is not None:
        out_h5ad = infercnv_out_path / f"{sample_name}.infercnv.h5ad"
        print(f"Writing output h5ad to {out_h5ad} ...")
        adata.write_h5ad(out_h5ad, compression="gzip")

    return sample_name

def _worker_wrapper(f_path, out_folder, cell_type_col, worker_msg, infercnv_kwargs):
    return _worker(f_path, out_folder, cell_type_col, worker_msg, **infercnv_kwargs)

def main(args=None):
    print("[Azure] CLI Arguments", args)
    
    os.makedirs(args.out_folder, exist_ok=True)

    # detect if we are in mtx mode or h5ad mode
    data_p = Path(args.data_folder)

    h5ads = sorted([p for p in data_p.iterdir() if p.is_file() and p.suffix.lower() == ".h5ad"])
    mtx_dirs = sorted({p.parent for p in data_p.rglob("*.matrix")})

    if h5ads:
        print("[Azure] Detected .h5ad files, proceeding in h5ad mode.")
        worker_msg = "h5ad"
        f_list = h5ads
    elif mtx_dirs:
        print("[Azure] Detected .mtx files, proceeding in mtx mode.")
        worker_msg = "mtx"
        f_list = mtx_dirs
    else:
        raise FileNotFoundError(
        "[Azure] No .h5ad or .mtx files found in input folder."
    )

    infercnv_kwargs = {
        "cutoff": args.cutoff,
        "denoise": args.denoise,
        "HMM": args.HMM,
        "ref_group_names": args.ref_group_names,
        "num_threads": args.n_threads,
        "analysis_mode": args.analysis_mode,
        "cluster_by_groups": args.cluster_by_groups,
        "tumor_subcluster_partition_method": args.tumor_subcluster_partition_method,
        "tumor_subcluster_pval": args.tumor_subcluster_pval,
        "sd_amplifier": args.sd_amplifier,
        "noise_logistic": args.noise_logistic,
        "BayesMaxPNormal": args.BayesMaxPNormal,
    }

    tasks = [
        (str(f), args.out_folder, args.cell_type_col, worker_msg, infercnv_kwargs)
        for f in f_list
    ]

    with Pool(processes=args.n_parallel) as pool:
        for name in pool.starmap(_worker_wrapper, tasks):
            print(f"✔ finished: {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="InferCNV Azure API")

    # path to local h5ad files
    parser.add_argument("--data_folder", type=str, required=True, help="Path to local data files folder")
    parser.add_argument("--out_folder", type=str, required=True, help="Path to local output folder")
    parser.add_argument("--ref_group_names", nargs="*", required=True, help="Reference group names for infercnv")
    parser.add_argument("--n_parallel", type=int, required=True, help="Number of samples to process in parallel")
    parser.add_argument("--n_threads", type=int, default=2, help="Number of threads allowed per sample")
    parser.add_argument("--cell_type_col", type=str, default="cell_type_infercnv", help="Column name for cell type")
    # infercnv parameters
    parser.add_argument("--cutoff", type=float, default=0.1, help="Cutoff for infercnv (default: 0.1)")
    parser.add_argument("--denoise", action="store_true", default=True)
    parser.add_argument("--no-denoise", dest="denoise", action="store_false")
    parser.add_argument("--HMM", action="store_true", default=True)
    parser.add_argument("--no-HMM", dest="HMM", action="store_false")
    parser.add_argument("--analysis_mode", type=str, default=None)
    parser.add_argument("--cluster_by_groups", action="store_true", default=None)
    parser.add_argument("--no-cluster_by_groups", dest="cluster_by_groups", action="store_false")
    parser.add_argument("--tumor_subcluster_partition_method", type=str, default=None)
    parser.add_argument("--tumor_subcluster_pval", type=float, default=None)
    parser.add_argument("--sd_amplifier", type=float, default=None)
    parser.add_argument("--noise_logistic", action="store_true", default=None)
    parser.add_argument("--no-noise_logistic", dest="noise_logistic", action="store_false")
    parser.add_argument("--BayesMaxPNormal", type=float, default=None)
    args = parser.parse_args()

    main(args)
