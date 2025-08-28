import os
import argparse
import matplotlib.pyplot as plt
import datetime
import numpy as np
import os.path as osp
import subprocess
import fnmatch
import scanpy as sc
import urllib.request
import glob
from tqdm import tqdm
from google.oauth2 import service_account
from google.cloud import storage
from pathlib import Path
from _util import *
from multiprocessing import Pool

def do_infercnv(adata, cell_type_col, out_path):
    adata = infercnv(
            adata,
            cell_type_col=cell_type_col,
            auto_infer=False,
            infercnv_out_path=out_path
        )
    return adata

def _worker(h5ad_path, out_folder, cell_type_col):
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    h5ad_path = Path(h5ad_path)
    sample_name = h5ad_path.stem

    out_dir = Path(out_folder) / sample_name
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = sc.read_h5ad(h5ad_path)
    adata = do_infercnv(adata, cell_type_col=cell_type_col, out_path=str(out_dir))

    out_h5ad = out_dir / f"{sample_name}.infercnv.h5ad"
    adata.write_h5ad(out_h5ad, compression="gzip")

    return

def main(args=None):
    print("CLI Arguments", args)
    
    os.makedirs(args.out_folder, exist_ok=True)
    h5ad_paths = [p for p in Path(args.h5ad_folder).iterdir() if p.is_file() and p.suffix.lower()==".h5ad"]
    if not h5ad_paths:
        print("No .h5ad files found.")
        return
    tasks = [(str(p), args.out_folder, args.cell_type_col) for p in h5ad_paths]
    with Pool(processes=args.n_parallel) as pool:
        for name in pool.starmap(_worker, tasks):
            print(f"✔ finished: {name}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="InferCNV Azure API")

    # path to local h5ad files
    parser.add_argument("--h5ad_folder", type=str, required=True, help="Path to local h5ad files folder")
    parser.add_argument("--out_folder", type=str, required=True, help="Path to local output folder")
    parser.add_argument("--n_parallel", type=int, required=True, help="Number of samples to process in parallel")
    parser.add_argument("--cell_type_col", type=str, default="cell_type_infercnv", help="Column name for cell type")
    args = parser.parse_args()

    main(args)