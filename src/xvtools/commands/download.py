import typer
import os
import os.path as osp
from google.oauth2 import service_account
from google.cloud import storage
from pathlib import Path
from rich import print as rprint
import fnmatch
from tqdm import tqdm

app = typer.Typer()

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
    blobs = bucket.list_blobs(prefix=prefix)

    # prep local
    if not osp.exists(local_dir):
        os.makedirs(local_dir, exist_ok=True)
    if pattern is None:
        pattern = "*"
    if excluded_suffixes is None:
        excluded_suffixes = ()
    else:
        excluded_suffixes = tuple(excluded_suffixes)
    matched_blobs = []
    for blob in blobs:
        # keep only those that actually lie under the prefix
        if not blob.name.startswith(prefix):
            continue
        # make path relative to the prefix
        relative = blob.name[len(prefix):].lstrip("/")

        # include filter: apply pattern to the *relative* path
        if not fnmatch.fnmatch(relative, pattern or "*"):
            continue

        # exclude by suffix (tuple of endings)
        if excluded_suffixes and relative.endswith(tuple(excluded_suffixes)):
            continue

        matched_blobs.append(blob)


    ####### debug
    rprint(f"[dim]Found under prefix[/]: {prefix}")
    rprint(f"[dim]Matched blobs:[/] {len(matched_blobs)}")
    for b in matched_blobs[:5]:
        rprint(f"  • {b.name}")

    # scan and download
    for blob in tqdm(matched_blobs, desc="Downloading files", unit="file"):
        tqdm.write(f"→ {blob.name}")
        relative_path = blob.name[len(prefix):]
        local_file_path = os.path.join(local_dir, relative_path)
        os.makedirs(os.path.dirname(local_file_path), exist_ok=True)
        blob.download_to_filename(local_file_path)

    return


@app.command("terra")
def terra_download(
    secret: Path = typer.Option(..., exists=True, help="Path to Terra service account JSON."),
    gcs_dir: str = typer.Option(..., help="gs://bucket/prefix to download"),
    dest: Path = typer.Option(..., help="Local destination directory"),
    pattern: str = typer.Option("*", help="fnmatch include (default: all)"),
    exclude: list[str] = typer.Option([], help="Suffixes to exclude, e.g. .bam .bai"),
):
    """Download from Terra-backed GCS using a service account."""
    rprint(f"[bold]XVTools[/] downloading [cyan]{gcs_dir}[/] → [green]{dest}[/] ...")
    dest.mkdir(parents=True, exist_ok=True)
    download_terra_data(
        SECRET=str(secret),
        gcs_dir=gcs_dir,
        local_dir=str(dest),
        pattern=pattern,
        excluded_suffixes=exclude,
    )
    rprint("[green]Done.[/]")