import os, datetime, subprocess, time, re
from pathlib import Path
import typer
from rich import print as rprint
import click
from typing import List, Optional
import shlex, ast

# Local smoke test calls your current driver:
from InfercnvAzureAPI import cli as infercnv_cli

app = typer.Typer()

def _parse_groups(s: str) -> List[str]:
    s = s.strip()
    # try Python/JSON list
    try:
        lit = ast.literal_eval(s)
        if isinstance(lit, (list, tuple)):
            return [str(x).strip() for x in lit if str(x).strip()]
    except Exception:
        pass
    # comma-separated
    if "," in s:
        return [p.strip() for p in s.split(",") if p.strip()]
    # whitespace fallback
    return [p.strip() for p in s.split() if p.strip()]


def _run(cmd: list[str], env=None):
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if p.returncode != 0:
        raise RuntimeError(
            f"$ {' '.join(cmd)}\n--- STDOUT ---\n{p.stdout}\n--- STDERR ---\n{p.stderr}"
        )
    return p.stdout


def _try_sync_outputs(g_src: str, out: Path):
    """
    Best-effort incremental download of new/changed files.
    Swallows 404/empty-container errors while the job is still running.
    """
    try:
        _run([
            "azcopy", "sync", g_src, str(out),
            "--recursive", "--delete-destination=false"
        ])
        rprint("[dim]Partial sync: pulled latest outputs from Azure[/]")
    except RuntimeError as e:
        msg = str(e)
        # Ignore expected 'not ready yet' conditions
        if ("ResourceNotFound" in msg or "The specified resource does not exist" in msg
            or "cannot list files" in msg or "StatusCode=404" in msg):
            return
        rprint(f"[yellow]Partial sync skipped:[/] {msg.splitlines()[-1]}")


# Prefer AML shortlink if present; otherwise fall back to first URL-ish thing.
_URL_RE_PREF = re.compile(r"https://aka\.ms/amlt\?[\w=-]+")
_URL_RE_ANY  = re.compile(r"https?://\S+")

def _extract_portal_url(text: str) -> str | None:
    m = _URL_RE_PREF.search(text)
    if m:
        return m.group(0)
    m = _URL_RE_ANY.search(text)
    return m.group(0) if m else None

@app.command("infercnv-aml")
def infercnv_aml(
    data: Path = typer.Option(..., exists=True, help="Input (h5ad dir or mtx directory tree)"),
    out: Path  = typer.Option(..., help="Local output directory"),
    job_name: str = typer.Option(None, help="Override job name (default: timestamped)"),
    n_parallel: int = typer.Option(4),
    sku: str = typer.Option("8C15"),
    # keep a repeatable flag if you like (user can do --ref-group-name X --ref-group-name Y)
    ref_group_name: List[str] = typer.Option(None, "--ref-group-name", help="Repeatable: --ref-group-name Hepatocyte --ref-group-name T_NK"),
    # single string, comma/JSON/whitespace acceptable
    ref_group_names_str: Optional[str] = typer.Option(None, "--ref-group-names-str", help="Comma/JSON list, e.g. 'Hepatocyte, T_NK' or '[\"Hepatocyte\",\"T_NK\"]'"),
    # pass-through optional knobs (add as needed)
    cutoff: float = typer.Option(0.1, help="Cutoff for CNV detection (default: 0.1)"),
    analysis_mode: str = typer.Option("subclusters", help="Analysis mode (default: subclusters)"),
    cluster_by_groups: bool = typer.Option(False, help="Cluster by groups (default: FALSE)"),
    denoise: bool = typer.Option(True, help="Denoise the data (default: TRUE)"),
    tumor_subcluster_partition_method: str = typer.Option("random_trees", help="Tumor subcluster partition method (default: random_trees)"),
    tumor_subcluster_pval: float = typer.Option(0.05, help="Tumor subcluster p-value (default: 0.05)"),
    n_threads: int = typer.Option(2, help="Number of threads to use (default: 2)"),
    sd_amplifier: float = typer.Option(1.5, help="Standard deviation amplifier (default: 1.5)"),
    noise_logistic: bool = typer.Option(True, help="Use noise logistic (default: TRUE)"),
    HMM: bool = typer.Option(False, help="Use HMM (default: FALSE)"),
    BayesMaxPNormal: float = typer.Option(0.2, help="Bayes maximum P(Normal) (default: 0.2)"),
    cell_type_col: str = typer.Option("cell_type_infercnv", help="Column name for cell type (default: cell_type_infercnv)"),
):
    """
    Resolve AML job YAML from template, upload data (azcopy), submit with 'amlt',
    poll, then sync outputs back.
    """

    if not job_name:
        t = datetime.datetime.now().astimezone().strftime("%Y%m%d_%H%M%S")
        job_name = f"infercnv_azure_{t}"
    out.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["JOB_NAME"] = job_name
    env["N_PARALLEL"] = str(n_parallel)
    env["N_THREADS"]  = str(n_threads)
    env["SKU"]        = sku

    ref_groups: List[str] = []
    if ref_group_name:
        ref_groups.extend([x.strip() for x in ref_group_name if x and x.strip()])
    if ref_group_names_str:
        ref_groups.extend(_parse_groups(ref_group_names_str))
    if not ref_groups:
        raise typer.BadParameter("Provide reference groups via --ref-group-name (repeatable) or --ref-group-names-str.")
    env["REF_ARG"] = "--ref_group_names " + " ".join(shlex.quote(x) for x in ref_groups)

    print("Arguments:")
    print(env['JOB_NAME'])
    print(env['N_PARALLEL'])
    print(env['N_THREADS'])
    print(env['SKU'])
    print(env['REF_ARG'])

    # optional args collection
    opts = []
    if cutoff is not None:
        opts += ["--cutoff", str(cutoff)]
    if analysis_mode:
        opts += ["--analysis_mode", analysis_mode]
    if cluster_by_groups is not None:
        opts += ["--cluster_by_groups" if cluster_by_groups else "--no-cluster_by_groups"]
    if denoise is not None:
        opts += ["--denoise" if denoise else "--no-denoise"]
    if tumor_subcluster_partition_method:
        opts += ["--tumor_subcluster_partition_method", tumor_subcluster_partition_method]
    if tumor_subcluster_pval is not None:
        opts += ["--tumor_subcluster_pval", str(tumor_subcluster_pval)]
    if n_threads is not None:
        opts += ["--n_threads", str(n_threads)]
    if sd_amplifier is not None:
        opts += ["--sd_amplifier", str(sd_amplifier)]
    if noise_logistic is not None:
        opts += ["--noise_logistic" if noise_logistic else "--no-noise_logistic"]
    if HMM is not None:
        opts += ["--HMM" if HMM else "--no-HMM"]
    if BayesMaxPNormal is not None:
        opts += ["--BayesMaxPNormal", str(BayesMaxPNormal)]
    if cell_type_col:
        opts += ["--cell_type_col", cell_type_col]
    env["OPTS_ARG"] = " ".join(opts)

    print(env['OPTS_ARG'])

    tpl_path = Path(__file__).resolve().parents[1] / "templates" / "job_infercnv.yml.tmpl"
    resolved = out / "job_infercnv.resolved.yml"

    txt = tpl_path.read_text()
    for k, v in env.items():
        txt = txt.replace("${"+k+"}", v)
    resolved.write_text(txt)
    rprint(f"[green]Generated[/] {resolved}")

    # sas_data = os.environ.get("XVTOOLS_SAS_DATA")
    # if not sas_data:
    #     raise RuntimeError("Set env XVTOOLS_SAS_DATA to a SAS token for the data container.")
    g_dst = f"https://exvivocoldeastus.blob.core.windows.net/data/broad_infercnv_data/{job_name}?<SAS_TOKEN>"
    _run(["azcopy", "sync", str(data), g_dst, "--recursive", "--delete-destination=false"], env=env)
    rprint("[cyan]Uploaded input to Azure blob[/]")

    # Submit AMLT job (requires amlt configured)
    submit_out = _run(["amlt", "run", str(resolved), "--replace", job_name, "-d", job_name], env=env)
    rprint(f"[green]Submitted[/] {job_name}")
    g_src = f"https://exvivocoldeastus.blob.core.windows.net/projects/Projects/broad_infercnv_out/{job_name}?<SAS_TOKEN>"

    url = _extract_portal_url(submit_out) or _extract_portal_url(_run(["amlt", "status", job_name], env=env))
    if url:
        rprint(f"[bold]Portal:[/bold] [link={url}]{url}[/link]")

    _STATUS_RE = re.compile(r'\b(queued|scheduling|preparing|starting|running|completed|failed|canceled)\b', re.I)
    prev_status = None
    while True:
        s = _run(["amlt", "status", job_name], env=env)
        m = _STATUS_RE.search(s)
        status = m.group(1).lower() if m else None

        if status and status != prev_status:
            rprint(f"[dim]status:[/] {status}")
            prev_status = status

        if status in {"queued","scheduling","preparing","starting","running"}:
            _try_sync_outputs(g_src, out)
            time.sleep(30)
            continue

        if status in {"failed","canceled"}:
            rprint("[red]Job failed or canceled[/]")
            _try_sync_outputs(g_src, out)
            return
        
        if status in {"completed"} or "pass" in s:
            rprint("[green]Job completed[/]")
            break

    # Download outputs
    _run(["azcopy", "sync", g_src, str(out), "--recursive", "--delete-destination=false"])
    rprint("[cyan]Synced outputs to local[/]")

@app.command("infercnv-local")
def infercnv_local(
    data: Path = typer.Option(..., exists=True),
    out: Path  = typer.Option(...),
    n_parallel: int = typer.Option(4),
    n_threads: int  = typer.Option(4),
    cell_type_col: str = typer.Option("cell_type"),
    ref_group_names: list[str] = typer.Option(None)
):
    """
    Run your existing Python driver locally (no AML). Good for smoke tests.
    """
    out.mkdir(parents=True, exist_ok=True)
    # Your driver uses argparse; easiest is to mimic CLI argv:
    import sys
    argv = [
        "--data_folder", str(data),
        "--out_folder",  str(out),
        "--n_parallel",  str(n_parallel),
        "--n_threads",   str(n_threads),
        "--cell_type_col", cell_type_col,
    ]
    if ref_group_names:
        argv += ["--ref_group_names", *ref_group_names]

    saved = sys.argv
    try:
        sys.argv = ["infercnv_cli.py"] + argv
        infercnv_cli.main()  # calls your existing program
    finally:
        sys.argv = saved