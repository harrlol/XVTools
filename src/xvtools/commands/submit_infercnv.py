import os, datetime, subprocess, time, re
from pathlib import Path
import typer
from rich import print as rprint

# Local smoke test calls your current driver:
from InfercnvAzureAPI import cli as infercnv_cli

app = typer.Typer()

def _run(cmd: list[str], env=None):
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if p.returncode != 0:
        raise RuntimeError(
            f"$ {' '.join(cmd)}\n--- STDOUT ---\n{p.stdout}\n--- STDERR ---\n{p.stderr}"
        )
    return p.stdout

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
    n_threads: int  = typer.Option(2),
    sku: str = typer.Option("8C15"),
    ref_group_names: list[str] = typer.Option(None, help="Space-separated, e.g. 'normal'")
    # pass-through optional knobs (add as needed)
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
    env["REF_ARG"]    = f"--ref_group_names {ref_group_names}" if ref_group_names else ""
    env["OPTS_ARG"]   = ""  # extend if you expose more switches

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
            time.sleep(30)
            continue
        if status in {"failed","canceled"}:
            rprint("[red]Job failed or canceled[/]")
            return
        if status in {"completed"} or "pass" in s:
            rprint("[green]Job completed[/]")
            break

    # Download outputs
    # sas_out = os.environ.get("XVTOOLS_SAS_PROJECTS")
    # if not sas_out:
    #     raise RuntimeError("Set env XVTOOLS_SAS_PROJECTS to a SAS token for the projects container.")
    g_src = f"https://exvivocoldeastus.blob.core.windows.net/projects/Projects/broad_infercnv_out/{job_name}?<SAS_TOKEN>"
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