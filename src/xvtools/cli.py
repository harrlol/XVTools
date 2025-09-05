import typer
from xvtools.commands.download import app as download_app
from xvtools.commands.submit_infercnv import app as submit_app

app = typer.Typer(help="XVTools: Terra↔Azure helpers for inferCNV & beyond.")
app.add_typer(download_app, name="download", help="Data download helpers (e.g., Terra/GCS).")
app.add_typer(submit_app,   name="submit",   help="Submit or run inferCNV jobs.")

def main():
    app()

if __name__ == "__main__":
    main()