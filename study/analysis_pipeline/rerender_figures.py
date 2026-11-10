"""Recreate the five figures from the package's fixed plot tables in a fresh directory."""
import argparse
from pathlib import Path
import shutil
import subprocess

here = Path(__file__).resolve().parent
root = here.parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("output",type=Path,help="New, previously nonexistent working directory")
args = parser.parse_args()
target = args.output.expanduser().resolve()
target.mkdir(parents=True,exist_ok=False)
shutil.copytree(here/"polished"/"code",target/"code")
shutil.copytree(root/"source_data",target/"source_data")
subprocess.run(["Rscript","code/render.R"],cwd=target,check=True)
print(target/"figures")

