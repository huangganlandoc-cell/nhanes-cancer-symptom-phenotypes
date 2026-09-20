"""
Download the raw NHANES and NCHS linked-mortality files used in this study.

Reproduces everything under data/raw_nhanes/ and data/mortality/ from the
public CDC endpoints. Run from the project root:

    python code/download_nhanes.py --out data

Notes
-----
* NHANES public data-file URLs follow the pattern
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{start_year}/DataFiles/{MODULE}{SUFFIX}.xpt
  An older pattern (.../Nchs/Nhanes/{cycle}/{MODULE}{SUFFIX}.XPT) now returns
  an HTML redirect page rather than an XPORT file; use the pattern below.
* The NCHS public-use linked mortality files are fixed-width .dat files with
  follow-up through 31 December 2019.
* Files are written as gzipped CSV to keep the archive small and readable from
  R, Python and Excel without a SAS/XPORT reader.
"""

import argparse
import io
import os

import pandas as pd
import requests

CYCLES = [(2005, "_D"), (2007, "_E"), (2009, "_F"), (2011, "_G"),
          (2013, "_H"), (2015, "_I"), (2017, "_J")]

MODULES = ["DEMO", "MCQ", "DPQ", "SLQ", "BMX", "SMQ", "BPQ", "DIQ", "PAQ", "ALQ", "HUQ"]

NHANES_URL = ("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/"
              "{year}/DataFiles/{module}{suffix}.xpt")

MORT_URL = ("https://ftp.cdc.gov/pub/HEALTH_STATISTICS/NCHS/datalinkage/"
            "linked_mortality/NHANES_{y}_{y2}_MORT_2019_PUBLIC.dat")

# Fixed-width layout of the NCHS public-use linked mortality file.
MORT_COLSPECS = [(0, 14), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21),
                 (21, 22), (22, 26), (26, 34), (34, 42), (42, 45), (45, 48)]
MORT_NAMES = ["SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "DIABETES",
              "HYPERTEN", "_r1", "_r2", "_r3", "_r4", "PERMTH_INT", "PERMTH_EXM"]


def fetch_xpt(year, module, suffix, timeout=240):
    url = NHANES_URL.format(year=year, module=module, suffix=suffix)
    r = requests.get(url, timeout=timeout)
    if r.status_code != 200:
        return None, url
    return pd.read_sas(io.BytesIO(r.content), format="xport"), url


def fetch_mortality(year, timeout=180):
    url = MORT_URL.format(y=year, y2=year + 1)
    r = requests.get(url, timeout=timeout)
    r.raise_for_status()
    return pd.read_fwf(io.StringIO(r.text), colspecs=MORT_COLSPECS, names=MORT_NAMES)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="data", help="output directory")
    args = ap.parse_args()

    raw_dir = os.path.join(args.out, "raw_nhanes")
    mort_dir = os.path.join(args.out, "mortality")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(mort_dir, exist_ok=True)

    manifest = []
    for year, suffix in CYCLES:
        for module in MODULES:
            df, url = fetch_xpt(year, module, suffix)
            if df is None:
                print(f"MISSING  {module}{suffix} {year}")
                continue
            name = f"{module}{suffix}_{year}-{year + 1}.csv.gz"
            df.to_csv(os.path.join(raw_dir, name), index=False, compression="gzip")
            manifest.append({"file": name, "module": module,
                             "cycle": f"{year}-{year + 1}",
                             "nhanes_suffix": suffix, "rows": len(df),
                             "cols": df.shape[1], "source_url": url})
            print(f"{name}  {df.shape}")

    pd.DataFrame(manifest).to_csv(os.path.join(raw_dir, "_manifest.csv"), index=False)

    frames = [fetch_mortality(year) for year, _ in CYCLES]
    mort = pd.concat(frames, ignore_index=True)
    mort["SEQN"] = pd.to_numeric(mort["SEQN"], errors="coerce")
    for col in ["ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "PERMTH_INT", "PERMTH_EXM"]:
        mort[col] = pd.to_numeric(mort[col], errors="coerce")
    out = os.path.join(mort_dir, "nchs_linked_mortality_2019.csv.gz")
    mort.to_csv(out, index=False, compression="gzip")
    print(f"mortality  {mort.shape}  ->  {out}")


if __name__ == "__main__":
    main()
