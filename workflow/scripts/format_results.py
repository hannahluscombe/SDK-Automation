"""Module to format otoole result data for MinFin"""

import sys
import pandas as pd
from pathlib import Path 

def format_data(df: pd.DataFrame) -> pd.DataFrame:
    """Format dataframe"""
    return df

if __name__ == "__main__":
    
    if "snakemake" in globals():
        in_dir = snakemake.params.in_dir
        out_dir = snakemake.params.out_dir
        minfin_files = snakemake.params.minfin_files
    else:
        if len(sys.argv) != 4:
            msg = "Usage: python {} <in_dir> <out_dir> <file_names>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            in_dir = sys.argv[1]
            out_dir = sys.argv[2]
            minfin_files = sys.argv[3]
    
    if not Path(out_dir).exists():
        Path(out_dir).mkdir()
        
    for f in Path(in_dir).iterdir():
        if f.is_file():
            if f.name in minfin_files:
                df = format_data(pd.read_csv(f))
                df.to_csv(Path(out_dir, f.name))