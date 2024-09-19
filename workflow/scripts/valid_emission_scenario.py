"""Extracts scenario results into a single file"""

import pandas as pd
import sys
from pathlib import Path
import os


def extract_emission_reduction_scenarios(results_dir: str, param: str) -> list[str]:
    """Extracts emission reduction scenarios for running script locally"""

    dirs = [f.path for f in os.scandir(results_dir) if f.is_dir()]
    emission_dirs = [d for d in dirs if d.endswith("_percent_reduction")]

    return [str(Path(x, "results", f"{param}.csv")) for x in emission_dirs]


def is_valid_files(parameter: str, csvs: list[str]) -> bool:
    """Checks that csv names match parameter name"""

    for csv in csvs:
        p = Path(csv)
        if p.stem != parameter:
            print(f"{p.stem} does not match parameter name {parameter}")
            return False
    return True


def extract_emission_reduction_percent(csv: str) -> str:
    """Extracts emission value percent from the input file path

    Assumes only one part of the file path ends with '_percent_reduction'
    """

    parts = Path(csv).parts
    for part in parts:
        if part.endswith("_percent_reduction"):
            return part.split("_percent_reduction")[0]


def is_result_indices_unique(df: pd.DataFrame) -> bool:
    """Checks that only 1 unique value exist in the indices"""

    for col in df.columns:
        if (col == "VALUE") or (col == "YEAR"):
            continue
        if len(df[col].unique()) != 1:
            print(f"Column {col} contains the unique values of {df[col].unique()}")
            return False
    return True


def get_valid_scenarios(csvs: list[str]) -> pd.DataFrame:
    """Determines valid scenarios in datastructure below

    Where the scenaro is the emission reducion targets.
    ie. 0 percent reducion, 10 percent reduction, ect...

    | reduction   | valid  |
    |-------------|--------|
    | 0           | True   |
    | 20          | True   |
    | 40          | True   |
    | ...         | ...    |
    | 100         | False  |
    """

    data = []

    for csv in csvs:

        valid = True

        reduction = extract_emission_reduction_percent(csv)

        df = pd.read_csv(csv)

        df = df[df.TECHNOLOGY.str.startswith("BST")]

        if not df.empty:
            if df.VALUE.sum() > 1:
                valid = False

        data.append([reduction, valid])

    return pd.DataFrame(data, columns=["reduction", "valid"])


if __name__ == "__main__":

    if "snakemake" in globals():

        in_csvs = snakemake.input.csvs
        out_csv = str(snakemake.output.csv)
        parameter = "CapitalInvestment"

    else:

        if len(sys.argv) != 3:
            msg = "Usage: python {} <parameter> <emission_red_dir> <save.csv>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            in_csv_dir = sys.argv[1]
            out_csv = sys.argv[2]

            parameter = "CapitalInvestment"

            in_csvs = extract_emission_reduction_scenarios(in_csv_dir, parameter)

    df = get_valid_scenarios(in_csvs)

    df.to_csv(out_csv, index=False)
