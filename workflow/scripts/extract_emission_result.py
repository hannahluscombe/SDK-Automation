"""Extracts scenario results into a single file"""

import pandas as pd
import sys
from pathlib import Path


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
        if col == "VALUE":
            continue
        if len(df[col].unique()) != 1:
            print(f"Column {col} contains the unique values of {df.unique()}")
            return False
    return True


def extract_data(csvs: list[str]) -> pd.DataFrame:
    """Extracts result data into format

    Where the index is the year and the column headers are the emission reducion
    targets. ie. 0 percent reducion, 10 percent reduction, ect...

    |      | 0   | 10  | 20  | ... | 100 |
    |------|-----|-----|-----|-----|-----|
    | 2020 | ### | ### | ### | ... | ### |
    | 2021 | ### | ### | ### | ... | ### |
    | 2022 | ### | ### | ### | ... | ### |
    | ...  | ... | ... | ... | ... | ... |
    | 2050 | ### | ### | ### | ### | ### |
    """

    dfs = []

    for csv in csvs:

        df = pd.read_csv(csv)
        assert is_result_indices_unique(df)

        reduction = extract_emission_reduction_percent(csv)

        df = (
            df[["YEAR", "VALUE"]].set_index("YEAR").rename(columns={"VALUE": reduction})
        )

        dfs.append(df)

    return pd.concat(dfs, axis=1)


if __name__ == "__main__":

    if "snakemake" in globals():

        parameter = snakemake.params.result
        in_csvs = snakemake.input
        out_csv = snakemake.output

    else:

        if len(sys.argv) != 4:
            msg = "Usage: python {} <parameter> <[parameter.csv]> <summary.csv>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            parameter = sys.argv[1]
            in_csvs = sys.argv[2]
            out_csv = sys.argv[3]

    assert is_valid_files(parameter, in_csvs)

    df = extract_data(in_csvs)

    df.to_csv(out_csv, index=True)
