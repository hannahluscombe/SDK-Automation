"""Updates SDK data"""

import pandas as pd
import sys
from otoole import read, write
import yaml


def is_correct_index(
    df: pd.DataFrame, parameter: str, config: dict[str, dict[str, str | int | float]]
) -> bool:
    """Check the index of the updating data"""

    valid_params = [x for x in config if config[x]["type"] == "param"]

    if parameter not in valid_params:
        raise KeyError(f"{parameter} is not a valid parameter")

    expected = config[parameter]["indices"]
    actual = df.index.names

    if sorted(expected) == sorted(actual):
        return True
    else:
        return False


def update_data(
    data: dict[str, pd.DataFrame], parameter: str, new: pd.DataFrame
) -> dict[str, pd.DataFrame]:
    """Corrects Input Data"""
    data[parameter].update(new)
    return data


if __name__ == "__main__":
    if "snakemake" in globals():
        config = snakemake.params.config
        parameter = snakemake.params.parameter
        new_data = snakemake.input.csv
        txt_in = snakemake.input.txt
        save_dir = snakemake.params.save_dir
    else:
        if len(sys.argv) != 6:
            msg = "Usage: python {} <otoole_config.yaml> <parameter> <new_data.csv> <in_datafile.txt> <out_datafile.txt>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            config = sys.argv[1]
            parameter = sys.argv[2]
            new_data = sys.argv[3]
            txt_in = sys.argv[4]
            save_dir = sys.argv[5]

    data, defaults = read(config, "datafile", txt_in)

    with open(config) as f:
        parsed_config = yaml.safe_load(f)

    df = pd.read_csv(new_data)
    df = df.set_index([x for x in df.columns if x != "VALUE"])

    assert is_correct_index(df, parameter, parsed_config)

    data = update_data(data, parameter, df)

    write(config, "csv", save_dir, data, defaults)
