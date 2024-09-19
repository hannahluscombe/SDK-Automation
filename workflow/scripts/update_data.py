"""Updates SDK data"""

import pandas as pd
import sys
from otoole import read, write
import yaml
from pathlib import Path


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
        data_input_type = snakemake.params.data_input_type
        new_data = snakemake.input.csv
        data_in = snakemake.input.data_in
        save_dir = snakemake.params.save_dir
    else:
        if len(sys.argv) != 7:
            msg = "Usage: python {} <otoole_config.yaml> <parameter> <new_data.csv> <in_datafile.txt> <out_datafile.txt>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            config = sys.argv[1]
            parameter = sys.argv[2]
            data_input_type = sys.argv[3]
            new_data = sys.argv[4]
            data_in = sys.argv[5]
            save_dir = sys.argv[6]

    if data_input_type == "csv":
        assert isinstance(data_in, list)
        data_in_dir = str(Path(data_in[0]).parent)
        data, defaults = read(config, "csv", data_in_dir)
    elif data_input_type == "datafile":
        data, defaults = read(config, "datafile", data_in)
    elif data_input_type == "excel":
        data, defaults = read(config, "excel", data_in)
    else:
        raise NotImplementedError

    with open(config) as f:
        parsed_config = yaml.safe_load(f)

    df = pd.read_csv(new_data)
    df = df.set_index([x for x in df.columns if x != "VALUE"])

    assert is_correct_index(df, parameter, parsed_config)

    data = update_data(data, parameter, df)

    write(config, "csv", save_dir, data, defaults)
