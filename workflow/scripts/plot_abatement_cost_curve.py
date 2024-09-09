"""Plots Emission Abatement Cost Curve"""

import pandas as pd
import matplotlib.pyplot as plt


def is_valid_year(year: int, df: pd.DataFrame) -> bool:
    """Checks for valid year"""

    if year in df.index:
        return True
    else:
        print(f"{year} not in index values of {df.index}")
        return False


def are_valid_dataframes(df1: pd.DataFrame, df2: pd.DataFrame) -> bool:
    """Checks for valid index, columns between dataframes"""

    if df1.index != df2.index:
        print("Indices between df1 and df2 are different\n")
        print(f"df1 index is {df1.index}\n")
        print(f"df2 index is {df2.index}\n")
        return False

    if df1.columns != df2.columns:
        print("Columns between df1 and df2 are different\n")
        print(f"df1 columns is {df1.columns}\n")
        print(f"df2 columns is {df2.columns}\n")
        return False

    return True


def plot_abatement_curve(
    year: int, emissions: pd.DataFrame, costs: pd.DataFrame
) -> tuple:

    year_emissions = emissions.loc[year]
    year_costs = costs.loc[year]

    fig, ax = plt.subplots(1, 1)

    return fig, ax


if __name__ == "__main__":

    if "snakemake" in globals():
        year = snakemake.params.year
        emissions_csv = snakemake.input.annual_emissions
        costs_csv = snakemake.input.total_cost
        plot_png = snakemake.output.plot

    else:

        if len(sys.argv) != 5:
            msg = "Usage: python {} <year> <annual_emissions.csv> <total_costs.csv> <plot.png>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            year = sys.argv[1]
            emissions_csv = sys.argv[2]
            costs_csv = sys.argv[3]
            plot_png = sys.argv[3]

    emissions = pd.read_csv(emissions_csv, index_col=[0])
    costs = pd.read_csv(costs_csv, index_col=[0])

    assert is_valid_year(year, emissions)
    assert is_valid_year(year, costs)
    assert are_valid_dataframes(emissions, costs)

    fig, ax = plot_abatement_curve(year, emissions, costs)

    fig.savefig(plot_png)
