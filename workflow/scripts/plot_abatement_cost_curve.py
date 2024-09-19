"""Plots Emission Abatement Cost Curve"""

import pandas as pd
import matplotlib.pyplot as plt
import sys


def get_valid_scenarios(csv: str) -> list[int]:
    df = pd.read_csv(csv, dtype={"reduction": int, "valid": bool})
    return df[df.valid == True].reduction.to_list()


def filter_on_valid_scenarios(
    df: pd.DataFrame, valid_scenarios: list[int | str]
) -> pd.DataFrame:
    """Filters results on valid scenarios"""

    valid_scenarios.sort()
    valid_scenarios = [str(x) for x in valid_scenarios]

    scenarios = df.columns

    return df[[x for x in scenarios if x in valid_scenarios]]


def is_valid_year(year: int, df: pd.DataFrame) -> bool:
    """Checks for valid year"""

    if int(year) in df.index:
        return True
    else:
        print(f"{year} not in index values of {df.index}")
        return False


def are_valid_dataframes(df1: pd.DataFrame, df2: pd.DataFrame) -> bool:
    """Checks for valid index, columns between dataframes"""

    if not df1.index.equals(df2.index):
        print("Indices between df1 and df2 are different\n")
        print(f"df1 index is {df1.index}\n")
        print(f"df2 index is {df2.index}\n")
        return False

    if not df1.columns.equals(df2.columns):
        print("Columns between df1 and df2 are different\n")
        print(f"df1 columns is {df1.columns}\n")
        print(f"df2 columns is {df2.columns}\n")
        return False

    return True


def plot_abatement_curve(
    year: int, emissions: pd.DataFrame, costs: pd.DataFrame
) -> tuple:

    year = int(year)

    year_emissions = emissions.loc[year]
    year_costs = costs.loc[year]

    df = pd.DataFrame({"Emissions": year_emissions, "Costs": year_costs}).set_index(
        "Emissions"
    )

    fig, ax = plt.subplots(1, 1)

    df.plot(kind="area", ax=ax)

    ax.set(xlabel="Emissions (MMT)", ylabel="Annual Discounted Cost (M$)")
    fig.suptitle(f"Pareto Abatement Cost Curve for {year}")

    return fig, ax


if __name__ == "__main__":

    if "snakemake" in globals():
        year = snakemake.params.year
        valid_scenarios_csv = snakemake.input.valid
        emissions_csv = snakemake.input.annual_emissions
        costs_csv = snakemake.input.total_cost
        plot_png = str(snakemake.output.plot)

    else:

        if len(sys.argv) != 6:
            msg = "Usage: python {} <year> <valid_scenarios.csv> <annual_emissions.csv> <total_costs.csv> <plot.png>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            year = sys.argv[1]
            valid_scenarios_csv = sys.argv[2]
            emissions_csv = sys.argv[3]
            costs_csv = sys.argv[4]
            plot_png = sys.argv[5]

    emissions = pd.read_csv(emissions_csv, index_col=[0])
    costs = pd.read_csv(costs_csv, index_col=[0])

    if valid_scenarios_csv:
        valid_scenarios = get_valid_scenarios(valid_scenarios_csv)
        emissions = filter_on_valid_scenarios(emissions, valid_scenarios)
        costs = filter_on_valid_scenarios(costs, valid_scenarios)

    assert is_valid_year(year, emissions)
    assert is_valid_year(year, costs)
    assert are_valid_dataframes(emissions, costs)

    fig, ax = plot_abatement_curve(year, emissions, costs)

    fig.savefig(plot_png)
