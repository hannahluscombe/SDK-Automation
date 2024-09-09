"""Creates emission curves"""

from typing import Optional
import pandas as pd
import numpy as np
import sys


def create_emission_curves(
    df: pd.DataFrame,
    reduction: float,
    base_year: int,
    region: Optional[str] = "RE1",
    emission: Optional[str] = "EMIC02",
    method: Optional[str] = "linear",
) -> pd.DataFrame:
    """Creates a dataframe of emission reduction curves

    Arguments:
        df: pd.DataFrame
            AnnualEmissions result dataframe
        reductions:float
            Emission reduction target as a percentage
        base_year: int
            referece year to take for '% reductions from x year'
        region: str
            Region to index over
        emission: str
            Emission to index over
        method: str
            'linear' or 'exp'

    Returns:
        Dataframe with years on the index, and reductions on the columns
    """

    assert df.index.names == ["REGION", "EMISSION", "YEAR"]

    start_year = df.index.get_level_values("YEAR").min()
    start_emissions = df.loc[(region, emission, start_year), "VALUE"]

    ref_emissions = df.loc[(region, emission, base_year), "VALUE"]

    curve = df.copy()
    curve["VALUE"] = np.nan
    curve.at[(region, emission, int(start_year)), "VALUE"] = start_emissions.round(4)
    curve.at[(region, emission, base_year), "VALUE"] = (
        ref_emissions * (100 - reduction) / 100
    )

    return curve.interpolate()


if __name__ == "__main__":

    if "snakemake" in globals():
        annual_emissions = snakemake.input.emissions
        reduction = float(snakemake.wildcards.emission_reduction)
        base_year = int(snakemake.params.base_year)
        method = snakemake.params.reduction_method
        region = snakemake.params.region
        emission = snakemake.params.emission
        output = snakemake.output
    else:
        if len(sys.argv) != 6:
            msg = "Usage: python {} <AnnualEmissions.csv> <reduction> <base_year> <reduction_method> <save_name>"
            print(msg.format(sys.argv[0]))
            sys.exit(1)
        else:
            annual_emissions = sys.argv[1]
            reduction = float(sys.argv[2])
            base_year = int(sys.argv[3])
            method = sys.argv[4]
            region = "RE1"
            emission = "EMIC02"
            output = sys.argv[5]

    emissions = pd.read_csv(annual_emissions, index_col=[0, 1, 2])
    
    if method not in ("linear", "exp"):
        print(f"{method} not a supported reduction method. Setting to 'linear'")
        method = "linear"

    curve = create_emission_curves(
        emissions, reduction, base_year, region, emission, method
    )

    curve.to_csv(output)
