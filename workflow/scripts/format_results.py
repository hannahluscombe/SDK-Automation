"""Module to format otoole result data for MinFin"""

import sys
import pandas as pd
from pathlib import Path 


def capital_investment_pwr(csv: str, round_to: int = 2) -> pd.DataFrame:
    capital_investment = pd.read_csv(csv) # r, t, y, value
    pwr_capex = capital_investment[capital_investment.TECHNOLOGY.str.startswith("PWR")].drop(columns=["REGION"])
    pwr_capex = pwr_capex.pivot(columns="TECHNOLOGY", index="YEAR", values="VALUE").fillna(0).round(round_to)
    return pwr_capex

def variable_cost_min(csv: str, round_to: int = 2) -> pd.DataFrame:
    variable_cost = pd.read_csv(csv) # r, t, y, value
    min_var_cost = variable_cost[
        variable_cost.TECHNOLOGY.str.startswith(("MIN", "IMP")) &
        variable_cost.TECHNOLOGY.str.endswith(("NGS", "COA", "OIL", "LFO", "HFO", "BIO"))
    ].drop(columns=["REGION"])
    return min_var_cost.pivot(columns="TECHNOLOGY", index="YEAR", values="VALUE").fillna(0).round(round_to)

def variable_cost_pwr(csv: str, round_to: int = 6) -> pd.DataFrame:
    variable_cost = pd.read_csv(csv) # r, t, y, value
    var_cost = variable_cost[variable_cost.TECHNOLOGY.str.startswith("PWR")].drop(columns=["REGION", "TECHNOLOGY"])
    return var_cost.groupby("YEAR").sum().fillna(0).round(round_to)

def fixed_cost_pwr(csv: str, round_to: int = 6) -> pd.DataFrame:
    fixed_cost = pd.read_csv(csv) # r, t, y, value
    fix_cost = fixed_cost[fixed_cost.TECHNOLOGY.str.startswith("PWR")].drop(columns=["REGION", "TECHNOLOGY"])
    return fix_cost.groupby("YEAR").sum().fillna(0).round(round_to)

def production_annual_pwr(csv: str, round_to: int = 2) -> pd.DataFrame:
    production = pd.read_csv(csv) # r, f, t, y, value
    prod = production[production.TECHNOLOGY.str.startswith("PWR")].drop(columns=["REGION", "TECHNOLOGY", "FUEL"])
    return prod.groupby("YEAR").sum().fillna(0).round(round_to)

# holds information on how to extract minfin results 
# {minfin_name: {variable: str, function: fn}}
MINFIN_DATA = {
    "PWR_CapitalInvestment": {"variable":"CapitalInvestment", "function":capital_investment_pwr},
    "MIN_AnnualVariableOperatingCost": {"variable":"AnnualVariableOperatingCost", "function":variable_cost_min},
    "PWR_AnnualVariableOperatingCost": {"variable":"AnnualVariableOperatingCost", "function":variable_cost_pwr},
    "PWR_AnnualFixedOperatingCost": {"variable":"AnnualFixedOperatingCost", "function":fixed_cost_pwr},
    "PWR_ProductionByTechnologyAnnual": {"variable":"ProductionByTechnologyAnnual", "function":production_annual_pwr},
}

# business logic
def format_data(in_dir: str, out_dir: str, minfin_file: str) -> pd.DataFrame:
    """Format dataframe
    
    Args:
        in_dir: path to otoole data dir
        out_dir: path to save dir of minfin data
        minfin_file: name of saved minfin file
    """
    try:
        csv = Path(in_dir, f"{MINFIN_DATA[minfin_file]['variable']}.csv")
    except KeyError:
        print(f"No configuration data for {minfin_file}")
        return 
    df = MINFIN_DATA[minfin_file]["function"](str(csv))
    df.to_csv(Path(out_dir, f"{minfin_file}.csv"))

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
        
    for minfin_file in minfin_files:
        format_data(in_dir, out_dir, Path(minfin_file).stem)