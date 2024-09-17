# SDK-Automation
This workflow facilitates the simultaneous execution of OSeMOSYS starter data kit (SDK) scenarios.

## Installation 

1. Clone the repository: 

```bash 
git clone https://github.com/hannahluscombe/SDK-Automation.git
```

2. Install [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html) 

3. Install the environment. Note, this may take a few minutes. 

```bash
conda env create -f workflow/envs/sdk.yaml
```

4. Activate the environment.

```bash
conda activate sdk
```

## Running options

Three reults configurations can be applied with this workflow. 

### Default SDK run 

This will run and solve all user provided SDK scenarios, and present results in otoole CSV format. 

#### Configuration

The following configuration options are available when running the default workflow. 

```yaml 
# path to sdk scenarios
run: "custom" # (zenodo|custom)
scenario:
  zenodo: config/SDK_zenodo.csv
  custom: config/SDK_custom.csv

# solver 
solver: cplex # (cbc|cplex|gurobi)

# osemosys setup 
model: "resources/osemosys_pp.txt"
otoole: "resources/otoole.yaml"
```

| Configuration     | Description                                               |
|-------------------|-----------------------------------------------------------|
| `run`             | Run a custom SDK, or run one from provided Zenodo deposit |
| `scenario.zenodo` | Path to scenarios from Zenodo deposits                    |
| `scneario.custom` | Path to scenarios from custom `*.xlsm` files              |
| `solver`          | Solver for the models                                     |
| `model`           | Path to OSeMOSYS model file                               |
| `otoole`          | Path to otoole configuration file                         |

The Zenodo SDK configuration file must be formatted as follows:

| country_scenario | data                                       |
|------------------|--------------------------------------------|
| scenario_name_1  | URL of Zenodo file containing the `*.xlsm` |
| scenario_name_2  | URL of Zenodo file containing the `*.xlsm` |
| ...              | ...                                        |
| scenario_name_x  | URL of Zenodo file containing the `*.xlsm` |

The Custom SDK configuration file must be formatted as follows:

| country_scenario | data                                                     |
|------------------|----------------------------------------------------------|
| scenario_name_1  | Name of custom `.xlsm` in `config/sdk_custom/` directory |
| scenario_name_2  | Name of custom `.xlsm` in `config/sdk_custom/` directory |
| ...              | ...                                                      |
| scenario_name_x  | Name of custom `.xlsm` in `config/sdk_custom/` directory |

#### Execution

This is the default running option and can be executed with the following command, replace `-j6` with the number of cores.

```bash 
snakemake -j6
```

### MinFin run

This will run do the same operations as the default scenario. Additionaly, data will be prepared to be imported into MinFin. 

#### Configuration

All configuration options from default runs are available. Additionaly, the following configuration options are available when running in MinFin mode. 

```yaml
# minfin data to extract
# NOTE: these names must match names in 'format_data.py' module!
minfin_files:
  - PWR_CapitalInvestment
  - MIN_AnnualVariableOperatingCost
  - PWR_AnnualVariableOperatingCost
  - PWR_AnnualFixedOperatingCost
  - PWR_ProductionByTechnologyAnnual
```

| Configuration  | Description                                                                    |
|----------------|--------------------------------------------------------------------------------|
| `minfin_files` | Functions in the `minfin_results.py` script to execute and create results for  |

#### Execution

Generating MinFin results can be done with the following command, replacing `-j6` with the number of cores.

```bash 
snakemake minfin -j6
```

### Emission Abatement Curves 

This will run do the same operations as the default scenario. Additionaly, scenario runs will be run at user defined emission constraint intervals. A Pareto Abatement Cost Curve figure will be generated. 


#### Configuration 

All configuration options from default runs are available. Additionaly, the following configuration options are available when running for emission abatement curves. 

```yaml 
# emission abatement config
emission_curves:
  config: "config/emission_curves.csv"
  base_year: 2050
  reduction_method: "linear" # (linear|exp)
  include_backstop: False
```

| Configuration | Description |
|---|---|
| `emission_curves.config` | Path to configuration file describing the emission abatement scenarios |
| `emission_curves.base_year` | Year to which apply the emission abatement curves against. For example, a 50% emission reduction objective against `2050` |
| `emission_curves.reduction_method` | Method to reduce emissions. Only "linear" is currently supported |
| `emission_curves.include_backstop` | Include scenario runs where backstop technologies are invested. If `False`, scenarios that include backstop are not reported in the Pareto Abatement Cost Curve figure - the scenarios will still be run and solved, just not included in the figure.  |

The Emission Abatement Scenario configuration file must be formatted as follows:

| emission_scenario | emission_reduction |
|-------------------|--------------------|
| 0                 | 0                  |
| 1                 | 20                 |
| 2                 | 40                 |
| ...               |                    |
| n                 | xx                 |

Where `emission_scenario` is a increasing integer. And `emission_reduction` is an integer between `0` and `100`. **You MUST provide the `scenario` `0` and `emission_reduction` `0` config option**  

#### Execution

Generating Emission Abatement Curve results can be done with the following command, replacing `-j6` with the number of cores.

```bash 
snakemake emission_curves -j6
```