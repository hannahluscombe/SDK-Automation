# SDK-Automation
This script facilitates the simultaneous execution of OSeMOSYS starter data kit scenarios across multiple countries in a single run. The input file is an Excel spreadsheet containing scenario names and corresponding country names, along with the URL links to download the OSeMOSYS sand files.

## Install and Run 
1. Clone the repository: 

```bash 
git clone https://github.com/hannahluscombe/SDK-Automation.git
```

2. Install [anaconda](https://www.anaconda.com/) or [miniconda](https://docs.anaconda.com/free/miniconda/index.html) 

3. Install the environment. Note, this may take a few minutes. 

```bash
conda env create -f workflow/envs/sdk.yaml
```

4. Populate the `config/SDK_Data.csv` file with your scenarios

5. Initiate the workflow 

```bash 
snakemake -j6
```