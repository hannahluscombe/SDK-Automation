"""Rules for working with SDK data"""

def get_url(wildcards):
    return URLS.at[wildcards.scenario, "URL"]

rule download_xlsm:
    message:"Downloading xlsm {wildcards.scenario} from Zenodo"
    params:
        url = get_url
    output:
        xlsm = "results/{scenario}/{scenario}.xlsm"
    script: 
        "../scripts/download_xlsm.py"

rule extract_data:
    message:"Extracting txt data for {wildcards.scenario}"
    input:
        xlsm = "results/{scenario}/{scenario}.xlsm"
    output:
        txt = "results/{scenario}/{scenario}.txt"
    script: 
        "../scripts/extract_data.py"