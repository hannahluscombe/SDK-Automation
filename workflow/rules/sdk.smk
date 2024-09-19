"""Rules for working with SDK data"""

def get_url(wildcards):
    return SDKS.at[wildcards.scenario, "data"]

def get_custom_file(wildcards):
    f = SDKS.at[wildcards.scenario, "data"]
    return f"config/sdk_custom/{f}"

if config["run"] == "custom":

    rule copy_xlsm:
        message:"Copying xlsm {wildcards.scenario} from Zenodo"
        input:
            xlsm = get_custom_file
        output:
            xlsm = "results/{scenario}/{scenario}.xlsm"
        shell: 
            "cp {input.xlsm} {output.xlsm}"

elif config["run"] == "zenodo":

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