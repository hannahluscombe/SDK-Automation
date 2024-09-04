"""General rules for model building and solving"""

def get_otoole_results(otoole_config: str) -> list[str]:
    """Gets result files to be created"""

    with open(otoole_config) as f:
        otoole = yaml.safe_load(f)

    results = [x for x in otoole if otoole[x]["type"] == "result"]

    missing = [
        "NewStorageCapacity",
        "NumberOfNewTechnologyUnits",
        "SalvageValueStorage",
        "StorageLevelDayTypeStart",
        "StorageLevelDayTypeFinish",
        "StorageLevelSeasonStart",
        "StorageLevelSeasonFinish",
        "StorageLevelYearStart",
        "StorageLevelYearFinish",
        "Trade",
    ]

    return [x for x in results if x not in missing]

OTOOLE_FILES = get_otoole_results(config["otoole"])

rule simplify_data:
    message:"Simplifying data for {wildcards.scenario}"
    input:
        txt = "results/{scenario}/{scenario}.txt"
    output:
        txt = "results/{scenario}/{scenario}_simple.txt"
    script: 
        "../scripts/simplify.py"

rule update_data:
    message:"Updating data for pre-processing"
    params: 
        config = "resources/otoole.yaml"
    input:
        txt = "results/{scenario}/{scenario}_simple.txt"
    output:
        txt = "results/{scenario}/{scenario}_updated_data.txt"
    script:
        "../scripts/update_data.py"

rule preprocess_data:
    message:"Pre-processing data for {wildcards.scenario}"
    params:
        data_format = "otoole" # (momani|otoole)
    input:
        txt = "results/{scenario}/{scenario}_updated_data.txt"
    output:
        txt = "results/{scenario}/{scenario}_preprocessed.txt"
    script: 
        "../scripts/preprocess.py"

rule build_model:
    message:"Building LP file for {wildcards.scenario}"
    params:
        model = "resources/osemosys_pp.txt"
    input:
        data = "results/{scenario}/{scenario}_preprocessed.txt"
    output:
        lp = temp("results/{scenario}/{scenario}.lp")
    shell: 
        "glpsol -m {params.model} -d {input.data} --wlp {output.lp} --check"


rule solve_model:
    message:"Solving model for {wildcards.scenario}"
    params:
        solver = config["solver"]
    input:
        lp = "results/{scenario}/{scenario}.lp"
    output:
        sol = temp("results/{scenario}/{scenario}.sol")
    shell: 
        """
        if [ {params.solver} = gurobi ]
        then
          gurobi_cl Method=2 ResultFile={output.sol} {input.lp}
        elif [ {params.solver} = cplex ]
        then
          cplex -c "read {input.lp}" "optimize" "write {output.sol}"
        else
          cbc {input.lp} solve -sec 1500 -solu {output.sol}
        fi
        """

rule process_results:
    message:"Processing results for {wildcards.scenario}"
    params:
        otoole_config = "resources/otoole.yaml",
        results_dir = "results/{scenario}/results", 
        solver = config["solver"]
    input:
        sol = "results/{scenario}/{scenario}.sol",
        data = "results/{scenario}/{scenario}_simple.txt" # do not use pre-processed! 
    output:
        expand("results/{{scenario}}/results/{csv}.csv", csv=OTOOLE_FILES)
    shell: 
        "otoole results {params.solver} csv {input.sol} {params.results_dir} datafile {input.data} {params.otoole_config}"