"""General rules for model building and solving"""


rule simplify_data:
    message:"Simplifying data for scenario {wildcards.scenario}"
    input:
        txt = "results/{scenario}/{scenario}.txt"
    output:
        txt = "results/{scenario}/{scenario}_simple.txt"
    script: 
        "../scripts/simplify.py"

rule create_new_discount_rate:
    message:"Creating new discount rate data for {wildcards.scenario}"
    params:
        region = "RE1",
        discount_rate = "0.05"
    output:
        csv = "results/{scenario}/updates/DiscountRate.csv"
    run:
        data = [[params.region, params.discount_rate]]
        df = pd.DataFrame(data, columns=["REGION", "VALUE"])
        df.to_csv(output.csv, index=False)

rule update_discount_rate:
    message:"Updating data for scenario {wildcards.scenario}"
    params: 
        config = "resources/otoole.yaml",
        parameter = "DiscountRate",
        save_dir = "results/{scenario}/data/",
        data_input_type = "datafile"
    input:
        csv = "results/{scenario}/updates/DiscountRate.csv",
        data_in = "results/{scenario}/{scenario}_simple.txt"
    output:
        csvs = expand("results/{{scenario}}/data/{csv}.csv", csv=OTOOLE_DATA)
    script:
        "../scripts/update_data.py"

rule create_datafile:
    message:"Creating datafile for scenario {wildcards.scenario}"
    params: 
        config = "resources/otoole.yaml",
        data_dir = "results/{scenario}/data/"
    input:
        csv = expand("results/{{scenario}}/data/{csv}.csv", csv=OTOOLE_DATA)
    output:
        txt = "results/{scenarios}/{scenario}_updated.txt"
    shell:
        """
        otoole convert csv datafile {params.data_dir} {output.txt} {params.config}
        """

rule preprocess_data:
    message:"Pre-processing data for scenario {wildcards.scenario}"
    params:
        data_format = "otoole" # (momani|otoole)
    input:
        txt = "results/{scenario}/{scenario}_updated.txt"
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
        data_dir = "results/{scenario}/data/",
        results_dir = "results/{scenario}/results", 
        solver = config["solver"]
    input:
        sol = "results/{scenario}/{scenario}.sol",
        data = expand("results/{{scenario}}/data/{csv}.csv", csv=OTOOLE_DATA)
    output:
        expand("results/{{scenario}}/results/{csv}.csv", csv=OTOOLE_RESULTS)
    shell: 
        "otoole results {params.solver} csv {input.sol} {params.results_dir} csv {params.data_dir} {params.otoole_config}"