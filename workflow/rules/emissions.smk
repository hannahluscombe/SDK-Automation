"""Generates Emission Abatement Curves"""


EMISSION_TARGETS = pd.read_csv(config["emission_config"], index_col="emission_reduction")
EMISSION_SCENARIOS = EMISSION_TARGETS.index.to_list()

rule copy_emission_scenario_base:
    message: "Copying base emission scenario and results"
    params:
        in_dir = "results/{scenario}/results",
        out_dir = "results/{scenario}/emission_curves/0_percent_reduction"
    input: 
        txt = "results/{scenario}/{scenario}_simple.txt",
        results = expand("results/{{scenario}}/results/{csv}.csv", csv=OTOOLE_FILES),
    output:
        txt = "results/{scenario}/emission_curves/0_percent_reduction/data.txt",
        results = expand("results/{{scenario}}/emission_curves/0_percent_reduction/results/{csv}.csv", csv=OTOOLE_FILES),
    shell:
        """
        cp {input.txt} {output.txt} && cp -r {params.in_dir} {params.out_dir}
        """

rule create_emission_scenario_templates:
    message: "Creating template emission reduction scenarios"
    wildcard_constraints:
        # emission_reduction=r"^(?!0$)\d+$"
        emission_reduction="|".join([str(x) for x in EMISSION_SCENARIOS if x != 0])
    input: 
        txt = "results/{scenario}/{scenario}_simple.txt"
    output:
        txt = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}.txt"
    shell:
        "cp {input.txt} {output.txt}"

rule create_emission_reduction_dataframe:
    message:"Creating Emission Reduction Scenarios for {wildcards.scenario}"
    wildcard_constraints:
        # emission_reduction=r"^(?!0$)\d+$"
        emission_reduction="|".join([str(x) for x in EMISSION_SCENARIOS if x != 0])
    params: 
        config = "resources/otoole.yaml",
        base_year = config["emission_base_year"],
        reduction_method = config["emission_reduction_method"],
        region = "RE1",
        emission = "EMIC02"
    input:
        base_emissions = "results/{scenario}/emission_curves/0_percent_reduction/results/AnnualEmissions.csv"
    output:
        emission_limits = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/updates/AnnualEmissionLimit.csv"
    script:
        "../scripts/create_emission_curve.py"

rule update_annual_emission_limit:
    message:"Updating Emission Limit for scenario {wildcards.scenario} and {wildcards.emission_reduction}%% C02 reduction"
    wildcard_constraints:
        # emission_reduction=r"^(?!0$)\d+$"
        emission_reduction="|".join([str(x) for x in EMISSION_SCENARIOS if x != 0])
    params: 
        config = "resources/otoole.yaml",
        parameter = "AnnualEmissionLimit"
    input:
        csv = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/updates/AnnualEmissionLimit.csv",
        txt = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}.txt"
    output:
        txt = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}_co2L.txt"
    script:
        "../scripts/update_data.py"


rule preprocess_emission_data:
    message:"Pre-processing data for scenario {wildcards.scenario} and {wildcards.emission_reduction}%% C02 reduction"
    wildcard_constraints:
        emission_reduction="|".join([str(x) for x in EMISSION_SCENARIOS if x != 0])
        # emission_reduction=r"^(?!0$)\d+$"
    params:
        data_format = "otoole" # (momani|otoole)
    input:
        txt = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}_co2L.txt"
    output:
        txt = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}_co2L_pp.txt"
    script: 
        "../scripts/preprocess.py"

rule build_emission_model:
    message:"Building LP file for scenario {wildcards.scenario} and {wildcards.emission_reduction}%% C02 reduction"
    wildcard_constraints:
        emission_reduction="|".join([str(x) for x in EMISSION_SCENARIOS if x != 0])
        # emission_reduction=r"^(?!0$)\d+$"
    params:
        model = "resources/osemosys_pp.txt"
    input:
        data = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}_co2L_pp.txt"
    output:
        lp = temp("results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}.lp")
    resources:
        mem_mb=2000
    shell: 
        "glpsol -m {params.model} -d {input.data} --wlp {output.lp} --check"


rule solve_emission_model:
    message:"Solving model for scenario {wildcards.scenario} and {wildcards.emission_reduction}%% C02 reduction"
    wildcard_constraints:
        emission_reduction="|".join([str(x) for x in EMISSION_SCENARIOS if x != 0])
        # emission_reduction=r"^(?!0$)\d+$"
    params:
        solver = config["solver"]
    input:
        lp = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}.lp"
    output:
        sol = temp("results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}.sol")
    resources:
        mem_mb=2000
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

rule process_emission_results:
    message:"Processing results for scenario {wildcards.scenario} and {wildcards.emission_reduction}%% C02 reduction"
    wildcard_constraints:
        emission_reduction="|".join([str(x) for x in EMISSION_SCENARIOS if x != 0])
        # emission_reduction=r"^(?!0$)\d+$"
    params:
        otoole_config = "resources/otoole.yaml",
        results_dir = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/results", 
        solver = config["solver"]
    input:
        sol = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}.sol",
        data = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/{scenario}_co2L.txt" # do not use pre-processed! 
    output:
        expand("results/{{scenario}}/emission_curves/{{emission_reduction}}_percent_reduction/results/{csv}.csv", csv=OTOOLE_FILES)
    resources:
        mem_mb=8000
    shell: 
        "otoole results {params.solver} csv {input.sol} {params.results_dir} datafile {input.data} {params.otoole_config}"

rule extract_annual_emission_results:
    message: "Extracting Annual Emissions for scenario {wildcards.scenario}"
    params:
        result = "AnnualEmissions"
    input:
        expand("results/{{scenario}}/emission_curves/{emission_reduction}_percent_reduction/results/AnnualEmissions.csv", emission_reduction=EMISSION_SCENARIOS)
    output:
        csv = "results/{scenario}/emission_curves/AnnualEmissions.csv"
    script: 
        "../scripts/extract_emission_result.py"

rule extract_discounted_cost_results:
    message: "Extracting Total Cost for scenario {wildcards.scenario}"
    params:
        result = "TotalDiscountedCost"
    input:
        expand("results/{{scenario}}/emission_curves/{emission_reduction}_percent_reduction/results/TotalDiscountedCost.csv", emission_reduction=EMISSION_SCENARIOS)
    output:
        csv = "results/{scenario}/emission_curves/TotalDiscountedCost.csv"
    script: 
        "../scripts/extract_emission_result.py"

rule plot_abatement_cost_curve:
    message: "Plotting Abatement Cost Curve for {wildcards.scenario}"
    params:
        year = config["emission_base_year"]
    input:
        annual_emissions = "results/{scenario}/emission_curves/AnnualEmissions.csv",
        total_cost = "results/{scenario}/emission_curves/TotalDiscountedCost.csv"
    output:
        plot = "results/{scenario}/emission_curves/abatement_cost_curve.png"
    script: 
        "../scripts/plot_abatement_cost_curve.py"