"""Generates Emission Abatement Curves"""


EMISSION_TARGETS = pd.read_csv(config["emission_config"], index_col="emission_reduction")
EMISSION_SCENARIOS = EMISSION_TARGETS.index.to_list() 

rule copy_emission_scenario_base:
    message: "Copying base emission scenario and results"
    input: 
        txt = "results/{scenario}/{scenario}_simple.txt"
        results = "results/{scenario}/results/"
    output:
        txt = "results/{scenario}/emission_curves/0_percent_reduction/data.txt"
        results = "results/{scenario}/emission_curves/0_percent_reduction/results"
    shell:
        """
        cp {input.txt} {output.txt} && cp -r {input.results} {output.results}
        """

rule create_emission_scenario_templates:
    message: "Creating template emission reduction scenarios"
    input: 
        txt = "results/{scenario}/{scenario}_simple.txt"
    output:
        txt = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/data.txt"
    shell:
        "cp {input.txt} {output.txt}"

rule create_emission_reduction_dataframe:
    message:"Creating Emission Reduction Scenarios for {wildcards.scenario}"
    params: 
        config = "resources/otoole.yaml"
        base_year = config["emission_base_year"]
        reduction_method = config["emission_reduction_method"]
        region = "RE1"
        emission = "EMIC02"
    input:
        base_emissions = "results/{scenario}/emission_curves/0_percent_reduction/results/AnnualEmissions.csv"
    output:
        emission_limits = "results/{scenario}/emission_curves/{emission_reduction}_percent_reduction/curve.csv"
    script:
        "../scripts/create_emission_curve.py"
