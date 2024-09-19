"""Rules for extracting minfin results"""

MINFIN_FILES = config["minfin_files"]
# quite an ugly implementation to track otoole files, but gets the job done i guess

rule minfin_results:
    message:"Formating result data for {wildcards.scenario}"
    params:
        minfin_files = expand("{csv}.csv", csv=MINFIN_FILES),
        in_dir = "results/{scenario}/results/",
        out_dir = "results/{scenario}/MinFin/"
    input:
        expand("results/{{scenario}}/results/{csv}.csv", csv=OTOOLE_RESULTS)
    output:
        expand("results/{{scenario}}/MinFin/{csv}.csv", csv=MINFIN_FILES)
    script:
        "../scripts/minfin_results.py"