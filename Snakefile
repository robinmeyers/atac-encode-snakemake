# The main entry point of your workflow.
# After configuring, running snakemake -n in a clone of this repository should successfully execute a dry-run of the workflow.


import os
import glob
import numpy as np
import pandas as pd
from snakemake.utils import validate, min_version


##### load config and sample sheets #####



configfile: "config.yaml"
# report: "report/workflow.rst"

validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_csv(config['samplesheet']).set_index("Sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")


# Allow users to fix the underlying OS via singularity.
singularity: "docker://continuumio/miniconda3"

# include: "rules/pdna.smk"

wildcard_constraints:
    directory=".+\/",
    condition="[^\/]+",
    sample="[^\/]+"

r1_fastq_suffix = config['fastq_suffix'].replace('%', '1')
r2_fastq_suffix = config['fastq_suffix'].replace('%', '2')


localrules: all, make_json_input, croo_collect_metadata, gather_qc


condition_list = np.unique(samples["Condition"]).tolist()

grouping_list = {}

if "grouping_columns" in config.keys():
    for grouping in config["grouping_columns"]:
        grouping_list[grouping] = np.unique(samples[grouping]).tolist()


def get_samples(condition, grouping="Condition", use_bioreps=True):
    biorep = 1
    
    if grouping == "Condition":
        condition_title = ""
        condition_description = ""
        condition_samples = samples[samples['Condition']==condition]
    elif grouping == "Consensus":
        condition_title = "Consensus"
        condition_description = "All samples combined"
        condition_samples = samples
    else:
        condition_title = condition
        condition_description = condition + " samples grouped by " + grouping
        condition_samples = samples[samples[grouping]==condition]

    condition_dict = {'fastqs' : {}}
    R1_fastqs = []
    R2_fastqs = []
    for sample, row in condition_samples.iterrows():

        for fastq_dir in config['fastq_dirs']:
            R1_fastqs = R1_fastqs + \
                glob.glob(os.path.join(fastq_dir, sample + r1_fastq_suffix)) + \
                glob.glob(os.path.join(fastq_dir, sample, sample + r1_fastq_suffix))
            R2_fastqs = R2_fastqs + \
                glob.glob(os.path.join(fastq_dir, sample + r2_fastq_suffix)) + \
                glob.glob(os.path.join(fastq_dir, sample, sample + r2_fastq_suffix))
        condition_title = row["Title"] if (condition_title == "" and row["Title"] != "") else condition_title
        condition_description = row["Description"] if condition_description == "" and row["Description"] != "" else condition_description

        if use_bioreps:
            condition_dict['fastqs']['rep' + str(biorep)] = {'R1' : R1_fastqs, 'R2' : R2_fastqs}
            biorep += 1
            R1_fastqs = []
            R2_fastqs = []

    condition_dict['title'] = condition_title
    condition_dict['description'] = condition_description
    if not use_bioreps:
        condition_dict['fastqs']['rep1'] = {'R1' : R1_fastqs, 'R2' : R2_fastqs}


    return(condition_dict)


conditions_dict = {}

for condition in condition_list:
    conditions_dict[condition] = get_samples(condition)
for grouping, groups in grouping_list.items():
    for group in groups:
        conditions_dict[group] = get_samples(group, grouping=grouping, use_bioreps=False)
conditions_dict["consensus"] = get_samples("consensus", grouping="Consensus", use_bioreps=False)








def get_target_files(wildcards):
    target_files = ["results/qc.tsv"]

    # target_files = target_files + [os.path.join("results", c, "croo_finished") for c in condition_list]

    return target_files


rule all:
    input: get_target_files
    run:
        print("workflow complete!")

rule gather_qc:
    input: [os.path.join("results", c, "qc/qc.json") for c in conditions_dict.keys()]
    output: "results/qc.tsv"
    # params: lambda (wildcards, input)
    shell:
        "qc2tsv {input} > {output}"


rule croo_collect_metadata:
    input: "results/{condition}/atac/metadata.json"
    output: "results/{condition}/qc/qc.json"
    log: "results/{condition}/croo.log"
    shell: "croo --out-dir results/{wildcards.condition} {input} > {log} 2>&1"


rule run_cromwell_workflow:
    input: "jsons/{condition}.json"
    output: "results/{condition}/atac/metadata.json"
    log: "results/{condition}/cromwell.log"
    shell:
        "caper run {config[wdl]} -i {input} --out-dir results/{wildcards.condition} -m {output} "
        "> {log} 2>&1"


def make_json_from_template(condition, json_out_file):
    json_in = open(config['template_json'], 'r')
    json_out = open(json_out_file, 'w')
    sample_json = json.load(json_in)
    condition_dict = conditions_dict[condition]

    sample_json['atac.title'] = condition_dict['title']
    sample_json['atac.description'] = condition_dict['description']
    for rep in condition_dict['fastqs'].keys():
        sample_json['atac.fastqs_' + rep + "_R1"] = condition_dict['fastqs'][rep]['R1']
        sample_json['atac.fastqs_' + rep + "_R2"] = condition_dict['fastqs'][rep]['R2']

    json.dump(sample_json, json_out, indent=4)
    json_in.close()



rule make_json_input:
    input:
    output: json="jsons/{condition}.json"
    run: make_json_from_template(wildcards.condition, output.json)




