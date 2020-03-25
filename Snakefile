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


condition_list = np.unique(samples["Condition"]).tolist() + ["consensus"]

def get_target_files(wildcards):
    target_files = ["results/qc.tsv"]

    # target_files = target_files + [os.path.join("results", c, "croo_finished") for c in condition_list]

    return target_files


rule all:
    input: get_target_files
    run:
        print("workflow complete!")

rule gather_qc:
    input: [os.path.join("results", c, "qc/qc.json") for c in condition_list]
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
    rep_i = 1
    condition_title = "Consensus" if condition == "consensus" else ""
    condition_description = "All samples combined" if condition == "consensus" else ""
    condition_samples = samples if condition == "consensus" else samples[samples['Condition']==condition]
    for sample, row in condition_samples.iterrows():
        R1_fastqs = []
        R2_fastqs = []
        for fastq_dir in config['fastq_dirs']:
            R1_fastqs = R1_fastqs + \
                glob.glob(os.path.join(fastq_dir, sample + r1_fastq_suffix)) + \
                glob.glob(os.path.join(fastq_dir, sample, sample + r1_fastq_suffix))
            R2_fastqs = R2_fastqs + \
                glob.glob(os.path.join(fastq_dir, sample + r2_fastq_suffix)) + \
                glob.glob(os.path.join(fastq_dir, sample, sample + r2_fastq_suffix))
        condition_title = row["Title"] if (condition_title == "" and row["Title"] != "") else condition_title
        condition_description = row["Description"] if condition_description == "" and row["Description"] != "" else condition_description
        sample_json['atac.fastqs_rep' + str(rep_i) + '_R1'] = R1_fastqs
        sample_json['atac.fastqs_rep' + str(rep_i) + '_R2'] = R2_fastqs
        if condition_title != "":
            sample_json['atac.title'] = condition_title
        if condition_description != "":
            sample_json['atac.description'] = condition_description
        rep_i += 1
    json.dump(sample_json, json_out, indent=4)
    json_in.close()
    # json_out.close()


rule make_json_input:
    input:
    output: json="jsons/{condition}.json"
    run: make_json_from_template(wildcards.condition, output.json)




