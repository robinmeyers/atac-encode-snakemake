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


localrules: all, make_json_input, croo_collect_metadata


condition_list = np.unique(samples["Condition"])


def get_target_files(wildcards):
    target_files = []

    target_files = target_files + [os.path.join("results", c, "croo_finished") for c in condition_list]

    return target_files


rule all:
    input: get_target_files
    run:
        print("workflow complete!")


rule croo_collect_metadata:
    input: "results/{condition}/atac/metadata.json"
    output: "results/{condition}/croo_finished"
    shell: "croo {input} && touch {output}"


rule run_cromwell_workflow:
    input: "jsons/{condition}.json"
    output: "results/{condition}/atac/metadata.json"
    log: "results/{condition}/cromwell.log"
    shell:
        "caper run {config[wdl]} -i {input} --out-dir results/{wildcards.condition} -m {output} "
        "> {log} 2>&1"


def make_json_from_template(condition, json_out):
    json_in = open(config['template_json'], 'r')
    json_out = open(json_out_file, 'w')
    json_spec = json.load(json_in)
    # for fastq_dir in config['fastq_dirs']:
        # for id, row in samples[samples['Sample']==sample].iterrows():
        # for fastq_dir in config['fastq_dirs']:
        #     if os.path.isdir(os.path.join(fastq_dir, id))
    json.dump(json_spec, json_out_file)
    json_in.close()
    json_out.close()


rule make_json_input:
    input:
    output: "jsons/{condition}.json"
    run: make_json_from_template(wildcards.condition, output)




