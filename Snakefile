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
    # directory=".+\/",
    condition="[^\/]+",
    group="[^\/]+",
    is_grouped="(groups\/)?"

r1_fastq_suffix = config['fastq_suffix'].replace('%', '1')
r2_fastq_suffix = config['fastq_suffix'].replace('%', '2')


localrules: all, make_json_input, make_grouped_json_input, merge_grouped_tagalign, collect_tag_align_files, croo_collect_metadata, gather_qc


# condition_list = np.unique(samples["Condition"]).tolist()

# grouping_list = {}



def get_samples(condition):
    
    condition_title = ""
    condition_description = ""
    condition_samples = samples[samples['Condition']==condition]
    condition_dict = {'fastqs' : {}}

    biorep = 1
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
        condition_dict['fastqs']['rep' + str(biorep)] = {'R1' : R1_fastqs, 'R2' : R2_fastqs}
        biorep += 1
        
    condition_dict['title'] = condition_title
    condition_dict['description'] = condition_description

    return(condition_dict)


conditions_dict = {}

for condition in np.unique(samples["Condition"]).tolist():
    conditions_dict[condition] = get_samples(condition)


groupings_dict = {
    'consensus' : {
        'title' : "Consensus",
        'description' : "All samples combined",
        'conditions' : list(conditions_dict)
        }
    }



if "grouping_columns" in config.keys():
    for grouping in config["grouping_columns"]:
        for group in np.unique(samples[grouping]).tolist():
# for grouping, groups in grouping_list.items():
    # for group in groups:
            groupings_dict[group] = {
                'title' : group,
                'description' : group + " samples grouped by " + grouping,
                'conditions' : np.unique(samples[samples[grouping]==group]["Condition"]).tolist()
                }


def get_target_files(wildcards):
    target_files = ["results/qc.tsv"]
    target_files = target_files + [os.path.join("results", c, "qc/qc.json") for c in list(conditions_dict)]
    target_files = target_files + [os.path.join("results/groups", g, "qc/qc.json") for g in list(groupings_dict)]

    # target_files = target_files + [os.path.join("results", c, "croo_finished") for c in condition_list]

    return target_files


rule all:
    input: get_target_files
    run:
        print("workflow complete!")

rule gather_qc:
    input:
        [os.path.join("results", c, "qc/qc.json") for c in list(conditions_dict)],
        [os.path.join("results", "groups", g, "qc/qc.json") for g in list(groupings_dict)]
    output: "results/qc.tsv"
    # params: lambda (wildcards, input)
    shell:
        "qc2tsv {input} > {output}"



rule collect_tag_align_files:
    input: "results/{condition}/atac/metadata.json"
    output: "results/{condition}/align/{condition}.tagAlign.gz"
    shell:
        "cat results/{wildcards.condition}/align/rep*/*.trim.merged.nodup.no_chrM_MT.tn5.tagAlign.gz > {output};"


rule croo_collect_metadata:
    input: "results/{is_grouped}{condition}/atac/metadata.json"
    output: "results/{is_grouped}{condition}/qc/qc.json",
    log: "results/{is_grouped}{condition}/croo.log"
    shell:
        "croo --out-dir results/{wildcards.is_grouped}{wildcards.condition} {input} > {log} 2>&1;"
        



def cromwell_inputs(wildcards):
    inputs = {'json' : os.path.join("jsons", wildcards.is_grouped + wildcards.condition + ".json")}
    if (wildcards.is_grouped):
        inputs['tagalign'] = os.path.join("results/groups/", wildcards.condition, wildcards.condition + ".grouped.tagAlign.gz")
    return inputs

rule run_cromwell_workflow:
    input: unpack(cromwell_inputs)
    output: "results/{is_grouped}{condition}/atac/metadata.json"
    log: "results/{is_grouped}{condition}/cromwell.log"
    shell:
        "caper run {config[wdl]} -i {input.json} --out-dir results/{wildcards.is_grouped}{wildcards.condition} -m {output} "
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


def make_grouped_json_from_template(grouping, tagalign_file, json_out_file):
    json_in = open(config['template_json'], 'r')
    json_out = open(json_out_file, 'w')
    sample_json = json.load(json_in)
    group_dict = groupings_dict[grouping]

    sample_json['atac.title'] = group_dict['title']
    sample_json['atac.description'] = group_dict['description']
    sample_json['atac.tas'] = [tagalign_file]

    json.dump(sample_json, json_out, indent=4)
    json_in.close()

rule make_grouped_json_input:
    input: tagalign="results/groups/{group}/{group}.grouped.tagAlign.gz"
    output: json="jsons/groups/{group}.json"
    run: make_grouped_json_from_template(wildcards.group, input.tagalign, output.json)

rule merge_grouped_tagalign:
    input: lambda wildcards: [os.path.join("results", c, "align", c + ".tagAlign.gz") for c in groupings_dict[wildcards.group]['conditions']]
    output: "results/groups/{group}/{group}.grouped.tagAlign.gz"
    shell: "cat {input} > {output}"



