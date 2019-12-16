# ENCODE ATAC-seq Workflow Wrapper 

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.com/robinmeyers/atac-encode-snakemake.svg?branch=master)](https://travis-ci.com/robinmeyers/atac-encode-snakemake)


## Authors

* Robin Meyers [@robinmeyers](https://github.com/robinmeyers)

## Purpose

**Under development**

This pipeline acts as a wrapper around the ENCODE ATAC-seq workflow for running multiple samples and gathering QC metrics.

## Usage

### Simple

#### Step 1: Install workflow


Install the [ATAC workflow](https://github.com/ENCODE-DCC/atac-seq-pipeline) as well as the ```croo``` and ```qc2tsv``` utilities.


Clone this repositiory and change into the new directory.

Activate the conda environment. Run the snakemake pipeline on provided test data.

```
$ conda activate encode-atac-seq-pipeline
$ snakemake --directory .test
```

Examine the outputs of the workflow in the directory ```.test/outs/```


#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

#### Step 3: Execute workflow

Ensure the correct conda environment is active with

```
$ conda activate encode-atac-seq-pipeline
```

Test your configuration by performing a dry-run via

```
$ snakemake -n
```

```
$ snakemake -j 999 -p  --cluster-config cluster.json --profile slurm
```

Before running in a cluster environment, edit the `cluter.json` config file.

# Step 4: Investigate results

The main outputs are in the  ```<sample>/outs/``` directories. 

After successful execution, you can create a self-contained interactive HTML report with workflow statistics.

```$ snakemake --report report.html```

This report can, e.g., be forwarded to your collaborators.

### Updating the workflow

If you installed the workflow by cloning the github repo, you can pull latest updates to workflow with 

```$ git pull --rebase```

This will require you to first commit any changes you made to your configuration file before pulling new updates.

Then simply rerun the `snakemake` command.


## Testing

Tests cases are in the subfolder `.test`. They are automtically executed via continuous integration with Travis CI.
