# VCANN-WDL
"Variant Calling Annotated, WDL" (VCANN-WDL)
#
Parallel SARS-CoV-2 Variant Calling Annotated with WDL orchestration.
This is a very simple workflow used for demonstration purposes.
This pipeline makes use of three docker images found in this repository: https://github.com/Cameron-Nguyen1/Docker_Images

## I. Introduction
This repository is just an example of WDL orchestration of a SARS-CoV-2 variant calling pipeline. In practice, this could be used to annotate variant calls from any organism by making the appropriate modifications to the variant calling docker image
or by adding other inputs to the WDL workflow.

## II. Input Arguments
Arguments are neatly formatted in "wf_params.json". Edit them as necessary, the pipeline is setup to run as-is.
```
sraCodes = SRA accessions to perform variant calling on.
cpu = CPUs per process.
mem = Memory per process.
read_len_required = Read length required to be considered for alignment.
```

## III. Example Usage
```
miniwdl run VCANN.wdl --input wf_params.json
```
