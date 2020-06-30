# Experimental branch of RASflow: RNA-Seq Analysis Snakemake Workflow
RASflow modified to meet some project-specific needs. 
Main differences with the master branch: creates a PCA plot, use local gene 
annotation files to skip slow and often failing internet queries. 
Main limitation: gene annotations limited to human genes (to be fixed later)

## To do

* Use gene annotation file for any species

## Workflow
<img src="https://github.com/zhxiaokang/RNA-Seq-analysis/blob/master/workflow/workflow_chart.jpg" width="450">

## Quick start
### Installation
Clone the repository:

`git clone https://github.com/zhxiaokang/RASflow.git`

Create the environment:

`conda env create -n rasflow -f env.yaml`

Activate the environment:

`conda activate rasflow`

### Set up configuration
Modify the metafile describing your data `configs/metadata.tsv`.

Customize the workflow based on your need in `configs/config_main.yaml`.


#### Example script to modify the config file

```bash

cd RASflow.git

DATAID=GSE126848
END=single
TRIM=no
CONTROL=[\"healthy\"]
TREAT=[\"NASH\"]

cp configs/config_main.template.yaml configs/config_main.yaml

sed -i "s/VAL_PROJECT_NAME/$DATAID/" configs/config_main.yaml
sed -i "s/VAL_TRIMMED/$TRIM/" configs/config_main.yaml
sed -i "s#VAL_READSPATH#data/datasets/$DATAID#" configs/config_main.yaml
sed -i "s#VAL_METAFILE#data/datasets/$DATAID/metadata.tsv#" configs/config_main.yaml
sed -i "s/VAL_END/$END/" configs/config_main.yaml
sed -i "s/VAL_NCORE/4/" configs/config_main.yaml
sed -i "s#VAL_FINALOUTPUT#../../RASflowResults#" configs/config_main.yaml
sed -i "s#VAL_TRANS#data/ref/GRCh38.99/Homo_sapiens.GRCh38.cdna.all.fa.gz#" configs/config_main.yaml
sed -i "s#VAL_GENOME#data/ref/GRCh38.99/Homo_sapiens.GRCh38.dna_sm.alt.fa.gz#" configs/config_main.yaml
sed -i "s#VAL_ANNOTATION#data/ref/GRCh38.99/Homo_sapiens.GRCh38.99.gtf.gz#" configs/config_main.yaml
sed -i "s/VAL_PAIR/FALSE/" configs/config_main.yaml
sed -i "s/VAL_CONTROL/$CONTROL/" configs/config_main.yaml
sed -i "s/VAL_TREAT/$TREAT/" configs/config_main.yaml
sed -i "s/VAL_EnsemblDataSet/hsapiens_gene_ensembl/" configs/config_main.yaml
```

### Run RASflow
`python main.py`

or the following to skip the question

`python main.py <<< y`

### Additional results

```bash
Rscript scripts/merge_deseq2.R # see all-samples* in output/project/trans/dea/coutGroup
Rscript scripts/merge_tpm.R    # see all-samples* in output/project/trans/tpmFile
```

## Tutorial
A more detailed tutorial of how to use this workflow can be found here: [Tutorial](https://github.com/zhxiaokang/RASflow/blob/master/Tutorial.pdf)

## Evaluation
RASflow has been evaluated on 4 datasets including two model organisms (human and mouse) and a non-model organism (Atlantic cod). To keep this repository as light as possible, the evaluation of RASflow on real datasets is deposited here: [RASflow_realData](https://git.app.uib.no/Xiaokang.Zhang/rasflow_realdata)

## References
RASflow is available as a preprint on bioRxiv: [Zhang X, Jonassen I. RASflow: An RNA-Seq Analysis Workflow with Snakemake. bioRxiv. 2019 Nov 839191.](https://www.biorxiv.org/content/10.1101/839191v1) Recently accepted by BMC Bioinformatics.
