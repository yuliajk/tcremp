# TCRemP

TCRemP is a package developed to perform T-cel receptor (TCR) embedding. The package utilizes prototypes, commonly encountered TCRs either sampled from a probabilistic V(D)J rearrangement model (see Murugan et al. 2012) or a large pool of individual TCR repertoires from the population.

The workflow is the following:
* TCRemP pipeline starts with a selection of ``k`` prototype TCR alpha and beta sequences, then it computes the distances from every of ``n`` input TCR alpha-beta pairs to ``2 * k`` prototypes for V, J and CDR3 regions, resulting in ``6 * k`` parameters (or ``3 * k`` for cases when only one of the chains is present).
* Resulting distances are treated as embedding co-ordinates and and are subject to principal component analysis (PCA). One can monitor the information conveyed by each PC, whether they are related to features such as Variable or Joining genes, CDR3 region length or a certain epitope.
> N.B. TCRemP is currently in active development, please submit the below for documentation and a proof-of-concept example.

Using TCRemP one can:
- perform an embedding for a set of T-cell clonotypes, defined by TCR’s Variable (V) and Joining (J) segment IDs and complementarity determining region 3 (CDR3, amino acid sequence placed at the V-J junction). The embedding is performed by mapping those features to real vectors using similarities to a set of **prototype** TCR sequences
- embed a set of clones, pairs of TCR alpha and beta chain clonotypes
- analyze the mapping by performing dimensionality reduction and evaluating principal components (PCs)
- cluster the embeddings using [DBSCAN](https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html) method with parameter selection using knee/elbow method
- visualize T-cell clone and clonotype embeddings using tSNE, coloring the visualization by user-specified clonotype labels, such as antigen specificities
- infer cluster that are significantly enriched in certain labels, e.g. TCR motifs belonging to CD8+ T-cell subset or specific to an antigen of interest

Planned features:
- (in progress) co-embed samples with  [VDJdb database](https://github.com/antigenomics/vdjdb-db) to predict TCRs associated with certain antigens, i.e. “annotate” TCR repertoires
- (in progress) perform imputation to correctly handle mixed single-/paired-chain data
- (in progress) implement B-cell receptor (BCR/antibody) prototypes to apply the method to antibody sequencing data

# Getting started

## Installation procedure

One can simply install the software using [pip](https://pypi.org/project/pip/):

```{bash}
pip install git+https://github.com/antigenomics/tcremp.git@0.0.1
```

Or clone the repository via git, make corresponding [conda](https://docs.conda.io/en/latest/) and install with requirements:

```{bash}
git clone https://github.com/antigenomics/mirpy.git
git clone https://github.com/antigenomics/tcremp.git
cd tcremp
conda create -n tcremp_env ipython python=3.11
conda activate tcremp_env   # or: "source activate tcremp_env" depending on your conda setup
pip install -r requirements.txt
```


## Preparing the input data

The input data typically consists of a table containing clonotypes as defined above, either TCR alpha, or beta, or both. One can additionally tag clonotypes/clones with user-defined ids, e.g. cell barcodes, and labels, e.g. antigen specificity or phenotype. One can also use a custom clonotype table instead of a pre-built set of prototypes (see data/data_prebuilt/VDJdb_data_paired_example.csv).

  
### Input format

#### Common requirements

1. V and J segment names should be provided based on [IMGT](https://www.imgt.org/) naming, e.g. ``TRAV35*03`` or ``TRBV11-2``. TCRemP will always use the major allele, so the alleles above will be transformed into ``TRBV*01`` and ``TRBV11*01``  
2. The data should not contain any missing data for any of the columns: V, J and CDR3. Although for paired-chain format one of the chains may be missing.
3. There should be no symbols except for 20 amino acids in CDR3s

#### Input columns
| Column name | Description | Required |
| ----------- | ----------- | ----------- |
| a_cdr3aa | TCR alpha chain cdr3 | required |
| a_v | TCR alpha V segment ID | required |
| a_j | TCR alphasegment ID | required |
| b_cdr3aa | TCR beta chain cdr3 | required |
| b_v | TCR beta V segment | required |
| b_j | TCR beta J segment | required |
| clonotype_index | user provided clonotype id, used in  | optional |
| label | user provided label, used in cluster post-analysis | optional |


#### Single chain table example example

Either wide with missing values

| clonotype_index | a_cdr3aa | a_v | a_j | b_cdr3aa | b_v | b_j | label |
| :---:   | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| GACTGCGCATCGTCGG-28   | CAGHTGNQFYF | TRAV35	 | TRAJ49 |   |  |  |  IVTDFSVIK |
| GACTGCGCATCGTCGG-28   |  |  |   | CASSWGGGSHYGYTF | TRBV11-2  | TRBJ1-2 |  IVTDFSVIK |


#### Paired chain example

A simple wide format

| data_id | a_cdr3aa | a_v | a_j | b_cdr3aa | b_v | b_j | label |
| :---:   | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| GACTGCGCATCGTCGG-28   | CAGHTGNQFYF | TRAV35	 | TRAJ49 | CASSWGGGSHYGYTF | TRBV11-2  | TRBJ1-2 | IVTDFSVIK |

 
## Running TCRemP

### Basic usage

Run the tool as 

```{bash}
python tcremp_run.py --input my_input_data.txt --output my_folder --chain TRA_TRB
```

The command above will:
- checks input data format and proofreads the dataset
- extracts TCR alpha and beta clonotypes from my_input_data.txt
- calculates distance scores from clonotypes for the built-in set of 3000 prototypes for each chain
- performs PCA and saves transformed data 
- runs DBSCAN clustering with parameters min_samples = 2, eps value detected by knee method,  and saves resulting clusters
All input will be saved in ``my_folder/``

If one runs

```{bash}
python tcremp_run.py --input my_input_data.txt --output my_folder --chain TRA_TRB --label epitope
```

the resulting clusters will be saved with both user-provided labels and cluster labels

The following command will skip the clustering and only save embeddings:

```{bash}
python tcremp_run.py --input my_input_data.txt --chain TRA_TRB --clstr_model none
```
 
### Command line parameters

The parameters for running ``tcremp`` main script are the following:

| parameter | short usage | description | available values | required | default value |
| --- | --- | --- | --- | --- | --- |
| --input | -i | input clonotype table  | path to file | yes | - |
| --chain | -c | single or paired clonotype chains | TRA, TRB, TRA_TRB | yes | - |
| --output | -o | pipeline output folder | path to directory | no | tcremp_{inputfilename}/ |
| --clstr_model | -m | clustering model to be used: dbscan, kmeans or run without clustering step (value 'none') | none, dbscan, kmeans | no | dbscan |
| --label | -l | name of the input file column with data classes for clustering. If provided, this value will be added to clustering output table and label of each cluster will be defined | str | no | - |
| --species | -s | species of built-in prototypes to be used | HomoSapiens, (MusMusculus - planned) | no | HomoSapiens |
| --n | -n | number of prototypes to be selected for embedding supplemented prototype table | integer | no | 3000 |
| --random | -r | random seed for random prototype selection. If not provided or 0, first n prototypes are selected | integer | no | 0 |
| --prototypes_path | -p | path to the custom input prototype table | path to file | no | - |
| --clonotype_index | -d | column containing user-provided clonotype id in input data, will be transfered to output tables | str | no | - |

 
### Output

All output files will contain the following **common columns**:
- tcremp_id - assigned identifier to each row of the input table
- cloneId - identifier of unique extracted clonotype. In case of single chain pipeline, stands for TRA or TRB clonotype, in case of paired chain pipeline stands for paired TRA and TRB clonotype.

The output folder will contain the following files:

| File name | description |
| --- | --- |
| tcremb_dists_{chain}.txt | Embeddings (the set of distances to prototypes) for each input clonotype row in the table, the chain can be: TRA, TRB or TRA_TRB for paired-chain format  |
| tcremb_pca_{chain}.txt | Principal components of the embedding, the chain can be: TRA, TRB or TRA_TRB for paired-chain format|
| tcremb_tsne_{chain}.txt | Calculated tSNE co-ordinates for each input clonotype row in the table, the chain can be: TRA, TRB or TRA_TRB for paired-chain format|
| tcremb_clstr_res_{chain}.txt | Cluster assignments for each input clonotype row in the table. If the cluster is -1, the clonotype does not belong to any cluster (can be treated as “noise” in some settings). Chain can be: TRA, TRB or TRA_TRB |
| filtered_out_data.txt | rows that were excluded during data cleansing with exclusion reason |
| clonotypes_{chain}.txt | Clonotypes used in the single-chain pipeline, chain can be: TRA or TRB |
| clonotypes_paired\_{chain}.txt | Clonotypes used in the paired-chain pipeline |
| prototypes{chain}_{n}.txt | prototypes on which the distances was calculated, where n - number of selected prototypes and the chainhain can be: TRA or TRB. The file is absent if default subset of prototypes was used for calculation|
| res_{chain}.txt | Raw distance table for single chain mapping. Chain can be: TRA or TRB |
| res_paired_{chain}.txt | Raw distance table for paired chain mapping. Chain can be: TRA or TRB |
