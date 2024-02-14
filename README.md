# TCRemb
TCRemb is a pakage depeloped for T cell gene receptors comparison.
TCRemb is in the in active development right now.

With TCRemb you can:
- cluster single chain clonotypes data (TRA or TRB)
- cluster paired chains clonotypes data (TRA and TRB)
- visualize similiarity of single chain or paired chains clonotypes with TSNE
- (In progress) train clustering ,odel and predict single chain or paired chains clonotypes data
- (In progress) train classificational model and predict single chain or paired chains clonotypes data

<img width="768" alt="image" src="https://github.com/yuliajk/tcr_emb/assets/74723905/5aeba670-2088-4c38-bfb7-21526a607c42">


# Instull
```
git clone https://github.com/yuliajk/tcr_emb.git
cd tcr_emb
conda create -n tcremb_env ipython python=3.11
conda activate tcremb_env   # or: "source activate tcremb_env" depending on your conda setup
pip install -r requirements.txt
```


# Input data
You can use as input data single chain or paired chains datasets:
- Rep-seq repertuares TRA or TRB data
- 10x VDJ data
- Clonotypes databases

## Input format
### Common requirements
**Single chain columns:**
- data_id - you can provide any column for identification. 
- cdr3aa - required
- v - required
- j - required
- chain - required
- label - either epitope either other labels for clustering or classification
  
**Double chain columns:**
- data_id - you can provide any column for identification. 
- a_cdr3aa - required
- TRAV - required
- TRAJ - required
- b_cdr3aa - required
- TRBV - required
- TRBJ - required
- label - either epitope either other labels for clustering or classification

**Other format requiremnts**
1. V and J genes must be provided in a format: TRAV35, TRBV11-2 (witout *01)
2. No missing data of any of the columns: v,j or cdr3
3. No ',' , '.' , ':' , ';' , '*' , '_' , '"' or other special simbols in V, J or CDR3

### Single chain example
| data_id | cdr3aa |v| j| chain | label |
| :---:   | :---: | :---: | :---: | :---: | :---: |
| GACTGCGCATCGTCGG-28   | CAGHTGNQFYF | TRAV35	 | TRAJ49 | TRA | IVTDFSVIK |
| GACTGCGCATCGTCGG-28   | CASSWGGGSHYGYTF | TRBV11-2  | TRBJ1-2 | TRB | IVTDFSVIK |

### Paired chains example
| data_id | a_cdr3aa | TRAV | TRAJ | b_cdr3aa | TRBV | TRBJ | labels |
| :---:   | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| GACTGCGCATCGTCGG-28   | CAGHTGNQFYF | TRAV35	 | TRAJ49 | CASSWGGGSHYGYTF | TRBV11-2  | TRBJ1-2 | IVTDFSVIK |

# Run tcremb

## From terminal
In tcr_emb repo directory
The command below will extrcat TRA clonotypes from my_input_data.txt, calculate distance scores and save clusters table with epitope as label in folder tcremb_outputs/my_run
```
python tcremb_run.py --input my_input_data.txt --runname my_run --chain TRA --label epitope
```

- **chain**
values: TRA, TRB, TRA_TRB
- **mode**
values: clstr, scores   In prgress: clsf, clstr_clsf, clstr_pre, clsf_pred, clstr_clsf_pred. clstr: will calculate dict scores, pca and perform clustering
scores: will only calculate dist scores

- **skip_scores**
if you have alreday calculted dist scored before and thay are in the aoutput folder (res_TRA.txt or res_TRB.txt), you can skip this stape by passing --skip_scores True
- **data_id**
you can define the column of input data that contains identificator of your data. This column will be included in output files
- **label**
if you provide the label for clustering, clustering metrics (purity, total cluster, cluster_matched) will be counted and labels of each cluster will be defined

## From Jupyter Notebook
Assign run_name and label for clustering. Load ypur data
```
run_name = 'my_run'
label = 'antigen.epitope'

data_preped = pd.read_csv('my_data.txt',sep='\t')

```

**Single chain pipeline** 

```
tcremb = TCRemb.TCRemb(run_name, data_preped)
tcremb.tcremb_clonotypes('TRA')
tcremb.tcremb_dists_count('TRA')
tcremb.tcremb_dists('TRA')
tcremb.tcremb_pca('TRA')
tcremb.tcremb_tsne('TRA')

tcremb.tcremb_clonotypes('TRB')
tcremb.tcremb_dists_count('TRB')
tcremb.tcremb_dists('TRB')
tcremb.tcremb_pca('TRB')
tcremb.tcremb_tsne('TRB')
```

**Data stored in object**
```
tcremb.input_data
tcremb.annot['TRA']
tcremb.clonotypes['TRA']
tcremb.pca['TRA']
tcremb.pca_clones['TRA']
tcremb.tsne['TRA']
tcremb.tsne_tsne['TRA']
```

**Paired chains pipeline** 
```
run_name = 'my_paired_run'
chain = 'TRA_TRB'
tcremb_paired = TCRemb.TCRemb(run_name, data_preped)
tcremb.tcremb_clonotypes(chain)
tcremb.tcremb_dists_count('TRA')
tcremb.tcremb_dists_count('TRB')
tcremb.tcremb_dists('TRA')
tcremb.tcremb_dists('TRB')
tcremb.tcremb_pca(chain)
tcremb.tcremb_tsne(chain)

```

**TSNE visualization**
Plot by all data instances with clonotypes duplicates. For example all cells or all rows from database
```
chain = 'TRA_TRB'
ml_utils.tsne_plot(pd.merge(tcremb.tsne[chain],tcremb.annot[chain]),  label, f'My_data, {chain}, TSNE colored by {label}')
```
Example
<img width="529" alt="image" src="https://github.com/yuliajk/tcr_emb/assets/74723905/296ab510-8308-4c54-8bc1-d73b600114c0">




Plot by clonotypes with clonotype size - coint of input instances with this clonotype.
```
chain = 'TRA_TRB'
ml_utils.tsne_plot(pd.merge(tcremb.tsne_clones[chain],tcremb.annot[chain]), label_high, f'My_data, {chain}, TSNE colored by {label}', to_size = 'clone_size')
```
Example
<img width="533" alt="image" src="https://github.com/yuliajk/tcr_emb/assets/74723905/6a11a2de-bd98-4d99-817d-5b1bf748a66f">



**Clustering**
```
kmeans = TCRemb.TCRemb_clustering('KMeans')
chain='TRA'
kmeans.clstr(chain,tcremb, label)
```
You can also pass ypur model with number of clusters
```
kmeans = TCRemb.TCRemb_clustering('KMeans')
chain='TRA'
n_clusters = 5295
random_state = 8
model =  KMeans(n_clusters=n_clusters, random_state=random_state)
kmeans.clstr(chain,tcremb, label, model)
```