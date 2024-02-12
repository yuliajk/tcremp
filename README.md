# TCRemb
TCRemb is a pakage depeloped for T cell gene receptors comparison.
TCRemb is in the in active development right now.

With TCRemb you can:
- cluster single chain clonotypes data (TRA or TRB)
- cluster paired chains clonotypes data (TRA and TRB)
- visualize similiarity of single chain or paired chains clonotypes with TSNE
- (In progress) train clustering ,odel and predict single chain or paired chains clonotypes data
- (In progress) train classificational model and predict single chain or paired chains clonotypes data

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
