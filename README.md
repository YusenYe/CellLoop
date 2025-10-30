# CellLoop:Identifying single-cell 3D genome chromatin loops #

----------

### Latest updates: January 15, 2025，version 0.0.1
## Contents：
1. [Overview](#Overviewn)
2. [Installation](#Installation)
3. [Usage](#Usage)
4. [Update Log](#UpdateLog)
5. [Maintainers](#Maintainers)
6. [Citation](#Citation)
## 1. Overview
![](/support/CellLoop.png)
CellLoop identifies chromatin loops in individual cells using a density-based center detection algorithm, which integrates intra-cellular and neighboring inter-cellular contact maps through a re-voting strategy (Fig. 1a). By capturing both local topological associations and contextual support from neighboring cells, CellLoop enables robust identification of chromatin loops in sparse single-cell contact maps. Fig. 1b illustrates the analytical workflow and potential applications of CellLoop, including the visualization and comparison of LFmap and CFmap, the identification and visualization of cell-type-specific chromatin loops, the characterization of chromatin loop features that define distinct cell states, and linking cell-type-specific loops with gene expression. Several advantages include:
(1). Loop Frequency Map (LFmap) obtained by CellLoop exhibited a significantly enhanced capacity. 

(2). CellLoop detected single-cell chromatin loops with an average time of merely 0.99s at 100kb resolution and ~13s at 10kb resolution in parallel with the support of 10 CPU processes, ensuring its capacity to process tens of thousands or even hundreds of thousands of single-cell data within a tolerable time. 

(3). CellLoop enables the signature of single cell from the perspective of chromatin loops. This, in turn, helps to facilitate the insight of subtle biological difference for single cells and further to reinterpret the biological functions of cell states or spatial domains from the perspective of 3D genomic chromatin loop-driven mechanisms.

## 2. Installation
please follow the following steps:
1. Install python version >= 3.8.
2. Make sure MPI (e.g. open-mpi) is installed on your system and is on your path (which mpicc can return a path).
3. Clone this repository and cd into it as below.
```
git clone https://github.com/YusenYe/CellLoop.git
cd CellLoop
```
4. Create a new conda environment and install the relevant packages according to the requirements.txt.
```
conda create -n CellLoop python=3.8
conda activate CellLoop
pip install -r requirements.txt
```
## 3. Usage
### Download the dataset that will be analyzed.
- Dip-C dataset: [Tan,2021](https://www.cell.com/cell/fulltext/S0092-8674(20)31754-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420317542%3Fshowall%3Dtrue) 
                 [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162511](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162511)
- HiRES dataset: [Liu,2023](https://www.science.org/doi/10.1126/science.adg3797)
                 [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223917](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223917)
- GAGE-seq dataset: [zhou,2024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238001)
                    [https://www.nature.com/articles/s41588-024-01745-3](https://www.nature.com/articles/s41588-024-01745-3)


### Required input files
- Single-cell 3D genome dataset: For the input of Dip-C data analysis, directly use the file with the suffix '.contacts.pairs.txt.gz' provided by the original data. For the input of HiRES data analysis, directly use the file with the suffix '.pairs.gz' provided by the original data. For the input of GAGE-seq, the data preprocessing process is provided in CellLoop package '/src/GAGE-seq_preprocess.py'. The generated single-cell files have a suffix of '.pairs'. All the single-cell data are placed in the same directory.
- Other initial features or embedding of single cells：**Single-modality 3D genomic data**: For Dip-C dataset, CellLoop applied Higashi, an integrative framework by the formation of hypergraph representation learning, to obtain initial cell embedding. The data preprocessing and operation of Higashi can be referred to [https://github.com/ma-compbio/Higashi](https://github.com/ma-compbio/Higashi).
                                                      **Double-modality with 3D genomic and other omics data**: CellLoop obtains the initial cell embeddings through the data of other omics by uniform manifold approximation and projection (UMAP). For HiRES dataset, the data preprocessing process is provided in CellLoop package '/src/Hires_preprocess.py'. For GAGE-seq dataset, the data preprocessing process is provided in CellLoop package '/src/GAGE-seq_preprocess.py'. The generated single-cell files is saved as 'adata_rna.h5ad'

### CellLoop Parameter Guide
#### 1. Required Parameters

Before running **CellLoop**, you must initialize several essential parameters according to your analysis task.

| **Parameter** | **Description** | **Example / Notes** |
|----------------|-----------------|---------------------|
| `INDIR` | Path to the directory containing the input files. | e.g., `INDIR="/mnt/nas/Datasets/2020_Dip-C/GSE162511"` |
| `FILE_SUFFIX` | Suffix of the raw input files. | e.g., `FILE_SUFFIX=".contacts.pairs.txt.gz"` |
| `CHR_COLUMNS` | Two integer column numbers indicating the chromosome columns in the raw file. | e.g.,`CHR_COLUMNS`=[1,3] |
| `POS_COLUMNS` | Two integer column numbers indicating the read position columns in the raw file. |e.g.,`POS_COLUMNS`=[2,4]  |
| `GENOME` | Genome name. | e.g., `hg38` or `mm10` |
| `MAXDIST` | Maximum genomic distance from the diagonal to consider. | e.g.,`MAXDIST`=10000000 |
| `MINCONCNUM` | Minimum contact number per single cell (used to filter out low-quality cells). |  e.g.,`MINCONCNUM`=250000 |
| `BINSIZE` | Bin size used for binning contacts; corresponds to the resolution of loop detection. | e.g., `10000` for 10 Kb |
| `LOW_CUTOFF` | Cutoff value to remove ultra-short-range contacts. |  e.g., '1000' for less than 1kb |
| `FEATURE` | Type of feature used for initializing the cell embeddings. | e.g., `FEATURE='higashi_embed'` |
| `KNN_CELL_NUM` | Maximum number of neighboring cells used to estimate interaction probabilities within a biological context. |  |
| `OUTDIR` | Path to the output directory where results will be saved. |  |
| `CHR_LEN` | Path to the chromosome length file. |  |

> 💡 **Note:**  
> Before setting parameters, you need to modify the `CellLoop.py` file to specify the location of the CellLoop package, for example:
> ```python
> main_dir = '/home/dell/Desktop/CellLoop_test/CellLoop'
> ```
> *(Path to the directory containing the `CellLoop.py` file.)*

---

#### 2. Optional Parameters

Additional parameters can be configured in the `create_parser()` function of `CellLoop.py` to control computational settings and local interaction definitions.

| **Parameter** | **Description** | **Notes** |
|----------------|-----------------|------------|
| `--threaded` and `-n` | Control multiprocessing. When `--threaded=True`, CellLoop utilizes multiprocessing on a single machine. The `-n` parameter specifies the number of parallel threads. The default is single-thread mode (`--threaded=False`). | |
| `--local-lower-limit` | Number of bins around the center (in each direction) to exclude from the local neighborhood. |  |
| `--local-upper-limit` | Number of bins around the center (in each direction) that define the local neighborhood. |  |
| `--filter file` | BED file specifying regions to be filtered out. |  |

>🧠 Tip
>
>For best performance:
>Use multi-core mode (--threaded True -n [CPU cores]) when analyzing large single-cell Hi-C datasets.
---

### Running CellLoop

#### Example: Dip-C Dataset Analysis
---
#####Overview

In this section, we take the **Dip-C dataset** as an example to demonstrate how to configure parameters and run **CellLoop** step-by-step.  
CellLoop is a scalable framework for identifying **chromatin loops at single-cell resolution**, integrating both intra-cell and inter-cell chromatin interaction signals.
---

- initial parameters

```
DATASET='Dip-C'
BINSIZE=100e3
MAXDIST=10000000
KNN_CELL_NUM=50
LOW_CUTOFF=1e3   
GENOME='mm10'
INDIR="/mnt/nas/Datasets/2020_Dip-C/GSE162511"##
FILE_SUFFIX=".contacts.pairs.txt.gz"
MINCONCNUM=250000
CHR_COLUMNS=[1,3]
POS_COLUMNS=[2,4]
FEATURE='higashi_embed'
OUTDIR=INDIR+'_'+str(int(BINSIZE//1000))+'kb'+'_knn='+str(KNN_CELL_NUM)
```
- Execute CellLoop

We can execute CellLoop in the following command-line way.
```
cd CellLoop
conda activate CellLoop 
python CellLoop.py
```
We can also execute CellLoop step1-4 by step in the IDE. CellLoop includes:

**Step 1.** Binning: we first constructed an undirected graph with bins as nodes from the raw contact map of the single cell with a specified resolution (20kb) with bin_sets() function.
```
bin_sets(args.indir, args.suffix, binsize = args.binsize, outdir = bin_dir, \
                      chr_columns = args.chr_columns, pos_columns = args.pos_columns, \
                      low_cutoff = args.low_cutoff, mincontactnum=args.mincontactnum,dataset=args.dataset,n_proc = n_proc, rank = rank, logger = logger)
```
**Step 2.** Generating the enhanced contact map of the single cells: we constructed a weighted undirected k-NN graph with cells as nodes by connecting the single cell to its k-NNs based on Euclidean distance calculated using the initial cell embedding space and generated the enhanced contact map of the single cells using the information of neighboring cells.
```
get_aggr_cells(bin_dir=bin_dir,outdir=SampleAggrCell_dir,feature_dir=feature_dir,feature=feature,chrom_lens=chrom_dict,\
                           n_proc=n_proc, rank=rank,binsize=args.binsize,k=knn_cell_num)
```
**Step 3.** Calling chromatin loops of single cells:   we defined an aggregating interaction strength matrix  based on approval of both intracell topology and intercell background strength of chromatin interactions. Next, CellLoop algorithm is used to distinguish chromatin loops from the aggregating interaction strength matrix by formulating our objective as a peak calling problem from a two-dimensional matrix.
```
scloops(indir = indir, outdir = scloops_outdir, chrom_lens = chrom_dict, \
                            binsize = args.binsize, dist = args.dist, \
                            neighborhood_limit_lower = args.local_lower_limit, \
                            neighborhood_limit_upper = args.local_upper_limit,alpha=alpha,\
                            rank = rank, n_proc = n_proc, max_mem = args.max_memory, logger = logger)
```
The scloops function can visually display results by specifying Figure = True. 

**Step 4.** Generating chromatin loop frequency map(LFmap): We further merge chromatin loops of single cells detected by CellLoop to generate LFmap. 
```
aggreloops(indir = scloops_outdir,outdir = aggreloops_dir,\
                       chrom_lens = chrom_dict, binsize = args.binsize,
                       rank = rank, n_proc = n_proc, logger = logger)
```
The aggreloops function can visually display results by specifying Figure = True. 

## 4. Update Log
## 5. Maintainers
Yusen Ye (ysye@xdiian.edu.cn)
## 6. Citation
