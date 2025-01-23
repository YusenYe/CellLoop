# CellLoop: identifying single-cell 3D genome chromatin loops #

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
CellLoop algorithm that applies a density-based center detection algorithm to recognize chromatin loops of a specific single-cell based on approval of both intracell topology and intercell background strength of chromatin interactions. Intracell topology aims to extract the proximity of the genomic spatial positions within a single cell, while intercell background strength aims to extract the interaction probability of genomic positions of single cells within a specific biological context, such as cell states or spatial domains. Density-based center detection algorithm only needs to define two simple metrics for each interaction in parallel, making it possible to efficiently identify chromatin loops in thousands of single-cell 3D maps. Several advantages include:

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
conda create -n CellLoop python=3.8 -r requirements.txt
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

### Running CellLoop

- Initial parameters

We will take the analysis of the Dip-C dataset as an example to explain the meaning and default values of the parameters.
```
#########################Initial parameters##############################################
main_dir='/home/dell/Desktop/CellLoop_test/CellLoop'
DATASET='Dip-C'
BINSIZE=100e3
if BINSIZE<=10e3:
    MAXDIST=5000000
    knn_cell_num=25
else:
    MAXDIST=10000000
    knn_cell_num=25
LOW_CUTOFF=1e3   
GENOME='mm10'
#dataset dir
INDIR="/mnt/nas/Datasets/2020_Dip-C/GSE162511"##
FILE_SUFFIX=".contacts.pairs.txt.gz"
MINCONCNUM=250000
CHR_COLUMNS=[1,3]
POS_COLUMNS=[2,4]
feature='higashi_embed'
OUTDIR=INDIR+'_'+str(int(BINSIZE//1000))+'kb'+'_knn='+str(knn_cell_num)
```
**A.** main_dir：the directory of CellLoop package. **B.** DATASET:The dataset currently being analyzed. **C.** BINSIZE：bin size used for binning the contacts. **D.**  MAXDIST: maximum distance from diagonal to consider.  **E.**  knn_cell_num: maximum cell number of KNN graph.  **F.**  LOW_CUTOFF: cut-off for removing short-range contacts. **G.** GENOME: genome name; hgxx or mmxx.  **G.** GENOME: genome name; hgxx or mmxx. **H.** INDIR: the directory of the input single-cell 3D genome data. **I.** FILE_SUFFIX: suffix of the input files. **J.** MINCONCNUM: minimum contact number for single cells. **K.** CHR_COLUMNS and POS_COLUMN:two integer column numbers for chromosomes and two integer column numbers for read positions in singel-cell 3D genome file.  **L.** feature: the type of feature used for initializing the embedding of cells. **M.** OUTDIR:output directory. 
You can also set more parameters in the create_parser() function of CellLoop package. For example, use "--threaded" and "--num-proc" to set the number of processes used in threaded mode. This mode utilizes multiprocessing on a single machine. 

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
