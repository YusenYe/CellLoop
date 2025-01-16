# CellLoop: identifying single-cell 3D genome chromatin loops #

----------

### Latest updates: January 15, 2025，version 0.0.1
## Contents：
1. [Overview](#1)
2. [Installation](#2)
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
- Dip-C dataset:[Tan,2021](https://www.cell.com/cell/fulltext/S0092-8674(20)31754-2?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867420317542%3Fshowall%3Dtrue) 
                [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162511](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162511)
- HiRES dataset:[Liu,2023](https://www.science.org/doi/10.1126/science.adg3797)
                [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223917](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE223917)
- GAGE-seq dataset:[zhou,2024](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE238001)
                   [https://www.nature.com/articles/s41588-024-01745-3](https://www.nature.com/articles/s41588-024-01745-3)
- DropleHiC dataset:[Chang,2024](https://www.nature.com/articles/s41587-024-02447-1)
                    [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253407](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253407)
### required input files
- Single-cell 3D genome dataset: For the input of Dip-C data analysis, directly use the file with the suffix '.contacts.pairs.txt.gz' provided by the original data. For the input of HiRES data analysis, directly use the file with the suffix '.pairs.gz' provided by the original data. For the input of GAGE-seq and DropletHiC, The data preprocessing process is provided in CellLoop package '/src/GAGE-seq_preprocess.py' and '/src/DropletHiC_processing'. The generated single-cell files have a suffix of '.pairs'.
- Other initial features or embedding of single cells：





## 4. Update Log
## 5. Maintainers
## 6. Citation
