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
4. Create a new conda environment using requirements.txt
```
conda create -n CellLoop python=3.8
```

6. 


## 3. Usage
## 4. Update Log
## 5. Maintainers
## 6. Citation
