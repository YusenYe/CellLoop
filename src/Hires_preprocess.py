#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 15:02:42 2023

@author: dell
"""
import pandas as pd 
import scanpy as sc



######## package dir##################
main_dir='/home/dell/Desktop/CellLoop_test/CellLoop'
#########path of the downloaded data
path2dir='/mnt/nas/Datasets/2023_Science'











###########################################################
RNAData_dir=path2dir+'/GSE223917_HiRES_emb.rna.umicount.tsv'
#read RNA-seq expression
R_Data=pd.read_csv(RNAData_dir,sep='\t',header=0,index_col=0)
#read cell info
cellinfo_dir=path2dir+'/GSE223917_HiRES_emb_metadata.xlsx'
cellinfo=pd.read_excel(cellinfo_dir,header=0,index_col=None)

#CellLoop package dir
using_cell_dir=main_dir+'/Dataset/HiRES/CellFeatures/RNA_embed/filelist.txt'
cellnames=pd.read_csv(using_cell_dir,header=None,index_col=None)

cellnames=[name.split('_',1)[1] for name in cellnames[0]]
cellinfo_2=pd.merge(pd.DataFrame(cellnames,columns=['Cellname']),cellinfo,how='left',on='Cellname')
cellinfo_2.index=cellinfo_2['Cellname']

R_Data=R_Data.loc[:,cellnames]

Markinfo=pd.DataFrame()
Markinfo['geneid']=R_Data.index.tolist()
Markinfo.index=R_Data.index.tolist()

adata3=sc.AnnData(R_Data.T,obs=cellinfo_2,var=Markinfo)

sc.pp.filter_genes(adata3, min_cells=5)
#adata3.raw=adata3
#sc.pl.highest_expr_genes(adata3, n_top=20, )
#sc.pp.filter_cells(adata3, min_genes=100)
sc.pp.normalize_total(adata3, target_sum=1e4)
sc.pp.log1p(adata3)

# to get gene promoter location
feature_dir=path2dir+'/CellFeatures/RNA_embed'
#adata3.write(feature_dir+'/adata_rna_raw.h5ad')

sc.pp.highly_variable_genes(adata3, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata3)

adata3.raw = adata3

adata3 = adata3[:, adata3.var.highly_variable]
#sc.pp.scale(adata3, max_value=10)
sc.tl.pca(adata3, svd_solver='arpack')

#sc.pl.pca(adata3, color='Son')

sc.pl.pca_variance_ratio(adata3, log=True)


sc.pp.neighbors(adata3, n_neighbors=15, n_pcs=20)
sc.tl.umap(adata3)
#sc.tl.tsne(adata,n_pcs=40)
sc.tl.leiden(adata3)
sc.pl.umap(adata3,color=['Celltype','leiden','Cellcycle phase'],palette='tab20',size=20,wspace=0.5)   


diff_analysis=True
if diff_analysis==True:
    #fig_add='/mnt/nas/Projects/Project18/SlidesOrFigures/Supplementary/Hierachical_Tree/'
    #read adata
    #adata=sc.read(fig_add+'adata.h5ad')
    # sc.tl.rank_genes_groups(adata,'leiden',method='t-test_overestim_var')
    # sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
    #adata.X[adata.X>0]=1
    ############case: speratly analyze subcelltypes:'Hippocampal Granule Cell'
    #######step 1:computing cell type-specific loops, first level analysis
    sc.tl.rank_genes_groups(adata3,'Celltype',method='t-test',use_raw=True)
    sc.pl.rank_genes_groups(adata3,n_genes=25,sharey=False)
    sc.pl.rank_genes_groups_heatmap(adata3,n_genes=100,cmap='bwr')
   
adata3.write(feature_dir+'/adata_rna.h5ad')
