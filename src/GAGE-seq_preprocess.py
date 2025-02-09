#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 17 00:56:51 2024

@author: dell
"""
import pandas as pd
import gzip



#########path of the downloaded data
path2dir='/mnt/data/Datasets/GAGE-seq_2024/'
###### single cell 3D data after preprocessing
sc3d2dir='/mnt/data/Datasets/GAGE-seq/SingleCells'
##### single cell Rna-seq data after preprocessing (CellLoop package dir)
main_dir='/home/dell/Desktop/CellLoop_test/CellLoop'
scRNA2dir=main_dir+'/Dataset/HiRES/CellFeatures'

samples=['GSM7657702_contact_mBC-HiC-0814-1.pairs','GSM7657700_contact_mBC-HiC-0716-2.pairs','GSM7657698_contact_mBC-HiC-0716-1.pairs']
RNA_samples=['GSM7657701_RNA_mBC-RNA-0814-1.tsv','GSM7657699_RNA_mBC-RNA-0716-2.tsv','GSM7657697_RNA_mBC-RNA-0716-1.tsv']
# #metacell dir
metacell_dir=path2dir+'/41588_2024_1745_MOESM3_ESM.xlsx'
clust_dat=pd.read_excel(metacell_dir,sheet_name='Supplementary Table 5')
max_open_files=1000
def oppf(filename, mode='r'):
    if filename.endswith('.gz'):
        return gzip.open(filename, mode) 
    else:
        return open(filename, mode)

# #sample=samples[1]
# open_files = {}
###################### preprecessing GAGE-seq #####################################
try:
    for s_n,sample in enumerate(samples):
        filedir=path2dir+sample
        pair_temp=pd.read_csv(filedir,sep='\t',header=None)
        lib_sub=sample.split('_')[2].rstrip('.pairs').split('-')
        lib=lib_sub[0]+'-'+lib_sub[2]+'-'+lib_sub[3]
        pair_temp=pair_temp.dropna()
        cellnames=pair_temp[4].unique()
        pair_temp[0]=[chrom.lstrip('mm10_') for chrom in pair_temp[0].values]
        pair_temp[2]=[chrom.lstrip('mm10_') for chrom in pair_temp[2].values]
        #cellname=cellnames[0]
        for c_n,cellname in enumerate(cellnames):
            print(s_n,c_n,cellname)
            temp=pair_temp.loc[pair_temp[4]==cellname]
            temp.to_csv(sc3d2dir+'/'+lib+'_'+cellname+'.csv',header=None,index=None)
except:
    
    clust_dat["uniq_barcode"] = clust_dat.apply(lambda row: f"{row['library']}:{row['well id']}", axis=1)
    
    clust_dat_dict = pd.Series(clust_dat["uniq_barcode"].values, index=clust_dat["uniq_barcode"]).to_dict()
    outf_dict = dict()
    open_files = {}

    for cls_id in pd.unique(clust_dat["uniq_barcode"]):
        outf_dict[cls_id] = "".join((sc3d2dir, "/", str(cls_id), ".pairs"))
    
    for s_n,sample in enumerate(samples):
        open_files = {}
        samplef = path2dir+sample
        lib_sub=sample.split('_')[2].rstrip('.pairs').split('-')
        lib=lib_sub[0]+'-'+lib_sub[2]+'-'+lib_sub[3]
        
        with oppf(samplef, 'rt') as infile:
            for d,dline in enumerate(infile):
                
                fields = dline.strip("\n").split("\t")
                
                cid = ":".join((lib, fields[-1]))
                
                if cid in clust_dat_dict:
                    fields[0]=fields[0].lstrip('mm10_')
                    fields[2]=fields[2].lstrip('mm10_')
                    
                    fields[-1] = cid
                    
                    wline = '\t'.join((fields))
                    cls = clust_dat_dict[cid]
                    if cls not in open_files:
                        open_files[cls] = open(outf_dict[cls], "a")
                    open_files[cls].write(wline + "\n")
                    # ofile = outf_dict[cls]
                    # with open(ofile, "a") as p:
                    #     p.write(wline + "\n")
                if len(open_files) >= max_open_files:
                    for file in open_files.values():
                        file.close()
                    open_files.clear()



######################preprocessing RNA-seq###########################################
#s_n=0;sample=RNA_samples[s_n]
all_cellgenes=pd.DataFrame()
for s_n,sample in enumerate(RNA_samples):
    filedir=path2dir+sample
    pair_temp=pd.read_csv(filedir,sep='\t',header=None)
    lib_sub=sample.split('_')[2].rstrip('.tsv').split('-')
    lib=lib_sub[0]+'-'+lib_sub[2]+'-'+lib_sub[3]
    pair_temp[1]=[lib+':'+cellname   for cellname in pair_temp[1].values]
    pair_temp[0]=[chrom.lstrip('mm10_') for chrom in pair_temp[0].values]
    if s_n==0:
        all_cellgenes=pair_temp
        
    else:
        all_cellgenes=pd.concat([all_cellgenes,pair_temp],axis=0)
        
all_cellgenes['val']=1 
cellgene_counts=all_cellgenes.groupby([1,2])['val'].count().reset_index()
cellgene_counts=cellgene_counts.sort_values(by=1)


from scipy.sparse import coo_matrix
import numpy as np
cellnames=cellgene_counts[1].unique().tolist()
genenames=cellgene_counts[2].unique().tolist()

cellgene_counts.columns=['cell','gene','val']

cellnames_pd=pd.DataFrame(cellnames)
cellnames_pd.columns=['cell']
cellnames_pd['cell_loc']=cellnames_pd.index

genenames_pd=pd.DataFrame(genenames)
genenames_pd.columns=['gene']
genenames_pd['gene_loc']=genenames_pd.index

cellgene_counts=pd.merge(cellgene_counts,cellnames_pd,on='cell',how='left')
cellgene_counts=pd.merge(cellgene_counts,genenames_pd,on='gene',how='left')
# cellgene_counts['gene_loc']=[genenames.index(gene) for gene in cellgene_counts[2].values]
# cellgene_counts['cellname_loc']=[cellnames.index(cellname) for cellname in cellgene_counts[1].values]
data=cellgene_counts['val'].values; row=cellgene_counts['gene_loc'].values.astype(np.int32); col=cellgene_counts['cell_loc'].values.astype(np.int32)
m= coo_matrix((data,(row,col)),shape=(max(row)+1,max(col)+1))
Ndata=m.toarray()

  
import scanpy as sc
data_pd=pd.DataFrame(Ndata.T, index=cellnames, columns=genenames)
adata_sc=sc.AnnData(data_pd.values,obs=cellnames_pd,var=genenames_pd,dtype=float)

sc.pp.filter_cells(adata_sc, min_genes=200)
sc.pp.filter_genes(adata_sc, min_cells=50)
adata_sc.layers["counts"] = adata_sc.X.copy()  
sc.tl.pca(adata_sc, svd_solver="arpack")
sc.pp.neighbors(adata_sc, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata_sc)
sc.tl.leiden(adata_sc)
sc.pl.umap(adata_sc,color=['n_genes','leiden'], palette='tab20',size=10,wspace=1,hspace=0.1)

#adata_sc.write(scRNA2dir+'/RNA_embed/adata_rna_raw.h5ad')

# Normalizing to median total counts
sc.pp.normalize_total(adata_sc, target_sum=1e4)
# Logarithmize the data
sc.pp.log1p(adata_sc)  

sc.pp.highly_variable_genes(adata_sc, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_sc.raw = adata_sc.copy()
adata_sc = adata_sc[:, adata_sc.var.highly_variable]
# sc.pp.scale(adata_sc, max_value=10)
sc.tl.pca(adata_sc, svd_solver="arpack")
# sc.pl.pca(adata_sc)
#sc.pl.highly_variable_genes(adata_sc)
# sc.tl.pca(adata_sc)
sc.pp.neighbors(adata_sc, n_neighbors=15, n_pcs=20)
sc.tl.umap(adata_sc)
sc.tl.leiden(adata_sc)
sc.pl.umap(adata_sc,color=['n_genes','leiden'], palette='tab20',size=10,wspace=1,hspace=0.1)
adata_sc.write(scRNA2dir+'/RNA_embed/adata_rna.h5ad')



# adata_sc=sc.read('/mnt/data/Datasets/GAGE-seq/CellFeatures/adata.h5ad')
# filelist_dir='/mnt/data/Datasets/GAGE-seq/SingleCells/filelist.txt'
# filenames=pd.read_csv(filelist_dir, index_col=None,header=None,sep='\t')
# adata_sc=adata_sc[[cellname in filenames[0].values for cellname in adata_sc.obs_names],:]


