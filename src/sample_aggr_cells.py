# -*- coding: utf-8 -*-
"""
Created on Wed Jun 29 14:46:12 2022

@author: user
"""
import os
import pandas as pd
import scanpy as sc
import random
import numpy as np
from plot.plot_units import get_bin_matrix,get_bin_sc_matrix
from scipy import sparse
import re

#from sklearn.neighbors import NearestNeighbors


def computing_knn(adata,k=20):
    #adata=sc.AnnData(df.values,dtype=float)
    sc.tl.pca(adata, n_comps=20)
    sc.pp.neighbors(adata, n_neighbors=k,method='umap')
    return adata.obsp['connectivities']
    #show neibor
    #temp= adata_conprob.obsp['connectivities'][1000,:]
    #temp[temp>0]
    # sc.tl.tsne(adata_conprob,n_pcs=20)
    # sc.pl.tsne(adata_conprob, color=['leiden','cell-type cluster','age','tissue'])

   
def generate_aggr_mat(cellid,bin_dir,outdir_u,binfiles,using_filenames,weights,binsize,chrom_list,suffix=".bedpe"):
    #recent filename
    sc_filename=binfiles[cellid]
    setname=re.sub("\.*"+suffix,"",os.path.basename(sc_filename))
    if os.path.exists(os.path.join(outdir_u, ".".join([setname,chrom_list[-1][0],"sc",'npz']))):
        pass
    else:
        print(cellid,setname)
        edgelist=pd.DataFrame([])
        for filename,weight in zip(using_filenames,weights.tolist()[0]):
            print(weight)
            #filename=using_filenames[0];weight=weights[0]
            #setname=re.sub("\.*"+suffix,"",os.path.basename(filename))
            filename=os.path.join(bin_dir,filename)
            temp=pd.read_csv(filename,header=None,index_col=None)
            temp[6]=temp[6]*weight
            #print(filename,len(temp))
            edgelist=pd.concat([edgelist,temp],axis=0)
            edgelist=edgelist.sort_values(by=[0,1,4])
            print(len(edgelist))
            #random select one chrom
        
        #read sc edgelist
        sc_filename=os.path.join(bin_dir,sc_filename)
        sc_edgelist=pd.read_csv(sc_filename,header=None,index_col=None)
        
        for chrom,chrom_len in chrom_list:
            print(chrom,chrom_len)
            #chrom='chr2';chrom_len=chrom_lens[chrom]
            m_all=get_bin_matrix(edgelist,binsize,chrom,chrom_len).toarray()
            m_all=m_all/np.sum(weights)
            #show_mat(m,count=1)
            #m= np.triu(m, 1)
            m=get_bin_sc_matrix(sc_edgelist,binsize,chrom,chrom_len)
            
            # if Figure==True:
            #     from plot.plot_units import show_mat
            #     start=500;end=800
            #     #start=6800;end=7200
                
            #     sub_matrix=m_all[start:end,start:end]
            #     show_mat(sub_matrix+sub_matrix.T,count=1)
                
            if m is not None:
                m=m.toarray()
                allmatrix_sp_all=sparse.csr_matrix(m_all) 
                output_filename = os.path.join(outdir_u, ".".join([setname,chrom, "aggr",'npz']))
                sparse.save_npz(output_filename,allmatrix_sp_all)
                allmatrix_sp=sparse.csr_matrix(m) 
                output_filename = os.path.join(outdir_u, ".".join([setname ,chrom, "sc",'npz']))
                sparse.save_npz(output_filename,allmatrix_sp)
            else:
                print(filename,chrom)
                pass
                
                

def get_proc_cellids(cellids, rank, n_proc):
    indices = list(range(rank,len(cellids), n_proc))
    proc_cellids = [cellids[i] for i in indices]
    return proc_cellids
        
        
        
def get_aggr_cells(bin_dir,outdir,feature_dir,feature,chrom_lens,n_proc,rank,binsize,k=20,suffix=".bedpe",Figure=False):
    """
    Sampling aggregated cells for candidated loops
    
    Parameters:
        bin_dir: bin pairs location
        CPD_dir: contact distribution loation
        chrom_lens: chrom lengths
    """
    
    if (feature=='higashi_embed'):
        adata=sc.read(feature_dir+'/adata_higashi.h5ad')
        sc.pp.filter_genes(adata,min_cells=np.floor(adata.shape[0]*0.01))
        sc.pp.normalize_total(adata, target_sum=adata.shape[0])
        Knn_mat=computing_knn(adata,k).todense()
        cellids=range(adata.shape[0])
        outdir_u=outdir+'_'+feature
    
    elif feature=='RNA_embed':
        adata=sc.read(feature_dir+'/adata_rna.h5ad')
        #sc.pp.filter_genes(adata,min_cells=np.floor(adata.shape[0]*0.01))
        #sc.pp.normalize_total(adata, target_sum=adata.shape[0])
        Knn_mat=computing_knn(adata,k).todense()
        cellids=range(adata.shape[0])
        outdir_u=outdir+'_'+feature    
    
    # if Figure==True:
    #     from plot.plot_units import show_heatmap,generate_cmap,white1,red
    #     show_heatmap(Knn_mat[100:300,:][:,100:300],generate_cmap([white1,red]))
    # sc.tl.umap(adata)
    # sc.tl.leiden(adata)
    # sc.pl.umap(adata,color=['leiden'],size=30,wspace=0.5)
    #import pandas as pd
    #for filelist
    #cellnames=adata.obs['cell'].values
    cellnames=pd.read_csv(feature_dir+'/filelist.txt',header=None,index_col=None,sep='\t')
    # cellnames.to_csv(feature_dir+'/filelist.txt',header=None,index=None)
    # cellnames=[str(line[0]).split('/')[-1].rstrip('.txt') for line in cellnames.values]
    # pd.DataFrame(cellnames).to_csv(feature_dir+'/filelist.txt',header=None,index=None)     
    try:
        os.makedirs(outdir_u)
    except:
        #print(e, 'excepting error in binning')
        pass    
    np.fill_diagonal(Knn_mat, 1)
    
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    
    
    #binfiles=os.listdir(bin_dir)
    #binfiles.sort() 
    binfiles=[cell+suffix  for cell in cellnames[0].tolist()]
    # if len(binfiles)!=len(cellnames):
    #     print("error!")
    #     #os.exit()
    #random sample number
    #cellids=random.sample(range(all_contact_p.shape[0]),Sample_cell_n)
    #cellids.sort()
    #cellids=range(using_compartments.shape[0])
    proc_cellids=get_proc_cellids(cellids,rank,n_proc)
    #proc_cellids=get_proc_cellids(cellids, 0, 1)
    #print(cellids)
    
    consider_cell_number=[]
    for cellid in proc_cellids:
        #cellid=proc_cellids[0]
        using_fids=np.where(Knn_mat[cellid,:]>0.1)[1]
        print(len(using_fids))
        consider_cell_number.append(len(using_fids))
        
        weights=Knn_mat[cellid,using_fids]
        #print(weights)
        using_filenames=[binfiles[fid] for fid in using_fids]
        #using_filename=binfiles[cellid]
        #if len(using_filenames)>min_cell_n:
        generate_aggr_mat(cellid,bin_dir,outdir_u,binfiles,using_filenames,weights,binsize,chrom_list)
        #generate corresponding single cells
        #using_fid=cellid
        
        
    
    
    # #computing knn graph
    # cell_corr=all_contact_p.T.corr()
    # #visual corr clusters
    # adata_cell_corr=sc.AnnData(cell_corr,dtype=float)
    # #cell target 1
    # Cell_label_1=DipC_CellTypes_F['cell-type cluster']
    # groups=Cell_label_1.values
    # adata_cell_corr.obs['cell-type cluster']=pd.Categorical(
    #     values=groups.astype('U'),
    #     categories=natsorted(map(str, np.unique(groups))),
    # )
    # #sc.tl.pca(adata_conprob, svd_solver='arpack')
    
    # sc.tl.pca(adata_cell_corr, n_comps=20)
    # sc.tl.tsne(adata_cell_corr,n_pcs=20)
    # sc.pl.tsne(adata_cell_corr, color='cell-type cluster')
    
    # sc.tl.umap(adata_cell_corr)
    # sc.pl.umap(adata_cell_corr, color='cell-type cluster')

def generate_bulk_mat(bin_dir,bulk_dir,binfiles,binsize,chrom_list,suffix=".bedpe",Figure=False):
    
    #recent filename
    #sc_filename=binfiles
    
    
    #edgelist=pd.DataFrame([])
    m_all=dict()
    for k,filename in enumerate(binfiles):
        #filename=binfiles[0]
        filename=os.path.join(bin_dir,filename)
        temp=pd.read_csv(filename,header=None,index_col=None)
        #print(filename,len(temp))
        #edgelist=pd.concat([edgelist,temp],axis=0)
        #edgelist=edgelist.sort_values(by=[0,1,4])
        print(k,len(temp))
        #random select one chrom
        for chrom,chrom_len in chrom_list:
            #print(chrom,chrom_len)
            #chrom='chr2';chrom_len=chrom_lens[chrom]
            m=get_bin_matrix(temp,binsize,chrom,chrom_len).toarray()
            NUM = int(np.ceil(chrom_len / binsize))
            add_n=NUM-m.shape[0]
            if add_n>0:
                m=np.pad(m,((0,add_n),(0,add_n)),'constant',constant_values = (0,0))
            if k==0:
                m_all[chrom]=m
            else:
                m_all[chrom]+=m
    
    for chrom,chrom_len in chrom_list:
        m=m_all[chrom]
        if Figure==True:
            start=1550;end=1750
            sub_matrix=m[start:end,start:end]
            from plot.plot_units import show_mat
            show_mat(sub_matrix+sub_matrix.T,count=1)
        allmatrix_sp_all=sparse.csr_matrix(m_all[chrom]) 
        output_filename = os.path.join(bulk_dir, ".".join([chrom, "bulk",'npz']))
        sparse.save_npz(output_filename,allmatrix_sp_all)
        
        


#(bin_dir,bulk_dir,binsize = args.binsize,chrom_lens=chrom_dict,suffix=args.suffix) 

#bin_dir,bulk_dir,binfiles,binsize,chrom_list,suffix=".bedpe",Figure=False       
def generate_bulk_hic(bin_dir,bulk_dir,binsize,chrom_lens,suffix=".bedpe"):
    try:
        os.makedirs(bulk_dir)
    except:
        #print(e, 'excepting error in binning')
        pass    
    binfiles=os.listdir(bin_dir)
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    #binfiles.sort()
    generate_bulk_mat(bin_dir,bulk_dir,binfiles,binsize,chrom_list,suffix)
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    

   



