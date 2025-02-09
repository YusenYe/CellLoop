# -*- coding: utf-8 -*-
"""
Created on Thu Jun 23 15:39:29 2022

@author: user
"""


import os
import random
import numpy as np
#import networkx as nx
#import scipy as sp
from colormap import Color, Colormap
#import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
#import umap
import scanpy as sc
from natsort import natsorted



def color(value):
  digit = list(map(str, range(10))) + list("ABCDEF")
  if isinstance(value, tuple):
    string = '#'
    for i in value:
      a1 = i // 16
      a2 = i % 16
      string += digit[a1] + digit[a2]
    return string
  elif isinstance(value, str):
    a1 = digit.index(value[1]) * 16 + digit.index(value[2])
    a2 = digit.index(value[3]) * 16 + digit.index(value[4])
    a3 = digit.index(value[5]) * 16 + digit.index(value[6])
    return (a1, a2, a3)

#single color
red=color((207,0,25));blue=color((0,42,111)); white1=color((255,255,255));white2=color((200,200,200));white3=color((230,230,230));black=color((0,0,0));
red2=color((150,38,38))
red3=color((99,27,27))

red1=color((204,57,57));blue1=color((19,78,119));purple1=color((184,137,219))
color_3=[red1,blue1,purple1]





#nature review article
d_red_n=color((145,12,28));red_n=color((203,33,50));blue_n=color((37,133,198));
green_n=color((0,173,174));purple_n=color((164,73,150));
brown_n=color((122,76,61));yellow_n=color((255,137,60));grey_n=color((141,152,175))
colors_nature=[d_red_n,red_n,blue_n,green_n,brown_n,yellow_n]

colors_diff=[d_red_n,red_n,blue_n,purple_n,brown_n,grey_n]

color_multicmap=[white3,white2,grey_n,yellow_n,red,black]

color_mult=['#B237B1','#B54273','#A156A2','#D6AEEC','#E290B6','#BE95E3',
                '#854DBC','#E08DD3','#386588','#855C19','#E77F88','#D0A62F',
                '#89553C','#DF8165','#A44B1C','#A78425','#B5731F','#B36F6C',
                '#E4AB97','#E5AA61','#E7AB63','#92D399','#E7AB63','#3D7647',
                '#7BDBAD','#57A359','#4A85BF','#6294AF','#4E9D81','#3CA122','#317929']
                
                
                




def generate_cmap(colors):
    #colors=[red1,red1_2,white1,blue1_2,blue1]
    red_list=list()
    green_list=list()
    blue_list=list()
    #['darkblue','seagreen','yellow','gold','coral','hotpink','red'],['white','green','blue','red']
    for color in colors:
        col=Color(color).rgb
        red_list.append(col[0])
        green_list.append(col[1])
        blue_list.append(col[2])
    c = Colormap()
    d=  {   'blue': blue_list,
            'green':green_list,
            'red':red_list}
    mycmap = c.cmap(d) 
    return mycmap 

# def get_matrix_from_edgelist(edgelist):
#     g = nx.from_pandas_edgelist(edgelist, source = 'x1', target = 'y1', edge_attr = ['weight'], create_using = nx.Graph())
#     m = sp.sparse.csc_matrix(nx.adjacency_matrix(g).astype(float))
#     del g
#     return m



def bin_matrix(df, bin_size):
    df.loc[:,'x1'] = (df.loc[:,'x1'] // bin_size).astype(int)
    df.loc[:,'y1'] = (df.loc[:,'y1'] // bin_size).astype(int)
    return df


def get_bin_matrix(edgelist, bin_size,chrom,chrom_len):
    
    edgelist.columns=['chr1','x1','x2','chr2','y1','y2','weight']
    edgelist = edgelist[(edgelist['chr1'] == chrom) & (edgelist['chr1'] == edgelist['chr2'])]
    edgelist = bin_matrix(edgelist, bin_size)
    # NUM = int(np.ceil(chrom_len / binsize))
    # #print('NUM', NUM)
    # edges = pd.DataFrame({'x1':list(range(0, NUM-1)), 'y1':list(range(1, NUM))})
    # edges = pd.concat([edges, edgelist[['x1', 'y1']]], axis = 0)
    #edgelist.loc[:,'weight'] = 1
    edges=edgelist.groupby(['chr1','x1','y1'])['weight'].sum()
    edges=edges.reset_index()
    edges=edges[['x1','y1','weight']]
    #edges=edges[edges['y1']-edges['x1']>=1]
    #print(len(edges))
    #del edgelist
    #m = get_matrix_from_edgelist(edges)
    from scipy.sparse import coo_matrix
    data=edges['weight'].values; row=edges['x1'].values.astype(np.int32); col=edges['y1'].values.astype(np.int32)
    m= coo_matrix((data,(row,col)),shape=(max(col)+1,max(col)+1))
    return m


#edgelist=sc_edgelist;bin_size=binsize
def get_bin_sc_matrix(edgelist, bin_size,chrom,chrom_len):
    
    edgelist.columns=['chr1','x1','x2','chr2','y1','y2','weight']
    edgelist = edgelist[(edgelist['chr1'] == chrom) & (edgelist['chr1'] == edgelist['chr2'])]
    edgelist = bin_matrix(edgelist, bin_size)
    # NUM = int(np.ceil(chrom_len / binsize))
    # #print('NUM', NUM)
    # edges = pd.DataFrame({'x1':list(range(0, NUM-1)), 'y1':list(range(1, NUM))})
    # edges = pd.concat([edges, edgelist[['x1', 'y1']]], axis = 0)
    #edgelist.loc[:,'weight'] = 1
    #edges=edgelist.groupby(['chr1','x1','y1'])['weight'].sum()
    #edges=edges.reset_index()
    edges=edgelist[['x1','y1','weight']]
    #edges=edges[edges['y1']-edges['x1']>=1]
    #print(len(edges))
    #del edgelist
    #m = get_matrix_from_edgelist(edges)
    from scipy.sparse import coo_matrix
    data=edges['weight'].values; row=edges['x1'].values.astype(np.int32); col=edges['y1'].values.astype(np.int32)
    
    if len(edges)!=0:
        m= coo_matrix((data,(row,col)),shape=(max(col)+1,max(col)+1))
        return m
    else:
        return None




def log2_norm(matrix):
	return np.log2(1+np.abs(matrix)) * np.sign(matrix)

def log10_norm(matrix):
	return np.log10(1+np.abs(matrix)) * np.sign(matrix)

def plot_scatter(X1,X2):
    fig = plt.figure(figsize=(6, 6))
    plt.scatter(X1,X2,color=red,marker='*',s=3)
    plt.show()


def show_mat(matrix,mycmap=generate_cmap([white1,red]),count=0,colorbar=0):
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    if count==1:
        c=ax.matshow(log10_norm(matrix), cmap=mycmap)
    else:
        c=ax.matshow(matrix, cmap=mycmap)
    
    if colorbar==1:
        fig.colorbar(c,ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    return ax


def show_mat_withax(matrix,ax,mycmap=generate_cmap([white1,red]),count=0,colorbar=0):
    if count==1:
        c=ax.matshow(log10_norm(matrix), cmap=mycmap)
    else:
        c=ax.matshow(matrix, cmap=mycmap)
    
    if colorbar==1:
        plt.colorbar(c,ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    return ax
    


def show_heatmap(heatmap,mycmap=generate_cmap([blue,white1,red]),count=0):
    
    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    #mycmap=generate_cmap([blue,white1,red])
    if count==1:
        ax.matshow(log10_norm(heatmap), cmap=mycmap)
    else:
        ax.matshow(heatmap, cmap=mycmap)
    #fig.colorbar(im,ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])



def show_compartment(compartments,calib):
    x=range(len(compartments))
    fig = plt.figure(figsize=(8,3))
    ax = fig.add_subplot(211)
    ax.plot(x,compartments)
    ax = fig.add_subplot(212)
    ax.plot(x,calib)
    
        
def visual_random_bin_maps(bin_dir,chrom_dict,bin_size):
    """
    Parameters
    ----------
    bin_dir : address of bin matrix.

    """
    binfiles=os.listdir(bin_dir)
    filename=os.path.join(bin_dir, binfiles[random.randint(0,len(binfiles))])
    edgelist=pd.read_csv(filename,header=None,index_col=None)
    
    print(len(edgelist))
    #random select one chrom
    chrom_list = [(k, chrom_dict[k]) for k in list(chrom_dict.keys())];chrom=chrom_list[1][0]
    #chrom=chrom_names[1]
    m=get_bin_matrix(edgelist, bin_size , chrom = chrom, chrom_len = chrom_dict[chrom])
    #m= triu(m,1, format='csc')
    #show_mat(m.todense()[500:1000,500:1000])
    show_mat(m.todense()[500:1000,500:1000],count=1)


def visual_random_multicell_bin_maps(bin_dir,chrom_dict,bin_size,k=10):
    """
    Parameters
    ----------
    bin_dir : address of bin matrix.

    """
    binfiles=os.listdir(bin_dir)
    edgelist=pd.DataFrame()
    rand=random.randint(0,len(binfiles)-k)
    for i in range(k):
        filename=os.path.join(bin_dir,binfiles[rand+i])
        temp=pd.read_csv(filename,header=None,index_col=None)
        print(filename,len(temp))
        edgelist=pd.concat([edgelist,temp],axis=0)
    edgelist=edgelist.sort_values(by=[0,1,4])
    print(len(edgelist))
    #random select one chrom
    chrom_list = [(k, chrom_dict[k]) for k in list(chrom_dict.keys())];chrom=chrom_list[1][0]
    #chrom=chrom_names[1]
    m=get_bin_matrix(edgelist,bin_size,chrom = chrom,chrom_len = chrom_dict[chrom])
    #m= triu(m,1, format='csc')
    #m[m>k/3]=k/3
    show_mat(m.todense()[500:1000,500:1000],count=1)


#matrix=sub_matrix+sub_matrix.T
def show_mat_loops(matrix,start,end,upper_limit,temp,mycmap=generate_cmap([white1,red]),count=1):
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    #mycmap=generate_cmap([white1,red])
    if count==1:
        ax.matshow( log10_norm(matrix), cmap=mycmap)
    else:
        ax.matshow(matrix, cmap=mycmap)
    #fig.colorbar(im,ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    #heatmap_2D(submatrix,mycmap)
    #plot loops
    for index in temp.index:
        loc_x=temp.loc[index,'i']-start
        loc_y=temp.loc[index,'j']-start
        print(loc_x,loc_y)
        ax.plot(loc_x,loc_y,color='g',marker='s', markerfacecolor='none',markersize=2.5)
        ax.plot(loc_y,loc_x,color='g',marker='s', markerfacecolor='none',markersize=2.5)
        
def show_mat_loops_v2(matrix,start,temp,mycmap=generate_cmap([white1,red]),count=1):
    fig = plt.figure(figsize=(4, 4))
    ax = fig.add_subplot(111)
    #mycmap=generate_cmap([white1,red])
    if count==1:
        ax.matshow( log10_norm(matrix), cmap=mycmap)
    else:
        ax.matshow(matrix, cmap=mycmap)
    #fig.colorbar(im,ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    #heatmap_2D(submatrix,mycmap)
    #plot loops
    
    for index in temp.index:
        loc_x=temp.loc[index,'i']-start
        loc_y=temp.loc[index,'j']-start
        print(loc_x,loc_y)
        ax.plot(loc_x,loc_y,color='g',marker='s', markerfacecolor='none',markersize=2.5)
        ax.plot(loc_y,loc_x,color='g',marker='s', markerfacecolor='none',markersize=2.5)
        

def show_mat_loops_v3(matrix,start,temp,temp2,mycmap=generate_cmap([white1,red]),count=1,colorbar=0):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)
    #mycmap=generate_cmap([white1,red])
    if count==1:
        im=ax.matshow( log10_norm(matrix), cmap=mycmap)
    else:
        im=ax.matshow(matrix, cmap=mycmap)
    #fig.colorbar(im,ax=ax)
    ax.set_xticks([])
    ax.set_yticks([])
    #heatmap_2D(submatrix,mycmap)
    #plot loops
    for index in temp.index:
        loc_x=temp.loc[index,'i']-start
        loc_y=temp.loc[index,'j']-start
        print(loc_x,loc_y)
        ax.plot(loc_x,loc_y,color='g',marker='s', markerfacecolor='none',markersize=2,linewidth=1,alpha=0.7)
        ax.plot(loc_y,loc_x,color='g',marker='s', markerfacecolor='none',markersize=2,linewidth=1,alpha=0.7)
    
    for index in temp2.index:
        loc_x=temp2.loc[index,'i']-start
        loc_y=temp2.loc[index,'j']-start
        print(loc_x,loc_y)
        ax.plot(loc_x,loc_y,color='b',marker='s', markerfacecolor='none',markersize=2,linewidth=1,alpha=0.7)
        ax.plot(loc_y,loc_x,color='b',marker='s', markerfacecolor='none',markersize=2,linewidth=1,alpha=0.7)
    
    if colorbar==1:
        fig.colorbar(im,ax=ax)
    






from typing import Sequence, Union, Mapping, List, Optional, Dict, Callable
def quick_preprocess(
        #adata: sc.AnnData,
        adata,
        hvgs: Optional[Sequence] = None,
        normalize_data: bool = True,
        target_sum: Optional[float] = 1e4,
        batch_key=None,
        n_top_genes: int = 30000,
        n_pcs: int = 30,  # if None, stop before PCA step
        nneigh: int = 10,  # 20 was used for clustering
        metric='cosine',
        copy = True):
    """
    Quick preprocess of the raw data.

    Notes
    -----
    if `normalize_data` is True, an adata with RAW counts is required!
    """
    # if copy:
    #     _adata = adata.copy()
    #     logging.info('A copy of AnnData made!')
    # else:
    #     _adata = adata
    # 0. filtering not qualified genes or cells
    # gene with 1% cell expression and cell with 5% gene expression
    sc.pp.filter_genes(adata,min_cells=np.floor(adata.shape[0]*0.005)) 
    #sc.pp.filter_cells(adata,min_genes=np.floor(adata.shape[1]*0.01))
    # 1: normalization
    # if normalize_data:
    #     normalize_default(_adata, target_sum=target_sum)
    # _adata.raw = _adata
    # 2: HVG selection (skipped if `hvgs` is not None)
    if hvgs is None:
        sc.pp.highly_variable_genes(
            adata, batch_key=batch_key, n_top_genes=n_top_genes)
        adata = adata[:, adata.var['highly_variable']].copy()
    else:
        adata = adata[:, hvgs].copy()  # detach form view-data
    # 3: z-score 
    #wrapper_scale(_adata, groupby=batch_key)
    #    sc.pp.scale(_adata)
    # 4: PCA
    if n_pcs is None:
        do_pca = False
    else:
        sc.tl.pca(adata, n_comps=n_pcs)
        do_pca = True
    # 5: k-nearest-neighbor graph
    if do_pca and nneigh is not None:
        sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=nneigh, metric=metric)
    # 6: leiden-clustering...(separated function)
    sc.tl.leiden(adata)
    sc.tl.tsne(adata)
    sc.tl.umap(adata)
    return adata


#Cell_Feat_pd=all_insusrtength.fillna(value=1).T
def plot_schic2017_cluster(Cell_Feat_pd):
    matrix=Cell_Feat_pd.values
    adata=sc.AnnData(matrix,dtype=float)
    sc.pp.filter_genes(adata,min_cells=np.floor(adata.shape[0]*0.01))
    cellnames=Cell_Feat_pd.index.tolist()
    celltypes=[name.split('_')[0] for name in cellnames]
    
    sc.pp.normalize_total(adata, target_sum=adata.shape[0])
    
    n_pcs=50
    sc.tl.pca(adata, n_comps=n_pcs)
    #do_pca = True
    nneigh=20;metric='euclidean';
    sc.pp.neighbors(adata, n_pcs=n_pcs, n_neighbors=nneigh, metric=metric)
    sc.tl.leiden(adata)
    sc.tl.tsne(adata)
    sc.tl.umap(adata)
    sc.pl.tsne(adata, color='leiden')
    sc.pl.umap(adata, color='leiden')
    #cell target 1
    groups=celltypes
    adata.obs['cell-type cluster']=pd.Categorical(
        values=groups,
        categories=natsorted(map(str, np.unique(groups))),
    )
    sc.pl.tsne(adata, color='cell-type cluster')
    sc.pl.umap(adata, color='cell-type cluster')
    
    
       
#matrix=all_chrom_mat
def cluster_cells(matrix,cellnames,cids,indir):
    adata=sc.AnnData(matrix,dtype=float)
    adata=quick_preprocess(
            adata,
            hvgs = None,
            batch_key=None,
            n_top_genes= 30000,
            n_pcs= 20,  # if None, stop before PCA step
            nneigh= 20,  # 20 was used for clustering
            metric='cosine',
            copy = True)
    
    
    #select high var features
    sc.pl.tsne(adata, color='leiden')
    sc.pl.umap(adata, color='leiden')
    
    
    #showing original paper labels
    parent=os.path.join(os.path.join(indir, os.pardir),os.pardir)
    DipC_CellTypes=pd.read_csv(os.path.abspath(parent)+'/DipC_celltypes.csv',sep='\t',
                                header=0,index_col=None)
    #restore qualified cells
    Cell_index=[DipC_CellTypes['cell'].tolist().index(cell) for cell in cellnames]
    DipC_CellTypes_F=DipC_CellTypes.loc[Cell_index,:]
    
    #cell target 1
    Cell_label_1=DipC_CellTypes_F['cell-type cluster']
    groups=Cell_label_1.values
    adata.obs['cell-type cluster']=pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )
    
    #cell target 2
    Cell_label_2=DipC_CellTypes_F['age']
    groups=Cell_label_2.values
    adata.obs['age']=pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )
    
    #cell target 3
    Cell_label_3=DipC_CellTypes_F['tissue']
    groups=Cell_label_3.values
    adata.obs['tissue']=pd.Categorical(
        values=groups.astype('U'),
        categories=natsorted(map(str, np.unique(groups))),
    )
    
    
    sc.pl.tsne(adata, color='cell-type cluster')
    sc.pl.umap(adata, color='cell-type cluster')
    
    sc.pl.tsne(adata, color=['age','tissue'])
    sc.pl.umap(adata, color=['age','tissue'])
          
# #x=cell_corr; target=cell_target
# def plot_umap(x,target):
#     reducer = umap.UMAP(random_state=42)
#     embedding = reducer.fit_transform(x)
#     print(embedding.shape)
#     plt.scatter(embedding[:, 0], embedding[:, 1],c=target,cmap='Spectral', s=5)



    



    