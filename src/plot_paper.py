#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 10:56:33 2023

@author: dell
"""


import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import seaborn as sns
import scanpy as sc


#looplists=candidated_allloops
def remove_filter_regions(looplists,binsize,genome='mm10'):
    
    if genome=='mm10':
        remove_file='/mnt/nas/Datasets/genome_size/mm10_filter_regions.txt'
    
    fliter_regions=pd.read_csv(remove_file,sep='\t',header=None,index_col=None)
    
    fliter_regions['i']=(fliter_regions[1]//binsize).astype(int)
    fliter_regions['j']=(fliter_regions[2]//binsize).astype(int)
    
    fliter_regions_1=['_'.join([fliter_regions.loc[k,0],(fliter_regions.loc[k,1]//binsize).astype(int).astype(str),\
                                             ])  for k in fliter_regions.index]
    fliter_regions_2=['_'.join([fliter_regions.loc[k,0],(fliter_regions.loc[k,2]//binsize).astype(int).astype(str),\
                                             ])  for k in fliter_regions.index]
        
    fliter_regions=set(fliter_regions_1).union(fliter_regions_2)
    
    looplists.index=range(len(looplists))
    region_lists=['_'.join([looplists.loc[k,'chrom'],str(looplists.loc[k,'i'])])  for k in looplists.index]
    
    is_ture_1=[region not in fliter_regions for region in region_lists]
    region_lists=['_'.join([looplists.loc[k,'chrom'],str(looplists.loc[k,'j'])])  for k in looplists.index]
    is_ture_2=[region not in fliter_regions for region in region_lists]
    return looplists[pd.Series(is_ture_1) & pd.Series(is_ture_2)]

def return_fliter_regions(binsize,chrom,genome='mm10'):
    if genome=='mm10':
        remove_file='/mnt/nas/Datasets/genome_size/mm10_filter_regions.txt'
    fliter_regions=pd.read_csv(remove_file,sep='\t',header=None,index_col=None)
    fliter_regions=fliter_regions.loc[fliter_regions[0]==chrom]
    fliter_regions=np.sort((fliter_regions[1]//binsize).astype(int).unique())
    
    return fliter_regions

def detect_dis(vector,rang=20):
    for k in range(rang):
        x=vector[rang-k];y=vector[rang+k]
        if  (x!=0):
            return -k
        elif (y!=0):
            return k
    return rang

def plot_3D_density(Z,loc_i,loc_j,rang,color,ax):
    x=np.arange(loc_i-rang,loc_i+rang+1,1)
    y=np.arange(loc_j-rang,loc_j+rang+1,1)
    X,Y=np.meshgrid(x,y)
    #from mpl_toolkits import mplot3d
    from plot.plot_units import generate_cmap,white1
    mycmap=generate_cmap([white1,color])
    ax.plot_surface(X,Y,Z,cmap=mycmap,edgecolor='none')

def plot_2D_density(vec,rang,color,ax):
    plt.plot(range(len(vec)),vec,c=color,linestyle='-',linewidth=2)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show() 
    
    
    
    
def computing_enrich_mat(looplists,Peak_add,binsize,chrom_lens,rang=20,Figure=0):
    AllBin_Peaks=dict()
    Peaks=pd.read_csv(Peak_add,header=None,index_col=None,sep='\t')
    for chrom in chrom_lens.keys():
        chrom_len=chrom_lens[chrom]
        bin_range=int(chrom_len//binsize)+1
        #read chiip-seq peaks
        Bin_Peak_count=np.zeros(bin_range)
        BinPeak=pd.DataFrame()
        Peaks_chr=Peaks.loc[Peaks[0]==chrom]
        BinPeak[1]=(Peaks_chr[1]//binsize).astype(int)
        BinPeak[2]=(Peaks_chr[2]//binsize).astype(int)
        BinPeak[3]=Peaks_chr[9]
        BinPeak.index=range(len(BinPeak))
        for ind in range(len(BinPeak)):
            temp1=BinPeak.loc[ind,1];temp2=BinPeak.loc[ind,2];val=BinPeak.loc[ind,3]
            if temp1!=temp2:
                Bin_Peak_count[temp2]+=val/2
                Bin_Peak_count[temp1]+=val/2
            else:
                Bin_Peak_count[temp1]+=val
        AllBin_Peaks[chrom]=Bin_Peak_count

    Mat_Enrich=np.zeros((rang*2+1,rang*2+1))
    for chrom in chrom_lens.keys():
        chr_looplists=looplists.loc[looplists['chrom']==chrom]
        Bin_Peak_count=AllBin_Peaks[chrom]
        bin_range=int(chrom_len//binsize)+1
        for ind in chr_looplists.index:
            loop=chr_looplists.loc[ind,:]
            i=loop['i'];j=loop['j']
            if ((i-rang)>=0)&((j+rang)<bin_range):
                vector1=Bin_Peak_count[i-rang:i+rang+1]
                dis1=detect_dis(vector1)
                vector2=Bin_Peak_count[j-rang:j+rang+1]
                dis2=detect_dis(vector2)
                print(dis1,dis2,vector1[dis1+rang],vector1[dis1+rang])
                Mat_Enrich[dis1+rang,dis2+rang]+=(vector1[dis1+rang]+vector2[dis2+rang])
    Mat_Enrich/=len(looplists)
    if Figure==1:
        from plot.plot_units import color_3
        Z=Mat_Enrich
        x=np.arange(-rang,rang+1,1)
        y=np.arange(-rang,rang+1,1)
        X,Y=np.meshgrid(x,y)
        from mpl_toolkits import mplot3d
        fig=plt.figure()
        ax=plt.axes(projection='3d')
        from plot.plot_units import generate_cmap,white1,red
        mycmap=generate_cmap([white1,red])
        ax.plot_surface(X,Y,Z,cmap=mycmap,edgecolor='none')
    return Mat_Enrich
    
                
                




#looplists=weak_scloops;Peak_add=sc_filenames[2];
def computing_enrich_vector(looplists,Peak_add,binsize,chrom_lens,rang=20,Figure=0):
    AllBin_Peaks=dict()
    Peaks=pd.read_csv(Peak_add,header=None,index_col=None,sep='\t')
    for chrom in chrom_lens.keys():
        chrom_len=chrom_lens[chrom]
        bin_range=int(chrom_len//binsize)+1
        #read chiip-seq peaks
        Bin_Peak_count=np.zeros(bin_range)
        BinPeak=pd.DataFrame()
        Peaks_chr=Peaks.loc[Peaks[0]==chrom]
        BinPeak[1]=(Peaks_chr[1]//binsize).astype(int)
        BinPeak[2]=(Peaks_chr[2]//binsize).astype(int)
        BinPeak[3]=Peaks_chr[9]
        BinPeak.index=range(len(BinPeak))
        for ind in range(len(BinPeak)):
            temp1=BinPeak.loc[ind,1];temp2=BinPeak.loc[ind,2];val=BinPeak.loc[ind,3]
            if temp1!=temp2:
                Bin_Peak_count[temp2]+=val/2
                Bin_Peak_count[temp1]+=val/2
            else:
                Bin_Peak_count[temp1]+=val
        AllBin_Peaks[chrom]=Bin_Peak_count
    Enrich_Vector=np.zeros(rang*2+1)
    #Enrich_Vector2=np.zeros(rang*2+1)
    for chrom in chrom_lens.keys():
        chr_looplists=looplists.loc[looplists['chrom']==chrom]
        Bin_Peak_count=AllBin_Peaks[chrom]
        bin_range=int(chrom_len//binsize)+1
        for ind in chr_looplists.index:
            loop=chr_looplists.loc[ind,:]
            i=loop['i'];j=loop['j']
            # peak1=Bin_Peak_count[i]
            # peak2=Bin_Peak_count[j]
            # ri=i-np.random.randint(0,20)
            # rj=j+np.random.randint(0,20)
            if ((i-rang)>=0) & ((j+rang)<bin_range):
                # vector1=Bin_Peak_count[i-rang:i+rang+1]
                # vector2=Bin_Peak_count[j-rang:j+rang+1]
                vector=np.zeros(2*rang+1)
                vector[0:rang]=Bin_Peak_count[i-rang:i]
                vector[rang+1:2*rang+1]=Bin_Peak_count[j+1:j+rang+1]
                vector[rang]=Bin_Peak_count[i]
                Enrich_Vector+=vector 
                #Enrich_Vector2+=vector2
    Enrich_Vector/=len(looplists)
    return Enrich_Vector

    
                   
def plot_histplot(Data,num_thr,color_3):
    fig,ax=plt.subplots(figsize=(4,4))
    #plt.hist(inter_counts,bins=45)
    sns.histplot(Data,x=0,binwidth=1,color=color_3[1],alpha=0.75)
    plt.plot([num_thr,num_thr],[0,2000],c=color_3[0],linestyle='--',linewidth=2)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.set_ylim(0)
    plt.show()   


#data=bulk_map;data=sub_mat;looplist=snaphic_loops
def plot_apa_matrix(looplist,data,Figure=1,ra=5):
    # observe_list=list()
    # expect_list=list()
    for k,ind in enumerate(looplist.index):
        i=looplist.loc[ind,'i']
        j=looplist.loc[ind,'j']
        temp=data[i-ra:i+ra+1,j-ra:j+ra+1]
        if k==0:
            apa_matrix=temp
        else:
            apa_matrix+=temp
        # all_loops=apa_matrix.sum()
        # center_loops=apa_matrix[ra-2:ra+2+1,ra-2:ra+2+1].sum()
        # average=(all_loops-center_loops)/((ra*2+1)*(ra*2+1)-(5*5))
        # observe_list.append(temp[ra,ra])
        # expect_list.append(average)
    apa_matrix_log=np.log10(apa_matrix/len(looplist)+1)
    if Figure==1:
        from plot.plot_units import generate_cmap,white1,red
        mycmap=generate_cmap([white1,red])
        fig,ax=plt.subplots(1,1,figsize=(8,8))
        sns.heatmap(apa_matrix_log,cmap=mycmap,ax=ax,square=True)
    
    all_loops=apa_matrix_log.sum()
    center_loops=apa_matrix_log[ra-2:ra+2+1,ra-2:ra+2+1].sum()
    average=(all_loops-center_loops)/((ra*2+1)*(ra*2+1)-(5*5))
    strength=apa_matrix_log[ra,ra]/average
    # from scipy import stats
    # tind=stats.ttest_ind(observe_list,expect_list,equal_var=False)
    print(strength,np.max(apa_matrix_log))
    return strength
 
#looplists=blue_loops_hiccups  ;datas=bulk_maps
def plot_apa_matrix_allchr(looplists,datas,proc_chroms,Figure=1,ra=5):
    for c,chrom in enumerate(proc_chroms):
        looplist=looplists.loc[looplists['chrom']==chrom]
        data=datas[chrom]
        for k,ind in enumerate(looplist.index):
            i=looplist.loc[ind,'i']
            j=looplist.loc[ind,'j']
            temp=data[i-ra:i+ra+1,j-ra:j+ra+1]
            if (k==0) & (c==0):
                apa_matrix=temp/10000
            else:
                apa_matrix+=temp/10000
        # all_loops=apa_matrix.sum()
        # center_loops=apa_matrix[ra-2:ra+2+1,ra-2:ra+2+1].sum()
        # average=(all_loops-center_loops)/((ra*2+1)*(ra*2+1)-(5*5))
        # observe_list.append(temp[ra,ra])
        # expect_list.append(average)
    apa_matrix_log=np.log10(apa_matrix/len(looplist)*10000+1)
    if Figure==1:
        from plot.plot_units import generate_cmap,white1,red
        mycmap=generate_cmap([white1,red])
        fig,ax=plt.subplots(1,1,figsize=(8,8))
        sns.heatmap(apa_matrix_log,cmap=mycmap,ax=ax,square=True)
    
    all_loops=apa_matrix_log.sum()
    center_loops=apa_matrix_log[ra-2:ra+2+1,ra-2:ra+2+1].sum()
    average=(all_loops-center_loops)/((ra*2+1)*(ra*2+1)-(5*5))
    strength=apa_matrix_log[ra,ra]/average
    # from scipy import stats
    # tind=stats.ttest_ind(observe_list,expect_list,equal_var=False)
    #print(strength,np.max(apa_matrix_log))
    return strength

           
#data=loop_num
def plot_displot(data,x='alldistance'):
    
    #fig = plt.figure(figsize =(4, 4))
    # ax = fig.add_subplot(111) 
    #fig,ax=plt.subplots(figsize=(4,4))
    from plot.plot_units import color_3
    plt.figure(figsize=(4,4))
    ax=sns.displot(data, x=x, fill=True,color=color_3[0],kde=True)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()

#data=loopnum_sc
def plot_displot_2(data,bins=100):
    
    #plot single cell loop number distribution
    from plot.plot_units import color_3
    #sns.set(rc={"figure.figsize":(4,4)})
    sns.displot(data,kde=False,bins=bins,color=color_3[0])
    #plt.figure(figsize=(4,4))
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()
    
    

    

#data=region_ratios
# def plot_boxplot_2(data,proc_chroms):
#     fig = plt.figure(figsize =(4, 4))
#     ax = fig.add_subplot(111) 
#     # Creating axes instance
#     bp = ax.boxplot(data, patch_artist = True,
#                       vert = 1,showfliers=False,widths=0.5, 
#                       medianprops={"color": "white", "linewidth": 0.5},
#                 boxprops={"facecolor": "grey", "edgecolor": "white",
#                            "linewidth": 0.5,"alpha":0.8},
#                 whiskerprops={"color": "grey", "linewidth": 1.0},
#                 capprops={"color": "grey", "linewidth": 1.0})
    
#     from plot.plot_units import color_3
#     colors = color_3 
#     for patch, color in zip(bp['boxes'], colors):
#         patch.set_facecolor(color)
    
#     import matplotlib.cm as cm
#     colors=cm.tab20.colors
#     Y1=data
#     X1=np.random.uniform(-0.35,0.35,len(Y1))+1
    
#     for ch in range(len(X1)):
#         ax.scatter(X1[ch],Y1[ch],c=colors[ch],alpha=0.9,s=10)        
    
#     ax = plt.gca()
#     ax.spines['right'].set_color('none')
#     ax.spines['top'].set_color('none')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.yaxis.set_ticks_position('left')
#     plt.show()




#data=cell_feadis_ratio
def plot_boxplot(data):
    fig = plt.figure(figsize =(4, 4))
    ax = fig.add_subplot(111) 
    # Creating axes instance
    bp = ax.boxplot(data, patch_artist = True,
                      vert = 1,showfliers=False,widths=0.5,
                      medianprops={"color": "white", "linewidth": 0.5},
                boxprops={"facecolor": "grey", "edgecolor": "white",
                           "linewidth": 0.5,"alpha":0.8},
                whiskerprops={"color": "grey", "linewidth": 1.0},
                capprops={"color": "grey", "linewidth": 1.0})
    
    if data.shape[1]<=3:
        from plot.plot_units import color_3
        colors = color_3 
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
    
    import matplotlib.cm as cm
    colors=cm.tab20.colors
    #add points
    for k,col in enumerate(data.columns):
        Y1=data[col].values
        X1=np.random.uniform(-0.35,0.35,size=len(Y1))+1+k
        for ch in range(len(X1)):
            ax.scatter(X1[ch],Y1[ch], c=colors[ch],alpha=0.9,s=5)        
    
    
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()
    
#data=looptype_ratio
def plot_boxplot_2(data):
    fig = plt.figure(figsize =(4, 4))
    ax = fig.add_subplot(111) 
    # Creating axes instance
    bp = ax.boxplot(data, patch_artist = True,
                      vert = 1,showfliers=False,widths=0.5,
                      medianprops={"color": "white", "linewidth": 0.5},
                boxprops={"facecolor": "grey", "edgecolor": "white",
                           "linewidth": 0.5,"alpha":0.8},
                whiskerprops={"color": "grey", "linewidth": 1.0},
                capprops={"color": "grey", "linewidth": 1.0})
    
    if len(data.columns)<=3:
        from plot.plot_units import color_3
        colors = color_3 
    else:
        import matplotlib.cm as cm
        colors=cm.tab20.colors    
    # for patch, color in zip(bp['boxes'], colors):
    #     patch.set_facecolor(color)
    # import matplotlib.cm as cm
    # colors=cm.tab20.colors
    #color='gray'
    #add points
    for k,col in enumerate(data.columns):
        Y1=data[col].values
        X1=np.random.uniform(-0.35,0.35,size=len(Y1))+1+k
        for ch in range(len(X1)):
            ax.scatter(X1[ch],Y1[ch], c=colors[k],alpha=0.2,s=3)        
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()

def plot_boxplot_xy(X,Y,colors=['#CC3939']):
    Tem=pd.DataFrame(X)
    Tem['type']='X'
    Tem2=pd.DataFrame(Y)
    Tem2['type']='Y'
    Tem=pd.concat([Tem,Tem2])
    Tem.columns=['val','type']
    
    fig,axes=plt.subplots(1,1,figsize=(4,4))
    sns.boxplot(Tem,x='type',y='val',palette=colors,boxprops={"alpha":0.9},fliersize=0)
    # for patch, color in zip(bp['boxes'], colors_3):
    #     patch.set_facecolor(color,alpha=0.75)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.legend('')
    plt.show() 



#data=Result_PD
def plot_boxplot_3(data,x='type',y='distance',colors=['#CC3939'],hue=None):
    
    #from plot.plot_units import color_3
    fig,axes=plt.subplots(1,1,figsize=(4,4))
    sns.boxplot(data,x=x,y=y,hue=hue,palette=colors,boxprops={"alpha":0.9},fliersize=0)
    # for patch, color in zip(bp['boxes'], colors_3):
    #     patch.set_facecolor(color,alpha=0.75)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.legend('')
    plt.show() 

def plot_umap_scanpy(adata,styles=['cell-type cluster','age','tissue','leiden'],color_list='tab20',size=20):
    
    #styles=['cell-type cluster','age','tissue','leiden']
    x=len(styles)*3+1;y=3
    fig,axes=plt.subplots(1,len(styles),figsize=(x,y))
    for k,style in enumerate(styles):
        ax=axes[k]
        sc.pl.umap(adata,color=style,palette=color_list,size=size,ax=ax,wspace=0)
        #ax.legend('',frameon=False) 
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        plt.show()




#adata=adata_100.copy()
def plot_umap_adata(adata,stage,color_list_temp,ax,size=15):
    X1=adata.obsm['X_umap'][:,0]
    Y1=adata.obsm['X_umap'][:,1]
    for c,celltype in zip(color_list_temp,np.sort(adata.obs[stage].unique())):
        X1_temp=X1[adata.obs[stage]==celltype]
        Y1_temp=Y1[adata.obs[stage]==celltype]
        ax.scatter(X1_temp,Y1_temp, color=c,alpha=1,s=size)


def plot_umap_adata_v2(adata,celltype,color_list_stages,ax,size=1):
    X1=adata.obsm['X_umap'][:,0]
    Y1=adata.obsm['X_umap'][:,1]
    from plot.plot_units import white3
    X1_temp=X1[(adata.obs['Celltype']!=celltype)]
    Y1_temp=Y1[(adata.obs['Celltype']!=celltype)]
    ax.scatter(X1_temp,Y1_temp, color=white3,alpha=1,s=size)
    
    for c,stage in zip(color_list_stages,adata.obs['Stage'].unique().sort_values()):
        X1_temp=X1[(adata.obs['Celltype']==celltype) & (adata.obs['Stage']==stage)]
        Y1_temp=Y1[(adata.obs['Celltype']==celltype) & (adata.obs['Stage']==stage)]
        ax.scatter(X1_temp,Y1_temp, color=c,alpha=1,s=size)
    
    
    
    


#data=cell_feadis_ratio
def plot_boxplot_4(data):
    
    fig = plt.figure(figsize =(4, 4))
    ax = fig.add_subplot(111) 
    from plot.plot_units import color_3
    # Creating axes instance
    bp = ax.boxplot(data, patch_artist = True,
                      vert = 1,showfliers=False,widths=0.5,
                      medianprops={"color": "black", "linewidth": 1},
                boxprops={"facecolor": color_3[0],"alpha":1 ,"edgecolor": "black",
                           "linewidth": 1},
                whiskerprops={"color": "black", "linewidth": 1.0},
                capprops={"color": "black", "linewidth": 1.0})
    
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()
    
def plot_boxplot_5(data):
    fig = plt.figure(figsize =(4,4))
    ax = fig.add_subplot(111) 
    # Creating axes instance
    bp = ax.boxplot(data, patch_artist = True,
                      vert = 1,showfliers=False,widths=0.5,
                      medianprops={"color": "black", "linewidth": 1,"alpha":0.7},
                boxprops={"facecolor": 'gray',"alpha":0.5 ,"edgecolor": "black",
                           "linewidth": 1,"alpha":0.7 },
                whiskerprops={"color": "black", "linewidth": 1.0,"alpha":0.7},
                capprops={"color": "black", "linewidth": 1.0,"alpha":0.7})
    
    if len(data.columns)<=3:
        from plot.plot_units import color_3
        colors = color_3 
    else:
        import matplotlib.cm as cm
        colors=cm.tab20.colors    
    # for patch, color in zip(bp['boxes'], colors):
    #     patch.set_facecolor(color)
    # import matplotlib.cm as cm
    # colors=cm.tab20.colors
    #color='gray'
    #add points
    for k,col in enumerate(data.columns):
        Y1=data[col].values
        X1=np.random.uniform(-0.35,0.35,size=len(Y1))+1+k
        ax.scatter(X1,Y1, c=colors[k],alpha=0.3,s=0.2)        
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()


    
    
    



#data=loop_num
def plot_vionlinplot(data,figsize=(4,4)):
    fig = plt.figure(figsize =figsize)
    ax = fig.add_subplot(111) 
    # Creating axes instance
    bp = ax.violinplot(data,
                     vert=1,widths=0.628,showmeans=False, showmedians=False, showextrema=False)
 
    if data.shape[1]<=3:
        from plot.plot_units import color_3
        colors = color_3 
        for patch, color in zip(bp['bodies'], colors):
            patch.set_facecolor(color)
    
    
    import matplotlib.cm as cm
    colors=cm.tab20.colors
    #add points
    for col in data.columns:
        Y1=data[col].values
        X1=np.random.uniform(-0.35,0.35,size=len(Y1))+1+col
        for ch in range(len(X1)):
            ax.scatter(X1[ch],Y1[ch], c=colors[ch],alpha=0.9,s=5)        
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()
    #plot_color_gradients('Qualitative',['tab20'])



# data=loop_num_dis
def plot_vionlinplot_2(data,figsize=(4,4)):
    
    fig = plt.figure(figsize =figsize)
    ax = fig.add_subplot(111) 
    # Creating axes instance
    bp = ax.violinplot(data,
                     vert=1,widths=0.628,showmeans=False, showmedians=False, showextrema=False)
    import matplotlib.cm as cm
    colors=cm.tab20.colors    
    for patch, color in zip(bp['bodies'], colors):
        patch.set_facecolor(color)
    
    #add points
    for col in data.columns:
        Y1=data[col].values
        X1=np.random.uniform(-0.4,0.4,size=len(Y1))+col-1
        for ch in range(len(X1)):
            ax.scatter(X1[ch],Y1[ch], c='grey',alpha=0.8,s=3)        
    
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()
    #plot_color_gradients('Qualitative',['tab20'])
    

# def plot_mergeloop_chrom(merge_loops,chrom,figsize=(10,4)):    
#     fig = plt.figure(figsize =figsize)
#     ax = fig.add_subplot(111) 
#     using_loops=merge_loops.loc[merge_loops['chrom']==chrom]
#     #data1
#     temp_loops=using_loops.loc[ (using_loops['scloops_flag']==1) & (using_loops['snaphic_flag']==0)]
#     Y1=temp_loops['log10(Freq_inter+1)'].values
#     X1=np.random.uniform(-0.4,0.4,size=len(Y1))+temp_loops['norm_distance_floor'].values 
#     for ch in range(len(X1)):
#         ax.scatter(X1[ch],Y1[ch], c='grey',alpha=0.8,s=3)      
    
#data=sub_matrix.sum(axis=1)
def plot_lineplot(data):
    from plot.plot_units import color_3
    pos_strength_pd=pd.DataFrame()
    pos_strength_pd['x']=range(len(data))
    pos_strength_pd['y']=data
    fig,ax=plt.subplots(figsize=(8,2))
    sns.lineplot(data=pos_strength_pd,x='x',y='y',color=color_3[0],ax=ax)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    plt.show()

#data=metrics_PD2[data.columns[0:6]]
def plot_multilineplot(data):   
    from plot.plot_units import colors_nature
    fig,ax=plt.subplots(figsize=(4,4))
    for k,col in enumerate(data.columns[0:6]):
        ax.plot(data.index,data[col],color=colors_nature[k],linestyle='--',marker='.',label=col)
        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        plt.show()
    ax.legend()
   
    
    
    
cmaps = {}
gradient = np.linspace(0, 1, 256)
gradient = np.vstack((gradient, gradient))
def plot_color_gradients(category, cmap_list):
    # Create figure and adjust figure height to number of colormaps
    nrows = len(cmap_list)
    figh = 0.35 + 0.15 + (nrows + (nrows - 1) * 0.1) * 0.22
    fig, axs = plt.subplots(nrows=nrows + 1, figsize=(6.4, figh))
    fig.subplots_adjust(top=1 - 0.35 / figh, bottom=0.15 / figh,
                        left=0.2, right=0.99)
    axs[0].set_title(f'{category} colormaps', fontsize=14)

    for ax, name in zip(axs, cmap_list):
        ax.imshow(gradient, aspect='auto', cmap=mpl.colormaps[name])
        ax.text(-0.01, 0.5, name, va='center', ha='right', fontsize=10,
                transform=ax.transAxes)

    # Turn off *all* ticks & spines, not just the ones with colormaps.
    for ax in axs:
        ax.set_axis_off()

    # Save colormap list for later.
    cmaps[category] = cmap_list
    

def plot_boxplot_interfreq(merge_loops,method='snaphic'):
    using_loops=merge_loops
    
    flag=method+'_flag'
    
    all_pd=pd.DataFrame()
    for k,dis in enumerate(np.sort(using_loops['norm_distance_floor'].unique())):
        
        # temp_pd0=pd.DataFrame()
        # Y2=using_loops.loc[(using_loops['snaphic_flag']==False)&(using_loops['scloops_flag']==True)& (using_loops['norm_distance_floor']==dis),
        #                    'log10(Freq_inter+1)'].values
        # temp_pd0['log10(Freq_inter+1)']=Y2
        # temp_pd0['method']='other'
        # temp_pd0['norm_distance_floor']=dis
        temp_pd2=pd.DataFrame()
        Y1=using_loops.loc[(using_loops[flag]==True)&(using_loops['scloops_flag']==True)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        temp_pd2['log10(Freq_inter+1)']=Y1
        temp_pd2['method']='overlap'
        temp_pd2['norm_distance_floor']=dis
        
        temp_pd1=pd.DataFrame()
        Y3=using_loops.loc[(using_loops[flag]==True)&(using_loops['scloops_flag']==False)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        temp_pd1['log10(Freq_inter+1)']=Y3
        temp_pd1['method']='nooverlap'
        temp_pd1['norm_distance_floor']=dis
        
        temp_pd=pd.concat([temp_pd1,temp_pd2])
        all_pd=pd.concat([all_pd,temp_pd])
    
    from plot.plot_units import generate_cmap,color_3
    fig,axes=plt.subplots(1,1,figsize=(4,4))
    bp=sns.boxplot(all_pd,x='norm_distance_floor',y='log10(Freq_inter+1)',hue='method',palette=color_3,\
                hue_order=['overlap','nooverlap'],fliersize=0,
                boxprops={"alpha":0.9})
    # for patch, color in zip(bp['boxes'], colors_3):
    #     patch.set_facecolor(color,alpha=0.75)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.legend('')
    plt.show()    

       
def plot_volcano(PD,p_thr=0.05,fd_change=0.05):
    #PD=Strong_Specelltype_loop_pd.copy()
    from plot.plot_units import color_3
    #colors = color_3 
    fig,ax=plt.subplots(figsize=(4,4))
    
    loop_ind1= (PD['pval']<=p_thr) & (PD['fd']>1+fd_change)
    
    X1=-np.log10(PD.loc[loop_ind1,'pval'].astype(float))
    Y1=PD.loc[loop_ind1,'fd']
    ax.scatter(Y1,X1, c=color_3[0],alpha=0.8,s=5)
    
    loop_ind2= (PD['pval']<=p_thr) & (PD['fd']<(-1-fd_change))
    X2=-np.log10(PD.loc[loop_ind2,'pval'].astype(float))
    Y2=PD.loc[loop_ind2,'fd']
    ax.scatter(Y2,X2, c=color_3[1],alpha=0.8,s=5)
    
    loop_ind3=  (~loop_ind1) & ( ~loop_ind2)
    
    X3=-np.log10(PD.loc[loop_ind3,'pval'].astype(float))
    Y3=PD.loc[loop_ind3,'fd']
    ax.scatter(Y3,X3, c='grey',alpha=0.2,s=1)

#PD=CT_temp_PD
def plot_volcano_2(PD,genelist=[],p_thr=0.05,fd_change=0.5):
    #PD=Strong_Specelltype_loop_pd.copy()
    
    
    from plot.plot_units import color_3
    #colors = color_3 
    fig,ax=plt.subplots(figsize=(4,4))
    loop_ind1= (PD['pval']<=p_thr) & (PD['fd']>1+fd_change)
    X1=-np.log10(PD.loc[loop_ind1,'pval'].astype(float))
    Y1=PD.loc[loop_ind1,'lgfd']
    ax.scatter(Y1,X1, c=color_3[0],alpha=0.8,s=10)
    
    
    loop_ind2= (PD['pval']<=p_thr) & (PD['fd']<(-1-fd_change))
    X2=-np.log10(PD.loc[loop_ind2,'pval'].astype(float))
    Y2=PD.loc[loop_ind2,'lgfd']
    ax.scatter(Y2,X2, c=color_3[1],alpha=0.8,s=10)
    
    loop_ind3=  (~loop_ind1) & ( ~loop_ind2)
    X3=-np.log10(PD.loc[loop_ind3,'pval'].astype(float))
    Y3=PD.loc[loop_ind3,'lgfd']
    ax.scatter(Y3,X3, c='grey',alpha=0.2,s=5)
    
    loop_ind4= (PD['pval']<=p_thr) & ((PD['fd']<-2)| (PD['fd']>2))
    X4=-np.log10(PD.loc[loop_ind4,'pval'].astype(float))
    Y4=PD.loc[loop_ind4,'lgfd']
    temp_PD=PD.loc[loop_ind4]
    #loop=temp_PD.index[0]
    for k,loop in enumerate(temp_PD.index):
        #loopvals=temp_PD[loop]
        if temp_PD.loc[loop,'gene'] in genelist:
            text=loop+' ('+temp_PD.loc[loop,'gene']+')'
            print(text,Y4[k],X4[k])
            if temp_PD.loc[loop,'direct']>0:
                xy=(Y4[k],X4[k]); xytext=(Y4[k]-0.5,X4[k]-0.2) 
            else:
                xy=(Y4[k],X4[k]); xytext=(Y4[k]-0.5,X4[k]+0.2) 
            plt.annotate(text, xy, xytext,fontsize=6)
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.legend('')
    plt.show()    
    return loop_ind1,loop_ind2,loop_ind3
    
    
    
def plot_scatterplots(merge_loops,chrom):
    
    from plot.plot_units import color_3
    colors = color_3 
    
    using_loops=merge_loops.loc[merge_loops['chrom']==chrom]
    import matplotlib.cm as cm
    colors=cm.tab20.colors   
    
    #only consider snaphic
    fig,axes=plt.subplots(1,7,figsize=(8,3))
    for k,dis in enumerate(np.sort(using_loops['norm_distance_floor'].unique())):
        # i=k%7
        # j=k//7
        ax=axes[k]
        Y2=using_loops.loc[(using_loops['snaphic_flag']==False)&(using_loops['scloops_flag']==True)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        X2=np.random.uniform(0,20,size=len(Y2))
        ax.scatter(X2,Y2, c='gray',alpha=0.3,s=3)
        
        Y3=using_loops.loc[(using_loops['snaphic_flag']==True)&(using_loops['scloops_flag']==False)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        X3=np.random.uniform(0,20,size=len(Y3))
        ax.scatter(X3,Y3, c=color_3[1],alpha=0.8,s=5)
        Y1=using_loops.loc[(using_loops['snaphic_flag']==True)&(using_loops['scloops_flag']==True)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        X1=np.random.uniform(0,20,size=len(Y1))
        ax.scatter(X1,Y1, c=color_3[0],alpha=0.8,s=5)
        ax.set_title(str(dis))
        ax.set_ylim(0.5,3.5)
        ax.set_xticks([])
        if k!=0:
            ax.set_yticks([])
    plt.subplots_adjust(wspace=0.05)
       
    
    #sns.scatterplot(x='norm_distance',y='Freq_inter')
    fig,axes=plt.subplots(1,7,figsize=(8,3))
    for k,dis in enumerate(np.sort(using_loops['norm_distance_floor'].unique())):
        # i=k%5
        # j=k//5
        ax=axes[k]
        Y2=using_loops.loc[(using_loops['hiccups_flag']==False)&(using_loops['scloops_flag']==True)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        X2=np.random.uniform(0,20,size=len(Y2))
        ax.scatter(X2,Y2, c='gray',alpha=0.3,s=3)
        
        Y3=using_loops.loc[(using_loops['hiccups_flag']==True)&(using_loops['scloops_flag']==False)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        X3=np.random.uniform(0,20,size=len(Y3))
        ax.scatter(X3,Y3, c=color_3[1],alpha=0.8,s=5)
        Y1=using_loops.loc[(using_loops['hiccups_flag']==True)&(using_loops['scloops_flag']==True)& (using_loops['norm_distance_floor']==dis),
                           'log10(Freq_inter+1)'].values
        X1=np.random.uniform(0,20,size=len(Y1))
        ax.scatter(X1,Y1, c=color_3[0],alpha=0.8,s=5)
        ax.set_title(str(dis))
        ax.set_ylim(0.5,3.5)
        ax.set_xticks([])
        if k!=0:
            ax.set_yticks([])
    plt.subplots_adjust(wspace=0.05)
        

        

def show_venn(snaphic_allloops,hiccups_loops,scloops_allloops):
    
    
    chrom='chr2'
    snaphic_loops_chr=snaphic_allloops[snaphic_allloops['chrom']==chrom]
    hiccups_loops_chr=hiccups_loops[hiccups_loops['chrom']==chrom]
    scloops_loops_chr=scloops_allloops[scloops_allloops['chrom']==chrom]
    
    
    start=500;end=7000
    snaphic_temp= snaphic_loops_chr.loc[( snaphic_loops_chr['i']>=start)&( snaphic_loops_chr['j']>=start)&( snaphic_loops_chr['i']<end)\
                        &( snaphic_loops_chr['j']<end)]   
    hiccups_temp=hiccups_loops_chr.loc[(hiccups_loops_chr['i']>=start)&(hiccups_loops_chr['j']>=start)&(hiccups_loops_chr['i']<end)\
                        &(hiccups_loops_chr['j']<end)]
    mine_loops_temp=scloops_loops_chr.loc[(scloops_loops_chr['i']>=start)&(scloops_loops_chr['j']>=start)&(scloops_loops_chr['i']<end)\
                        &(scloops_loops_chr['j']<end)]   
    
        
    fig,axes=plt.subplots(1,5)
    for k,dis in enumerate([4,5,6,7,8]):
        print(dis)
        from matplotlib_venn import venn3
        set1=set(mine_loops_temp.loc[mine_loops_temp['norm_distance_floor']==dis,'fea_names'])
        set2=set(snaphic_temp.loc[snaphic_temp['norm_distance_floor']==dis,'fea_names'])
        set3=set(hiccups_temp.loc[hiccups_temp['norm_distance_floor']==dis,'fea_names'])
        subsets=[set1,set2,set3]
        g=venn3(subsets,('scloops','snaphic','hiccups'),ax=axes[k])
        
        
    set2=set(snaphic_temp['fea_names'])
    set3=set(hiccups_temp['fea_names'])
    x=[]
    y=[]
    for thr in np.arange(0,20,2):
        temp=mine_loops_temp[mine_loops_temp['weight']>thr]
        set1=set(temp['fea_names'])
        x.append(len(set1.intersection(set2))/len(set2))
        y.append(len(set1.intersection(set3))/len(set3))
    
#counts=temp    
def plot_barplot(counts):
   
    from plot.plot_units import color_3
    colors=color_3
    fig,ax=plt.subplots(figsize=(4,4))
    ax.bar(range(len(counts)),counts,color=colors[0])
    ax = plt.gca()
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    #ax.legend('')
    plt.show()    
       
    
    
    
    
    
    
    
    
    
        
    
    
    
    
    