# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 17:04:19 2021

@author: user
"""

import os
import numpy as np
import scipy as sp
#from scipy import stats
import pandas as pd
import skimage
#import subprocess
#from statsmodels.stats.multitest import multipletests
import gc
#import sys
#import h5py
import seaborn as sns
from scipy import sparse
import re
import matplotlib.pyplot as plt

#import time
def get_proc_chroms(chrom_lens, rank, n_proc):
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    chrom_list.sort(key=lambda x: x[1])
    chrom_list.reverse()
    chrom_names = [i[0] for i in chrom_list]
    #chrom_names = list(chrom_lens.keys())
    #chrom_names.sort()
    
    indices = list(range(rank, len(chrom_names), n_proc))
    proc_chroms = [chrom_names[i] for i in indices]
    return proc_chroms


def get_nth_diag_indices(mat, offset):
    rows, cols_orig = np.diag_indices_from(mat)
    cols = cols_orig.copy()
    if offset > 0:
        cols += offset
        rows = rows[:-offset]
        cols = cols[:-offset]
    return rows, cols

from scipy.special import expit
def normalize_along_diagonal_from_numpy_v2(d,max_bin_distance, trim = 0.01):
    #d=strengthes.copy()
    for offset in range(1, d.shape[0]):
        r, c = get_nth_diag_indices(d, offset)
        if offset<=max_bin_distance:
            vals = d[r,c].tolist()
            #vals_orig=[val for val in vals_orig if val > 0]
            #vals = vals_orig.copy()
            #vals.sort()
            #vals.reverse()
            #trim_value = vals[(round(trim * len(vals)) - 1)]
            #trim_index = round(trim * len(vals)) - 1
            #remaining = vals[(trim_index):]
            #remaining=vals
            count=np.sum(np.array(vals)!= 0)
            #print(offset, count)
            if count!=0:
                mu= np.sum(vals)/count
                #mu = np.mean(remaining)
                temp=np.array(vals)
                sd = np.std(temp[temp>0])
                values=(np.array(vals)-mu)/sd
                #values=expit(values)
                d[r,c]=values
                print(offset, mu)
            else:
                d[r,c]=0
        else:
            d[r,c]=0
    return d           

#d=strengthes,local_den,global_den
def normlization_v2(d, dist, binsize):
    max_bin_distance = int((dist) // binsize)
    #df = df[(df['x1'] >= 50000) & (df['y1'] <= last_bin * binsize - 50000)]
    #output_filename = os.path.join(outdir, ".".join([setname, chrom, "normalized", "rwr", "bedpe"]))
    d=normalize_along_diagonal_from_numpy_v2(d,max_bin_distance)
    #d=d+d.T
    return d


def sqrt_norm(matrix):
    coverage = (np.sqrt(np.sum(matrix, axis=-1)))
    with np.errstate(divide='ignore', invalid='ignore'):
        matrix = matrix / coverage.reshape((-1, 1))
        matrix = matrix / coverage.reshape((1, -1))
    matrix[np.isnan(matrix)] = 0.0
    matrix[np.isinf(matrix)] = 0.0
    return matrix

def kth_diag_indices(a, k):
    rows, cols = np.diag_indices_from(a)
    if k < 0:
        return rows[-k:], cols[:k]
    elif k > 0:
        return rows[:-k], cols[k:]
    else:
        return rows, cols

def oe(matrix, expected = None):
    new_matrix = np.zeros_like(matrix)
    for k in range(len(matrix)):
        rows, cols = kth_diag_indices(matrix, k)
        diag = np.diag(matrix,k)
        if expected is not None:
            expect = expected[k]
        else:
            expect = np.sum(diag) / (np.sum(diag != 0.0) + 1e-15)
        if expect == 0:
            new_matrix[rows, cols] = 0.0
        else:
            new_matrix[rows, cols] = diag / (expect)
    new_matrix = new_matrix + new_matrix.T
    return new_matrix

#import numpy as np
#mat=np.random.rand(50,50,10)
def computing_interactions_strengthes(sub_mat,upper_limit=5,lower_limit=2):
    # mat_size=sub_mat.shape[0]
    # RANGE=range(upper_limit,mat_size-upper_limit)
    # sub_mat_pd=sp.sparse.coo_matrix(sub_mat[RANGE,:][:,RANGE])
    # sub_mat_pd=pd.DataFrame({'i': sub_mat_pd.row+upper_limit, 'j': sub_mat_pd.col+upper_limit,'value':sub_mat_pd})
    # sub_mat_pd['strength']=0
    # for index,key in sub_mat_pd.iterrows():
    #     i=key['i'];j=key['j'];#value=key['value']
    #     print(index,i,j)
    #     big=sub_mat[i-upper_limit:i+upper_limit+1,j-upper_limit:j+upper_limit+1]
    #     small=sub_mat[i-lower_limit:i+lower_limit+1,j-lower_limit:j+lower_limit+1]
    #     big_count = np.sum(~np.isnan(big))
    #     small_count = np.sum(~np.isnan(small))
    #     local=(big.sum()-small.sum())/(big_count-small_count)
    #     sub_mat_pd.loc[index,'strength']=sub_mat[i,j]/local 
    sub_mat= np.triu(sub_mat, 1)
    a = np.empty(sub_mat.shape)
    a[:] = np.nan
    sub_mat= np.tril(a,0)+sub_mat
    #sub_mat=np.float16(sub_mat)
    
    big_neighborhood = skimage.util.view_as_windows(sub_mat, (2*upper_limit+1,2*upper_limit+1), step = 1)
    small_neighborhood = skimage.util.view_as_windows(sub_mat, (2*lower_limit+1,2*lower_limit+1), step = 1)
    min_neighborhood = skimage.util.view_as_windows(sub_mat, (2*1+1,2*1+1), step = 1)
    #reshape to (matsize, matsize, numcells, num_neighbors+1)
    big_neighborhood = np.squeeze(big_neighborhood.reshape(big_neighborhood.shape[0],\
                                               big_neighborhood.shape[0], 1, -1))
    
    small_neighborhood = np.squeeze(small_neighborhood.reshape(small_neighborhood.shape[0],\
                                               small_neighborhood.shape[0], 1, -1))  
    min_neighborhood = np.squeeze(min_neighborhood.reshape(min_neighborhood.shape[0],\
                                               min_neighborhood.shape[0], 1, -1))    
    
    big_neighborhood_counts = np.sum(~np.isnan(big_neighborhood), axis = -1)
    small_neighborhood_counts = np.sum(~np.isnan(small_neighborhood), axis = -1)
    min_neighborhood_counts = np.sum(~np.isnan(min_neighborhood), axis = -1)
     
    min_nonzero_neighborhood_counts = np.sum((~np.isnan(min_neighborhood)) & (min_neighborhood!=0), axis = -1)
    
    
    
    big_neighborhood = np.nansum(big_neighborhood, axis = -1)
    small_neighborhood = np.nansum(small_neighborhood,axis = -1)
    min_neighborhood = np.nansum(min_neighborhood,axis = -1)
    #remove edge cases that are used only as neighbors
    trim_size = upper_limit - lower_limit
    small_neighborhood = small_neighborhood[trim_size:-trim_size, trim_size:-trim_size]
    small_neighborhood_counts = small_neighborhood_counts[trim_size:-trim_size, trim_size:-trim_size]
    
    trim_size = upper_limit - 1
    min_neighborhood = min_neighborhood[trim_size:-trim_size, trim_size:-trim_size]
    min_neighborhood_counts = min_neighborhood_counts[trim_size:-trim_size, trim_size:-trim_size]
    
    min_nonzero_neighborhood_counts=min_nonzero_neighborhood_counts[trim_size:-trim_size, trim_size:-trim_size]
    
    local_den=min_neighborhood/min_neighborhood_counts
    local_den[min_nonzero_neighborhood_counts<=1]=0
    
    local_neighborhood = big_neighborhood - small_neighborhood
    local_neighborhood_counts = big_neighborhood_counts - small_neighborhood_counts    
    del small_neighborhood, big_neighborhood, big_neighborhood_counts, small_neighborhood_counts,min_neighborhood,min_neighborhood_counts
    gc.collect()
    local_neighborhood /= local_neighborhood_counts
    del local_neighborhood_counts
    gc.collect()
    
    sub_mat = sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    #local_interaction_strengthes=sub_mat/local_neighborhood
    #strengthes=(sub_mat/((sub_mat>0)*(local_neighborhood>0.01)*local_neighborhood))
    
    local_neighborhood[local_neighborhood<np.nanmean(local_neighborhood)]=0
    
    strengthes=(sub_mat)/((sub_mat>0)*(local_neighborhood>0)*local_neighborhood)
    #strengthes=(small_neighborhood/(2*lower_limit+1)/((small_neighborhood>4*sub_mat)*(local_neighborhood>0.01)*local_neighborhood))
    #small_neighborhood/(2*lower_limit+1)/(2*lower_limit+1)
    strengthes[np.isnan(strengthes)] = 0
    strengthes[np.isinf(strengthes)] = 0
    strengthes[strengthes<1.5]=0
    local_den[strengthes==0] = 0
    #strengthes_pd=sp.sparse.coo_matrix(strengthes)
    #strengthes_pd=pd.DataFrame({'i': strengthes_pd.row+upper_limit, 'j': strengthes_pd.col+upper_limit,'den':strengthes_pd.data})
    #sns.displot(strengthes_pd["den"]) 
    return strengthes, local_den


def computing_sc_interactions_den(sc_sub_mat,upper_limit):
    sc_sub_mat= np.triu(sc_sub_mat, 1)
    a = np.empty(sc_sub_mat.shape)
    a[:] = np.nan
    sc_sub_mat= np.tril(a,0)+sc_sub_mat
    
    min_neighborhood = skimage.util.view_as_windows(sc_sub_mat, (2*1+1,2*1+1), step = 1)
    min_neighborhood = np.squeeze(min_neighborhood.reshape(min_neighborhood.shape[0],\
                                               min_neighborhood.shape[0], 1, -1)) 
    min_neighborhood_counts = np.sum(~np.isnan(min_neighborhood), axis = -1)
    min_neighborhood = np.nansum(min_neighborhood,axis = -1)
    trim_size = upper_limit - 1
    min_neighborhood = min_neighborhood[trim_size:-trim_size, trim_size:-trim_size]
    min_neighborhood_counts = min_neighborhood_counts[trim_size:-trim_size, trim_size:-trim_size]
    sc_local_den=min_neighborhood/min_neighborhood_counts
    sc_local_den[np.isnan(sc_local_den)] = 0
    sc_local_den[np.isinf(sc_local_den)] = 0
    return sc_local_den
    
#boundaries_record=sc_lis_pd_v2    
def cluster_result(boundaries_record,alpha=1):
    boundaries_record['strength*hdmd']=boundaries_record['strength']*boundaries_record['hdmd']
    boundaries_record['rank']=boundaries_record['strength*hdmd'].rank(ascending = False, method = 'dense')

    temp_rank = boundaries_record['rank'] / np.nanmax(boundaries_record['rank'])
    #print(candidates['eta'].describe())
    temp_eta = boundaries_record['strength*hdmd'] / np.nanmax(boundaries_record['strength*hdmd'])
    
    boundaries_record['transformed_strength*hdmd'] =  (temp_eta + temp_rank)/np.sqrt(2)
    
    breakpoint = boundaries_record.iloc[boundaries_record['transformed_strength*hdmd'].idxmin()]['strength*hdmd']
    breakpoint= alpha*breakpoint
    
    boundaries_record['cluster']=-1
    # boundaries_record.loc[ ((boundaries_record['strength*hdmd']>breakpoint) & (boundaries_record['hdmd']>np.sqrt(8))),\
    #                       'cluster'] = 1
    
    boundaries_record.loc[ (boundaries_record['strength*hdmd']>breakpoint),\
                          'cluster'] = 1
    
    #boundaries_record=boundaries_record.loc[boundaries_record['cluster']==1,:]
    return boundaries_record,breakpoint

#boundaries_record=sc_lis_pd    
# def cluster_result_v2(boundaries_record):
#     #boundaries_record['strength*hdmd']=boundaries_record['strength']*boundaries_record['hdmd']
#     boundaries_record['cluster']=-1
#     boundaries_record.loc[ ((boundaries_record['strength']>=1) & (boundaries_record['hdmd']>np.sqrt(50))) ,\
#                           'cluster'] = 1
#     return boundaries_record
def plot_dist_mat(temp_mat):
    temp_csc=sp.sparse.coo_matrix(temp_mat)
    temp_pd=pd.DataFrame({'i': temp_csc.row, 'j': temp_csc.col,'strength':temp_csc.data})
    sns.displot(temp_pd['strength'])

def identify_loops(All_strengthes,MAXR,DIS,upper_limit,binsize,Max_bin_distance,alpha=1):
    dist=Max_bin_distance
    #MAXR=100; dist=20000000;
    # strengthes_norm=strengthes.copy()
    # strengthes_norm=normlization_v2(strengthes_norm, dist, binsize
    #mycmap=generate_cmap([white1,red]);heatmap_2D(local_den_norm[100:200,100:200],mycmap)
    #strength_norm=normlization(strengthes,dist,binsize)
    strengthes_pad=np.pad(All_strengthes,((MAXR,MAXR),(MAXR,MAXR)),'constant',constant_values = (0,0))
    A=skimage.util.view_as_windows(strengthes_pad, (2*MAXR+1,2*MAXR+1), step=1)
    strengthes_pd=sp.sparse.coo_matrix(All_strengthes)
    sc_lis_pd=pd.DataFrame({'i': strengthes_pd.row+upper_limit, 'j': strengthes_pd.col+upper_limit,'strength':strengthes_pd.data})
    sc_lis_pd=sc_lis_pd[(sc_lis_pd['j']-sc_lis_pd['i'])<dist//binsize]
    #print(len(sc_lis_pd))
    #sns.displot(sc_lis_pd['strength']
    #A = A.reshape(A.shape[0],A.shape[0], 1,-1)
    locs_x=(sc_lis_pd['i']-upper_limit).tolist();
    locs_y=(sc_lis_pd['j']-upper_limit).tolist()
    
    #sc_lis_pd=pd.DataFrame(locs)
    hdmd_list=[]
    for m in range(len(locs_x)):
        #print(m)
        loc_x=locs_x[m];loc_y=locs_y[m]; value=All_strengthes[loc_x,loc_y]
        sub_A=A[loc_x,loc_y,:,:]
        sub_dis=(sub_A>value)*DIS
        temp=sub_dis[sub_dis>0]
        if len(temp)==0:
            hdmd_list.append(DIS[0,0])
        else:
            hdmd_list.append(np.min(temp))
    sc_lis_pd['hdmd']=hdmd_list
    
    sc_lis_pd_v2=sc_lis_pd.copy()
    sc_lis_pd_v2.index=range(len(sc_lis_pd_v2))
    
    sc_interactions_pd_raw,breakpoint=cluster_result(sc_lis_pd_v2,alpha)
    sc_interactions_pd=sc_interactions_pd_raw.loc[sc_interactions_pd_raw['cluster']==1,:]
    #return sc_interactions_pd,sc_interactions_pd_raw
    return sc_interactions_pd
    
#mat=strengthes.copy()
def minmax_array(mat):
    mat_norm=(mat-np.min(mat))/(np.max(mat)-np.min(mat))
    return mat_norm



# sub_mat=sub_matric;similarity=sub_similarity
def compute_significance_interaction_bulk(sub_mat,similarity,DIS,MAXR,Max_bin_distance,binsize,upper_limit=5,lower_limit=2,alpha=1,\
                                      Figure=False):
    #sub_mat = oe(sub_mat,expected=None)
    # from plot.plot_units import show_mat
    # show_mat(sub_mat,count=1)
    #strengthes,local_den=computing_interactions_strengthes(sub_mat,upper_limit,lower_limit)
    
    # strengthes_norm=strengthes.copy()
    # strengthes_norm=(expit(minmax_array(strengthes_norm)*5)-0.5)*2
    
    # # local_den_norm=local_den.copy()
    # # local_den_norm=(expit(minmax_array(local_den_norm)*5)-0.5)*2
    
    point_den=sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    point_den_norm=np.log10(1+np.abs(point_den)) * np.sign(point_den)
    
    #All_strengthes=np.multiply(point_den_norm,strengthes_norm)
    #All_strengthes=np.multiply(local_den_norm,strengthes_norm,point_den_norm)
    #All_strengthes=point_den_norm
    All_strengthes=np.multiply(point_den_norm,similarity[upper_limit:-upper_limit, upper_limit:-upper_limit])
    
    if Figure==True: 
        from plot.plot_units import show_mat
        start=575;end=775
        
        
        sub_matrix=sub_mat+sub_mat.T
        show_mat(sub_matrix[start:end,:][:,start:end],count=0)
        
    
        sub_matrix=All_strengthes+All_strengthes.T
        show_mat(sub_matrix[start:end,:][:,start:end],count=0)

        #start=500;end=800
        # sub_matrix=strengthes_norm+strengthes_norm.T
        # show_mat(sub_matrix[start:end,:][:,start:end],count=1)
        
        
    # from plot.plot_units import show_mat
    # show_mat(All_strengthes_aggr,count=1)
    #point_strengthes=(expit(point_strengthes)-0.5)*2
    #All_strengthes_aggr=np.multiply(All_strengthes_aggr,point_strengthes)
    aggr_interactions_pd=identify_loops(All_strengthes,MAXR,DIS,upper_limit,binsize,Max_bin_distance,alpha=0.5)
    
    if Figure==True:
        #show sc cells
        start=575;end=775
        #sc_interactions_pd['i']=sc_interactions_pd['i']
        #sc_interactions_pd['j']=sc_interactions_pd['j']
        temp=aggr_interactions_pd.loc[(aggr_interactions_pd['i']>=start)&(aggr_interactions_pd['j']>=start)&(aggr_interactions_pd['i']<end)\
                            &(aggr_interactions_pd['j']<end)]
        
        sub_matrix=sub_mat[start:end,start:end]
        sub_matrix[np.isnan(sub_matrix)] = 0
        
        #sub_matrix[sub_matrix>30]=30
        from plot.plot_units import show_mat_loops
        show_mat_loops(sub_matrix+sub_matrix.T,start,end,upper_limit,temp,count=1) 
    
    return aggr_interactions_pd


def get_localden_singlecell(sc_sub_mat,upper_limit):
    min_neighborhood = skimage.util.view_as_windows(sc_sub_mat, (2*1+1,2*1+1), step = 1)
    min_neighborhood = np.squeeze(min_neighborhood.reshape(min_neighborhood.shape[0],\
                                               min_neighborhood.shape[0], 1, -1))  
        
        
    min_neighborhood_counts = np.sum(~np.isnan(min_neighborhood), axis = -1)
    min_nonzero_neighborhood_counts = np.sum((~np.isnan(min_neighborhood)) & (min_neighborhood!=0), axis = -1)
    min_neighborhood = np.nansum(min_neighborhood,axis = -1)
    
    
    #remove edge cases that are used only as neighbors
    trim_size = upper_limit - 1
    min_neighborhood = min_neighborhood[trim_size:-trim_size, trim_size:-trim_size]
    min_neighborhood_counts = min_neighborhood_counts[trim_size:-trim_size, trim_size:-trim_size]
    
    min_nonzero_neighborhood_counts=min_nonzero_neighborhood_counts[trim_size:-trim_size, trim_size:-trim_size]
    
    local_den=min_neighborhood/min_neighborhood_counts
    local_den[min_nonzero_neighborhood_counts<=1]=0
    return local_den
        
from scipy.special import softmax    
#sub_mat=sub_matric;sc_sub_mat=sc_sub_matric;similarity=sub_similarity
def compute_significance_interaction(sub_mat,sc_sub_mat,similarity,DIS,MAXR,Max_bin_distance,binsize,upper_limit=5,lower_limit=2,alpha=0.5,\
                                     Figure=False):
    #computing interaction strengthes
    strengthes,local_den=computing_interactions_strengthes(sub_mat,upper_limit,lower_limit)
    #plot_dist_mat(strengthes); plot_dist_mat(local_den); plot_dist_mat(global_den)
    #strengthes_norm=strengthes.copy();strengthes_norm=normlization_v2(strengthes_norm, dist, binsize)
    #local_den_norm=strengthes.copy();local_den_norm=normlization_v2(local_den_norm, dist, binsize)
    #global_den_norm=strengthes.copy();global_den_norm=normlization_v2(global_den_norm, dist, binsize)

    #strengthes_norm=strengthes.copy(); strengthes_norm=(expit(strengthes_norm)-0.5)*2
    #plot_dist_mat(strengthes_norm)
    #point_strength=sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    # strengthes[strengthes!=0]
    # import matplotlib.pyplot as plt
    # plt.hist(strengthes[strengthes!=0])
    # import seaborn as sns
    # sns.distplot(strengthes[strengthes>1],rug=True)
    # try:
    #     strengthes[strengthes>np.percentile(strengthes[strengthes>1],95)]=0
    # except:
    #     passrank = rank; n_proc = n_proc; logger = logger
    strengthes_norm=strengthes.copy()
    strengthes_norm=(expit(minmax_array(strengthes_norm)*5)-0.5)*2
    
    local_den_norm=local_den.copy()
    local_den_norm=(expit(minmax_array(local_den_norm)*5)-0.5)*2
    
    # point_den=sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    # point_den_norm=(expit(minmax_array(point_den)*5)-0.5)*2
    
    #local_den_norm=minmax_array(local_den_norm)
    
    #Computing background supporting signal 
    All_strengthes=np.multiply(strengthes_norm,local_den_norm)
    # All_strengthes=(All_strengthes-All_strengthes.mean())/All_strengthes.std()
    # All_strengthes[All_strengthes<1.96]=0
    #All_strengthes=minmax_array(All_strengthes)
    All_strengthes_1=All_strengthes.copy()
    
    
    if Figure==True: 
        from plot.plot_units import show_mat
        start=1000;end=1200
        sub_matrix=sc_sub_mat+sc_sub_mat.T
        show_mat(sub_matrix[start:end,:][:,start:end],count=1)
        
        sub_matrix=sub_mat+sub_mat.T
        show_mat(sub_matrix[start:end,:][:,start:end],count=1)
        
    
        sub_matrix=All_strengthes+All_strengthes.T
        #sub_matrix=All_strengthes+All_strengthes.T
        show_mat(sub_matrix[start:end,:][:,start:end],count=1)
        
        
        #start=500;end=800
        sub_matrix=local_den_norm+local_den_norm.T
        show_mat(sub_matrix[start:end,:][:,start:end],count=1)
        
        #start=500;end=800
        sub_matrix=strengthes_norm+strengthes_norm.T
        from plot.plot_units import generate_cmap,red,blue,white1,red1,blue1,black,red2,red3
        mycmap=generate_cmap([white1,red,red2,red3,black])  
        show_mat(sub_matrix[start:end,:][:,start:end],mycmap,count=1)
        
        # from plot.plot_units import show_mat
        # show_mat(point_den[start:end,:][:,start:end],count=0)
        
    #similarity_norm=minmax_array(softmax(similarity))
    # import seaborn as sns
    # sns.distplot(similarity_norm[similarity_norm>0],rug=True)
    
    #np.percentile(local_den[local_den>0],95)
    # try:
    #     local_den[local_den<np.percentile(local_den[local_den>0],10)]=0
    # except:
    #     pass
    
    #Computing topological structure info from single cell
    #similarity_norm=similarity.copy()
    #similarity_norm=(expit(minmax_array(similarity_norm)*5)-0.5)*2
    #similarity_norm=softmax(similarity[upper_limit:-upper_limit, upper_limit:-upper_limit])
    All_strengthes=np.multiply(All_strengthes,similarity[upper_limit:-upper_limit, upper_limit:-upper_limit])
    #All_strengthes=All_strengthes+similarity[upper_limit:-upper_limit, upper_limit:-upper_limit]
    All_strengthes_2= All_strengthes.copy()
    
    
    if Figure==True:
        from plot.plot_units import show_mat,generate_cmap,red,blue,white1,red1,blue1
        mycmap=generate_cmap([blue1,white1,red1])
        show_mat(similarity[start:end,:][:,start:end],mycmap,count=0)
        #show_mat(similarity[start:end,:][:,start:end],mycmap,count=0)
        
    
    
    if Figure==True:
        #from plot.plot_units import show_mat
        
        from plot.plot_units import show_mat,generate_cmap,red,blue,white1,red1,blue1,black
        mycmap=generate_cmap([white1,red,black])
        
        sub_matrix=All_strengthes_2+All_strengthes_2.T
        show_mat(sub_matrix[start:end,:][:,start:end],mycmap,count=0)

    sc_point_strengthes=sc_sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    #+minmax_array(point_den)/100000
    All_strengthes=np.multiply(All_strengthes,sc_point_strengthes>0)
    #computing aggr loops
    All_strengthes_aggr=All_strengthes.copy()
    All_strengthes_aggr=minmax_array(All_strengthes_aggr)
    
    if Figure==True:
        #from plot.plot_units import show_mat
        
        from plot.plot_units import show_mat,generate_cmap,red,blue,white1,red1,blue1,black,red2,red3
        mycmap=generate_cmap([white1,red,red2,red3,black])
        
        
        sub_matrix=All_strengthes_aggr+All_strengthes_aggr.T
        show_mat(sub_matrix[start:end,:][:,start:end],mycmap,count=1)
        
        # sub_matrix=sc_point_strengthes+sc_point_strengthes.T
        # show_mat(sub_matrix[start:end,:][:,start:end]>0,count=0)
    #point_strengthes=sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    
    #All_strengthes_aggr[point_strengthes==0]=
    #point_strengthes=(expit(point_strengthes)-0.5)*2
    #All_strengthes_aggr=np.multiply(All_strengthes_aggr,point_strengthes)
    #sc_interactions_pd,sc_interactions_pd_raw=identify_loops(All_strengthes_aggr,MAXR,DIS,upper_limit,binsize,Max_bin_distance,alpha)
    sc_interactions_pd=identify_loops(All_strengthes_aggr,MAXR,DIS,upper_limit,binsize,Max_bin_distance,alpha)
    #sc_point_strengthes=sc_sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    #generate single cell loops
    # using_index=[ind for ind in aggr_interactions_pd.index if sc_point_strengthes[aggr_interactions_pd.loc[ind,'i']-upper_limit,aggr_interactions_pd.loc[ind,'j']-upper_limit]!=0]    
    # sc_interactions_pd=aggr_interactions_pd.loc[using_index]
    print(sc_interactions_pd.shape)
    
    
    if Figure==True:
        #show sc cells
        #start=800;end=1000
        #sc_interactions_pd['i']=sc_interactions_pd['i']
        #sc_interactions_pd['j']=sc_interactions_pd['j']
        temp=sc_interactions_pd.loc[(sc_interactions_pd['i']>=start)&(sc_interactions_pd['j']>=start)&(sc_interactions_pd['i']<end)\
                            &(sc_interactions_pd['j']<end)]
            
        
        from plot.plot_units import generate_cmap,red,blue,white1,red1,blue1,black,red2,red3
        mycmap=generate_cmap([white1,red,red2,red3,black])    
        
        sub_matrix=All_strengthes_aggr[start-upper_limit:end-upper_limit,start-upper_limit:end-upper_limit]
        sub_matrix[np.isnan(sub_matrix)] = 0
        from plot.plot_units import show_mat_loops
        show_mat_loops(sub_matrix+sub_matrix.T,start,end,upper_limit,temp,mycmap,count=0)
        
        # sub_matrix=All_strengthes_2[start-upper_limit:end-upper_limit,start-upper_limit:end-upper_limit]
        # sub_matrix[np.isnan(sub_matrix)] = 0
        # from plot.plot_units import show_mat_loops
        # show_mat_loops(sub_matrix+sub_matrix.T,start,end,upper_limit,temp,count=1)
        
        
        sub_matrix=similarity[start:end,start:end]
        sub_matrix[np.isnan(sub_matrix)] = 0
        from plot.plot_units import show_mat_loops
        
        from plot.plot_units import generate_cmap,red,blue,white1,red1,blue1
        mycmap=generate_cmap([blue1,white1,red1])
        
        show_mat_loops(sub_matrix+sub_matrix.T,start,end,upper_limit,temp,mycmap,count=0) 
        
        
        
        sub_matrix=sub_mat[start:end,start:end]
        sub_matrix[np.isnan(sub_matrix)] = 0
        from plot.plot_units import show_mat_loops
        show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=1) 
        
        
        
        sub_matrix=sc_sub_mat[start:end,start:end]
        sub_matrix[np.isnan(sub_matrix)] = 0
        from plot.plot_units import show_mat_loops
        show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=1) 
        
        
        
        
            
            
         
       
        #show aggr cells
        #sc_interactions_pd['i']=sc_interactions_pd['i']
        #sc_interactions_pd['j']=sc_interactions_pd['j']
        # temp=aggr_interactions_pd.loc[(aggr_interactions_pd['i']>=start)&(aggr_interactions_pd['j']>=start)&(aggr_interactions_pd['i']<end)\
        #                     &(aggr_interactions_pd['j']<end)]
        
        
        
        
        # sub_matrix=similarity[start:end,start:end]
        # sub_matrix[np.isnan(sub_matrix)] = 0
        # from plot.plot_units import show_mat_loops
        # show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=1) 
    
    
    #sc_interactions_pd_raw=boundaries_record
    #computing sc loops
    # All_strengthes_sc=All_strengthes.copy()
    # sc_point_strengthes=sc_sub_mat[upper_limit:-upper_limit, upper_limit:-upper_limit]
    # #point_strengthes=(expit(point_strengthes)-0.5)*2
    # #sc_strengthes=np.multiply(sc_local_den_norm,point_info)
    # #All_strengthes=np.multiply(All_strengthes,sc_local_den_norm)
    # All_strengthes_sc[sc_point_strengthes==0]=0
    # #point_strengthes=(expit(point_strengthes)-0.5)*2
    # #All_strengthes=np.multiply(All_strengthes,point_strengthes)
    # sc_interactions_pd=identify_loops(All_strengthes_sc,MAXR,DIS,upper_limit,binsize,Max_bin_distance,alpha)
    
    
    # if Figure==True:
        
    #     import matplotlib.pyplot as plt
    #     from mpl_toolkits.mplot3d import Axes3D
         
       
       
        
    #     sc_point_strengthes[sc_point_strengthes>0]=1
    #     #xs;    #single_cell strengthes
    #     #ys;    # strengthes
    #     #zs;    # local_den
    #     #ss;    # hdmd
    #     #clus   # clusters
        
    #     xs=[sc_point_strengthes[sc_interactions_pd_raw.loc[ind,'i']-upper_limit, sc_interactions_pd_raw.loc[ind,'j']-upper_limit] for ind in sc_interactions_pd_raw.index]
    #     ys=[strengthes_norm[sc_interactions_pd_raw.loc[ind,'i']-upper_limit, sc_interactions_pd_raw.loc[ind,'j']-upper_limit] for ind in sc_interactions_pd_raw.index]
    #     zs=[local_den_norm[sc_interactions_pd_raw.loc[ind,'i']-upper_limit, sc_interactions_pd_raw.loc[ind,'j']-upper_limit] for ind in sc_interactions_pd_raw.index]
        
    #     # ws=sc_interactions_pd_raw['strength']
    #     # ss=sc_interactions_pd_raw['hdmd']
    #     # clus=sc_interactions_pd_raw['cluster']
        
        
    #     sc_interactions_pd_raw['sc_v']=xs
    #     sc_interactions_pd_raw['stre_v']=ys
    #     sc_interactions_pd_raw['local_v']=zs
        
    #     from plot.plot_units import generate_cmap,color
    #     red1=color((207,0,25));blue1=color((0,42,111)); sha_bule=color((230,238,248));
    #     colors
        
    #     markers=['.','o']
    #     mycmap=generate_cmap([blue1,red1])
        
        
        
    #     sc_interactions_pd_raw.loc[sc_interactions_pd_raw['cluster']==-1,'cluster']=0
         
    # from plot.plot_units import red1,blue1,purple1    
    # colors=[red1,blue1];markers=['.','o']
    # ax=plt.figure(figsize=(4,4)).add_subplot(111)
    # for c in sc_interactions_pd_raw['cluster'].unique():
    #     color=colors[c]
    #     for sc in [0,1]:
    #         marker=markers[sc]
    #         indexes= (sc_interactions_pd_raw['cluster']==c) & (sc_interactions_pd_raw['sc_v']==sc)
    #         xs=sc_interactions_pd_raw.loc[indexes,'stre_v']
    #         ys=sc_interactions_pd_raw.loc[indexes,'local_v']
    #         #zs=sc_interactions_pd_raw.loc[indexes,'sc_v']
    #         if sc==0:
    #             ax.scatter(xs,ys,s=5,c=color,marker=marker)
    #         else:
    #             ax.scatter(xs,ys,s=15,c=color,marker=marker)
    
    # ax.set_xlabel('')
    # ax.set_ylabel('')
    
    # ax.grid(False)
    # plt.gcf()
    # ax.set_facecolor('w')
    
    # ax.set_xlim([0.01,1])
    # ax.set_ylim([0.07,1])
    
    # ax = plt.gca()
    # ax.spines['right'].set_color('none')
    # ax.spines['top'].set_color('none')
    # ax.xaxis.set_ticks_position('bottom')
    # ax.yaxis.set_ticks_position('left')
        
        
    #     # fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],"2d_decision",chrom, "pdf"]))
    #     # try:
    #     #     os.makedirs(os.path.join(outdir,'fig'))
    #     # except:
    #     #     pass 
    #     # #import matplotlib.pyplot as plt
    #     # plt.savefig(fig_address)
    #     # plt.close()
        
        
    #     ax=plt.figure(figsize=(4,4)).add_subplot(111,projection='3d')
    #     for c in sc_interactions_pd_raw['cluster'].unique():
    #         color=colors[c]
    #         for sc in [0,1]:
    #             marker=markers[sc]
    #             indexes= (sc_interactions_pd_raw['cluster']==c) & (sc_interactions_pd_raw['sc_v']==sc)
    #             xs=sc_interactions_pd_raw.loc[indexes,'stre_v']
    #             ys=sc_interactions_pd_raw.loc[indexes,'local_v']
    #             zs=sc_interactions_pd_raw.loc[indexes,'hdmd']
    #             if sc==0:
    #                 ax.scatter(xs,ys,zs,s=5,c=color,marker=marker)
    #             else:
    #                 ax.scatter(xs,ys,zs,s=15,c=color,marker=marker)
        
    #     #设置坐标轴
    #     ax.set_xlabel('')
    #     ax.set_ylabel('')
    #     ax.set_zlabel('')
        
    #     # ax.xaxis.pane.set_edgecolor('w')
    #     # ax.yaxis.pane.set_edgecolor('w')
    #     # ax.zaxis.pane.set_edgecolor('w')
        
    #     ax.grid(False)
        
        
    #     ax.w_xaxis.set_pane_color('w')
    #     ax.w_yaxis.set_pane_color('w')
    #     ax.w_zaxis.set_pane_color('w')
        
        
    #     # ax.w_xaxis.set_pane_color(sha_bule)
    #     # ax.w_yaxis.set_pane_color(sha_bule)
    #     # ax.w_zaxis.set_pane_color(sha_bule)
        
    #     ax.set_xticklabels([])
    #     ax.set_yticklabels([])
    #     ax.set_zticklabels([])
        
    #     # fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],"3d_decision",chrom, "pdf"]))
    #     # try:
    #     #     os.makedirs(os.path.join(outdir,'fig'))
    #     # except:
    #     #     pass 
    #     # #import matplotlib.pyplot as plt
    #     # plt.savefig(fig_address)
    #     # plt.close()
        
        
    #     #plot rank graph
        
    #     ax=plt.figure(figsize=(4,4)).add_subplot(111)
    #     for c in sc_interactions_pd_raw['cluster'].unique():
    #         color=colors[c]
    #         for sc in [0,1]:
    #             marker=markers[sc]
    #             indexes= (sc_interactions_pd_raw['cluster']==c) & (sc_interactions_pd_raw['sc_v']==sc)
    #             xs=sc_interactions_pd_raw.loc[indexes,'strength']
    #             ys=sc_interactions_pd_raw.loc[indexes,'hdmd']
    #             if sc==0:
    #                 ax.scatter(xs,ys,s=5,c=color,marker=marker)
    #             else:
    #                 ax.scatter(xs,ys,s=15,c=color,marker=marker)
    #     ax = plt.gca()
    #     ax.spines['right'].set_color('none')
    #     ax.spines['top'].set_color('none')
    #     ax.xaxis.set_ticks_position('bottom')
    #     ax.yaxis.set_ticks_position('left')
        
        
        
    #     # fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],"2d_decision_2",chrom, "pdf"]))
    #     # try:
    #     #     os.makedirs(os.path.join(outdir,'fig'))
    #     # except:
    #     #     pass 
    #     # #import matplotlib.pyplot as plt
    #     # plt.savefig(fig_address)
    #     # plt.close()
        
        
    #     ax=plt.figure(figsize=(4,4)).add_subplot(111)
    #     for c in sc_interactions_pd_raw['cluster'].unique():
    #         color=colors[c]
    #         for sc in [0,1]:
    #             marker=markers[sc]
    #             indexes= (sc_interactions_pd_raw['cluster']==c) & (sc_interactions_pd_raw['sc_v']==sc)
    #             xs=sc_interactions_pd_raw.loc[indexes,'rank']
    #             ys=sc_interactions_pd_raw.loc[indexes,'transformed_strength*hdmd']
    #             if sc==0:
    #                 ax.scatter(xs,ys,s=5,c=color,marker=marker)
    #             else:
    #                 ax.scatter(xs,ys,s=5,c=color,marker=marker)
        
        
        
    #     ax = plt.gca()
    #     ax.spines['right'].set_color('none')
    #     ax.spines['top'].set_color('none')
    #     ax.xaxis.set_ticks_position('bottom')
    #     ax.yaxis.set_ticks_position('left')
        
        
        
    #     # fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],"rank_decision",chrom, "pdf"]))
    #     # try:
    #     #     os.makedirs(os.path.join(outdir,'fig'))
    #     # except:
    #     #     pass 
    #     # #import matplotlib.pyplot as plt
    #     # plt.savefig(fig_address)
    #     # plt.close()
        
        
        
        
        
        
        
        
    #     ax=plt.figure(figsize=(4,4)).add_subplot(111)
    #     for c in sc_interactions_pd_raw['cluster'].unique():
    #         color=colors[c]
    #         for sc in [0,1]:
    #             marker=markers[sc]
    #             indexes= (sc_interactions_pd_raw['cluster']==c) & (sc_interactions_pd_raw['sc_v']==sc)
    #             xs=sc_interactions_pd_raw.loc[indexes,'rank']
    #             ys=sc_interactions_pd_raw.loc[indexes,'strength*hdmd']
    #             if sc==0:
    #                 ax.scatter(xs,ys,s=5,c=color,marker=marker)
    #             else:
    #                 ax.scatter(xs,ys,s=5,c=color,marker=marker)
                    
        
    #     ax = plt.gca()
    #     ax.spines['right'].set_color('none')
    #     ax.spines['top'].set_color('none')
    #     ax.xaxis.set_ticks_position('bottom')
    #     ax.yaxis.set_ticks_position('left')
        
        
        
    #     # ax.set_xlim([0,600])
    #     # ax.set_ylim([0,2.5])
        
    #     # fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],"rank_decision_2_local",chrom, "pdf"]))
    #     # try:
    #     #     os.makedirs(os.path.join(outdir,'fig'))
    #     # except:
    #     #     pass 
    #     # #import matplotlib.pyplot as plt
    #     # plt.savefig(fig_address)
    #     # # plt.close()
                    
    return sc_interactions_pd
                    
# def computing_allc_interactions(outdir,upper_limit=5, lower_limit=2):
#     # cell_lablels_add=os.path.join(prior_knowledge,'cell_4238_meta_cluster.txt')
#     # cell_info=pd.read_csv(cell_lablels_add,sep='\t',index_col=None,header=0)
#     # #get cell labels
#     # real_index=[cell_info['cell_id'].tolist().index(cell[11:]) for k, cell in enumerate(setnames)]
#     # cell_names=cell_info['clusters'].unique()
#     #distance table
#     MAXR=200
#     DIS=np.zeros((2*MAXR+1,2*MAXR+1))
#     for i in range(MAXR+1):
#         for j in range(MAXR+1):
#             dis=np.sqrt(i*i+j*j)
#             DIS[MAXR+i,MAXR+j]=dis;DIS[MAXR+i,MAXR-j]=dis;DIS[MAXR-i,MAXR+j]=dis;DIS[MAXR-i,MAXR-j]=dis
    
#     #computing global interaction strength
#     # mat_filenames=os.listdir(os.path.join(outdir, "TestData"))
#     # for ck, mat_filename in enumerate(mat_filenames):
#     #     mat_filename=mat_filenames[ck]
#     #     mat=np.load(os.path.join(outdir, "TestData",mat_filename))
#     #     mat = np.stack(mat,axis =-1)
#     #     if ck==0:
#     #         globle_sum_mat=np.sum(mat, axis = -1)
#     #     else:
#     #         globle_sum_mat+=np.sum(mat, axis = -1)
#     mat_filenames=os.listdir(os.path.join(outdir, "TestData"))
#     for ck, mat_filename in enumerate(mat_filenames):
#         mat_filename=mat_filenames[ck]
#         mat=np.load(os.path.join(outdir, "TestData",mat_filename))
#         mat = np.stack(mat,axis =-1)
#         cell_num=mat.shape[2]
#         for k in np.arange(0,cell_num,10):
#             sub_mat=mat[:,:,k:k+1].reshape((mat.shape[0],mat.shape[0]))
#             #globle_sum_mat=(mat>0).sum(axis=-1)/cell_num
#             #np.sum(mat, axis = -1)/cell_num
#             #extracting significant interactions of single cells
#             #comb_sc_interactions=compute_significance_interaction(sub_mat,upper_limit,lower_limit)

# def covert_pd_to_mat(df,binsize,mat_size):
#     df=df[df[6]>0]
#     df=df.iloc[:,[1,4,6]]
#     df.columns=["x1","y1",'value']
#     df['x1']=(df['x1']//binsize).astype(int)
#     df['y1']=(df['y1']//binsize).astype(int)
#     cell_mat = sp.sparse.csr_matrix( (df['value'], (df['x1'],df['y1'])), shape=(mat_size,mat_size)) 
#     cell_mat=cell_mat.toarray()
#     #start=500;end=2000
#     #sub_matrix=cell_mat[start:end,start:end]
#     #heatmap_2D(sub_matrix,mycmap)

def get_submatric_bulk(sub_mat,similarity,upper_limit,Max_bin_distance,binsize,mat_size):
    max_distance_bin = Max_bin_distance//binsize
    chrom_bins=sub_mat.shape[0]
    #mat_size more two times than max_distance_bin
    for i in range(0,chrom_bins,int(mat_size-max_distance_bin-upper_limit)):
        print(i)
        if (i+mat_size)<chrom_bins:
            sub_matric=sub_mat[i:i+mat_size,i:i+mat_size]
            sub_similarity=similarity[i:i+mat_size,i:i+mat_size]
            end_flag=0
            yield sub_matric,sub_similarity,i,end_flag
        else:
            sub_matric=sub_mat[i:chrom_bins,i:chrom_bins]
            #sc_sub_matric=sc_sub_mat[i:i+chrom_bins,i:i+chrom_bins]
            sub_similarity=similarity[i:chrom_bins,i:chrom_bins]
            end_flag=1
            yield sub_matric,sub_similarity,i,end_flag
            break;
        #yield sub_matric,sub_similarity,i,end_flag


    

def get_submatric(sub_mat,sc_sub_mat,similarity,upper_limit,Max_bin_distance,binsize,mat_size):
    max_distance_bin = Max_bin_distance//binsize
    chrom_bins=sub_mat.shape[0]
    #mat_size more two times than max_distance_bin
    end_flag=0
    for i in range(0,chrom_bins,int(mat_size-max_distance_bin-upper_limit)):
        print(i)
        if end_flag==0:
            if (i+mat_size)<chrom_bins:
                sub_matric=sub_mat[i:i+mat_size,i:i+mat_size]
                sc_sub_matric=sc_sub_mat[i:i+mat_size,i:i+mat_size]
                sub_similarity=similarity[i:i+mat_size,i:i+mat_size]
                end_flag=0
                yield sub_matric,sc_sub_matric,sub_similarity,i,end_flag
            else:
                sub_matric=sub_mat[i-100:chrom_bins,i-100:chrom_bins]
                #sc_sub_matric=sc_sub_mat[i:i+chrom_bins,i:i+chrom_bins]
                sc_sub_matric=sc_sub_mat[i-100:chrom_bins,i-100:chrom_bins]
                sub_similarity=similarity[i-100:chrom_bins,i-100:chrom_bins]
                end_flag=1
                yield sub_matric,sc_sub_matric,sub_similarity,i-100,end_flag
                break;
        #print(i, sub_matric.shape)
        



def bin_matrix(df, binsize):
    df.loc[:,'x1'] = df.loc[:,'x1'] // binsize
    df.loc[:,'y1'] = df.loc[:,'y1'] // binsize
    return df



import time
import time
def statics_time_multikernal(indir,outdir,chrom_lens, binsize, dist, neighborhood_limit_lower, \
            neighborhood_limit_upper,alpha=0.5,reso=100000, rank = 0, n_proc = 1, max_mem = 2, logger = None,Figure=False):
    #table for loop distance
    MAXR=100
    DIS=np.zeros((2*MAXR+1,2*MAXR+1))
    for i in range(MAXR+1):
        for j in range(MAXR+1):
            dis=np.sqrt(i*i+j*j)
            DIS[MAXR+i,MAXR+j]=dis;DIS[MAXR+i,MAXR-j]=dis;DIS[MAXR-i,MAXR+j]=dis;DIS[MAXR-i,MAXR-j]=dis
    
    upper_limit=neighborhood_limit_upper;
    lower_limit=neighborhood_limit_lower; 
    Max_bin_distance=dist
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    completed_filenames = os.listdir(indir)
    
    binsize=reso
    indir="/mnt/nas/Datasets/2020_Dip-C/GSE162511_"+str(int(binsize//1000))+'kb_knn=50/sampleaggrcell_dir_higashi_embed'
    mat_size=int((dist//binsize)*2+100)
    for chrom in proc_chroms:
        using_suffix=".".join([chrom, "sc", "npz"])
        sc_filenames=[os.path.join(indir,name) for name in completed_filenames if name.split('.',1)[1]==using_suffix]
        using_suffix2=".".join([chrom, "aggr", "npz"])
        #filenames=[name.rstrip(using_suffix)+'.'+using_suffix2 for name in sc_filenames] 
        filenames=[os.path.join(indir,re.sub("\.*"+using_suffix,'.'+using_suffix2,os.path.basename(name)) ) for name in sc_filenames]
        setnames = [os.path.basename(fname)[:-(len(using_suffix2)+1)] for fname in filenames]
        logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
        #k=1000;filename=filenames[k]  #1,2, 10,100,1000
        chr_time=[]
        t = time.perf_counter()
        chr_loop_number=[]
        for k, filename in enumerate(filenames[0:1812:364]):
            #address=os.path.join(outdir1, ".".join([setnames[k],"interactions",chrom, "csv"]))
            #sc_address=os.path.join(outdir, ".".join([setnames[k],'sc',"interactions",chrom, "csv"]))
            #raw_address=os.path.join(outdir, ".".join([setnames[k],'sc_raw',"interactions",chrom, "csv"]))
            print( filename)
            #sub_mat=np.zeros((NUM,NUM))
            csr_mat = sparse.load_npz(filename)
            print("the chromosome: %s the th%d cell:, the mat size %d" % (chrom,k,csr_mat.shape[0]) )
            sub_mat = csr_mat.toarray()
            #sub_mat=sqrt_norm(sub_mat+sub_mat.T)
            sub_mat = sqrt_norm(sub_mat+sub_mat.T)
            #sub_mat=np.log2(1+np.abs(sub_mat))*np.sign(sub_mat)
            sub_mat = oe(sub_mat,expected=None)
            #show_mat(sub_mat,count=1)
            #read sc mat
            csr_mat = sparse.load_npz(sc_filenames[k])
            sc_sub_mat = csr_mat.toarray()
            sc_sub_mat = sc_sub_mat + sc_sub_mat.T
            sc_sub_mat=sc_sub_mat+np.eye(sc_sub_mat.shape[0])
            
            #Computing bin's embedding in single cell interaction networks
            g = nx.from_numpy_array(sc_sub_mat,create_using = nx.Graph())
            #g = nx.from_numpy_array(sub_mat,create_using = nx.Graph())
            model = Node2Vec(g, dimensions=32,walk_length=3,num_walks=50,workers=4,p=1,q=1) 
            #model = Node2Vec(g, dimensions=32,walk_length=10,num_walks=100,workers=4,p=1,q=1,weight_key='weight') 
            # Use temp_folder for big graphs
            model = model.fit(window=3)  # Any keywords acceptable by gensim.Word2Vec can be passed, `dimensions` and `workers` are automatically passed (from the Node2Vec constructor)
            #model = model.fit(window=5)
            #model.wv.most_similar('200.0')
            embeddings = {}
            for word in g.nodes():
                embeddings[word] = model.wv[str(word)] 
            similarity=cosine_similarity(pd.DataFrame(embeddings).T)
            
            del g
            gc.collect()
            
            #sub_mat=sqrt_norm(sub_mat+sub_mat.T)
            #sc_sub_mat = sqrt_norm(sc_sub_mat+sc_sub_mat.T)
            #sub_mat=np.log2(1+np.abs(sub_mat))*np.sign(sub_mat)
            #show_mat(sc_sub_mat)
            #sc_sub_mat = oe(sc_sub_mat,expected=None)
            add_n=sub_mat.shape[0]-sc_sub_mat.shape[0]
            if add_n>0:
                sc_sub_mat=np.pad(sc_sub_mat,((0,add_n),(0,add_n)),'constant',constant_values = (0,0))
                similarity=np.pad(similarity,((0,add_n),(0,add_n)),'constant',constant_values = (0,0))
            
            if sub_mat.shape[0]<=mat_size:
                #sub_mat=sub_mat[:,600:1800][600:1800,:];sc_sub_mat=sc_sub_mat[:,600:1800][600:1800,:]
                sc_interactions_pd_a=compute_significance_interaction(sub_mat,sc_sub_mat,similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                    upper_limit,\
                                                                    lower_limit,alpha) 
            else:
                submatrices=get_submatric(sub_mat,sc_sub_mat,similarity,upper_limit,Max_bin_distance,binsize,mat_size)
                for i,(sub_matric,sc_sub_matric,sub_similarity,start_index,end_flag) in enumerate(submatrices):
                    print(i,sub_matric.shape,end_flag)
                    # if i==0:
                    #     break
                    try:
                        if sc_sub_matric.sum()>np.diagonal(sc_sub_matric).sum():
                            
                            sc_interactions_pd=compute_significance_interaction(sub_matric,sc_sub_matric,sub_similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                            upper_limit,\
                                                                            lower_limit,alpha) 
                            sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']>upper_limit) & (sc_interactions_pd['j']>upper_limit)]
                            #sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']>upper_limit) & (sc_interactions_pd_raw['j']>upper_limit)]   
                            
                            if end_flag==0:
                                sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd['j']<(mat_size-Max_bin_distance//binsize))]
                                #sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd_raw['j']<(mat_size-Max_bin_distance//binsize))]
                            else:
                                sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']>upper_limit+100) & (sc_interactions_pd['j']>upper_limit+100)]
                        else:
                                sc_interactions_pd=pd.DataFrame()
                                
                        if start_index==0:
                            sc_interactions_pd_a=sc_interactions_pd
                            #sc_interactions_pd_raw_a=sc_interactions_pd_raw
                        else:
                            sc_interactions_pd['i']=sc_interactions_pd['i']+start_index
                            sc_interactions_pd['j']=sc_interactions_pd['j']+start_index
                            #sc_interactions_pd_raw['i']=sc_interactions_pd_raw['i']+start_index
                            #sc_interactions_pd_raw['j']=sc_interactions_pd_raw['j']+start_index
                            sc_interactions_pd_a=pd.concat([sc_interactions_pd_a,sc_interactions_pd])
                            #sc_interactions_pd_raw_a=pd.concat([sc_interactions_pd_raw_a,sc_interactions_pd_raw])
                    except:
                        sc_interactions_pd_a=pd.DataFrame([])
                        pass



def statics_time(indir,outdir,chrom_lens, binsize, dist, neighborhood_limit_lower, \
            neighborhood_limit_upper,alpha=0.5, rank = 0, n_proc = 1, max_mem = 2, logger = None,Figure=False):
    
    resolutions=[10000,20000,50000,100000]
    
    logger.set_rank(rank)
    try:
        os.makedirs(outdir)
    except:
        pass
    #table for loop distance
    MAXR=100
    DIS=np.zeros((2*MAXR+1,2*MAXR+1))
    for i in range(MAXR+1):
        for j in range(MAXR+1):
            dis=np.sqrt(i*i+j*j)
            DIS[MAXR+i,MAXR+j]=dis;DIS[MAXR+i,MAXR-j]=dis;DIS[MAXR-i,MAXR+j]=dis;DIS[MAXR-i,MAXR-j]=dis
    
    upper_limit=neighborhood_limit_upper;
    lower_limit=neighborhood_limit_lower; 
    Max_bin_distance=dist
    proc_chroms = get_proc_chroms(chrom_lens, rank, 1)
    completed_filenames = os.listdir(indir)
    
    all_time={}
    all_loop_number={}
    for reso in resolutions:
        binsize=reso
        indir="/mnt/nas/Datasets/2020_Dip-C/GSE162511_"+str(int(binsize//1000))+'kb_knn=50/sampleaggrcell_dir_higashi_embed'
        
        mat_size=int((dist//binsize)*2+100)
        
        for chrom in proc_chroms:
            using_suffix=".".join([chrom, "sc", "npz"])
            sc_filenames=[os.path.join(indir,name) for name in completed_filenames if name.split('.',1)[1]==using_suffix]
            using_suffix2=".".join([chrom, "aggr", "npz"])
            #filenames=[name.rstrip(using_suffix)+'.'+using_suffix2 for name in sc_filenames] 
            filenames=[os.path.join(indir,re.sub("\.*"+using_suffix,'.'+using_suffix2,os.path.basename(name)) ) for name in sc_filenames]
            setnames = [os.path.basename(fname)[:-(len(using_suffix2)+1)] for fname in filenames]
            logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
            #k=1000;filename=filenames[k]  #1,2, 10,100,1000
            chr_time=[]
            t = time.perf_counter()
            chr_loop_number=[]
            for k, filename in enumerate(filenames[0:1812:364]):
                #address=os.path.join(outdir1, ".".join([setnames[k],"interactions",chrom, "csv"]))
                #sc_address=os.path.join(outdir, ".".join([setnames[k],'sc',"interactions",chrom, "csv"]))
                #raw_address=os.path.join(outdir, ".".join([setnames[k],'sc_raw',"interactions",chrom, "csv"]))
                print( filename)
                #sub_mat=np.zeros((NUM,NUM))
                csr_mat = sparse.load_npz(filename)
                print("the chromosome: %s the th%d cell:, the mat size %d" % (chrom,k,csr_mat.shape[0]) )
                sub_mat = csr_mat.toarray()
                #sub_mat=sqrt_norm(sub_mat+sub_mat.T)
                sub_mat = sqrt_norm(sub_mat+sub_mat.T)
                #sub_mat=np.log2(1+np.abs(sub_mat))*np.sign(sub_mat)
                sub_mat = oe(sub_mat,expected=None)
                #show_mat(sub_mat,count=1)
                #read sc mat
                csr_mat = sparse.load_npz(sc_filenames[k])
                sc_sub_mat = csr_mat.toarray()
                sc_sub_mat = sc_sub_mat + sc_sub_mat.T
                sc_sub_mat=sc_sub_mat+np.eye(sc_sub_mat.shape[0])
                
                #Computing bin's embedding in single cell interaction networks
                g = nx.from_numpy_array(sc_sub_mat,create_using = nx.Graph())
                #g = nx.from_numpy_array(sub_mat,create_using = nx.Graph())
                model = Node2Vec(g, dimensions=32,walk_length=3,num_walks=50,workers=4,p=1,q=1) 
                #model = Node2Vec(g, dimensions=32,walk_length=10,num_walks=100,workers=4,p=1,q=1,weight_key='weight') 
                # Use temp_folder for big graphs
                model = model.fit(window=3)  # Any keywords acceptable by gensim.Word2Vec can be passed, `dimensions` and `workers` are automatically passed (from the Node2Vec constructor)
                #model = model.fit(window=5)
                #model.wv.most_similar('200.0')
                embeddings = {}
                for word in g.nodes():
                    embeddings[word] = model.wv[str(word)] 
                similarity=cosine_similarity(pd.DataFrame(embeddings).T)
                
                del g
                gc.collect()
                
                #sub_mat=sqrt_norm(sub_mat+sub_mat.T)
                #sc_sub_mat = sqrt_norm(sc_sub_mat+sc_sub_mat.T)
                #sub_mat=np.log2(1+np.abs(sub_mat))*np.sign(sub_mat)
                #show_mat(sc_sub_mat)
                #sc_sub_mat = oe(sc_sub_mat,expected=None)
                add_n=sub_mat.shape[0]-sc_sub_mat.shape[0]
                if add_n>0:
                    sc_sub_mat=np.pad(sc_sub_mat,((0,add_n),(0,add_n)),'constant',constant_values = (0,0))
                    similarity=np.pad(similarity,((0,add_n),(0,add_n)),'constant',constant_values = (0,0))
                
                if sub_mat.shape[0]<=mat_size:
                    #sub_mat=sub_mat[:,600:1800][600:1800,:];sc_sub_mat=sc_sub_mat[:,600:1800][600:1800,:]
                    sc_interactions_pd_a=compute_significance_interaction(sub_mat,sc_sub_mat,similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                        upper_limit,\
                                                                        lower_limit,alpha) 
                else:
                    submatrices=get_submatric(sub_mat,sc_sub_mat,similarity,upper_limit,Max_bin_distance,binsize,mat_size)
                    for i,(sub_matric,sc_sub_matric,sub_similarity,start_index,end_flag) in enumerate(submatrices):
                        print(i,sub_matric.shape,end_flag)
                        # if i==0:
                        #     break
                        try:
                            if sc_sub_matric.sum()>np.diagonal(sc_sub_matric).sum():
                                
                                sc_interactions_pd=compute_significance_interaction(sub_matric,sc_sub_matric,sub_similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                                upper_limit,\
                                                                                lower_limit,alpha) 
                                sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']>upper_limit) & (sc_interactions_pd['j']>upper_limit)]
                                #sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']>upper_limit) & (sc_interactions_pd_raw['j']>upper_limit)]   
                                
                                if end_flag==0:
                                    sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd['j']<(mat_size-Max_bin_distance//binsize))]
                                    #sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd_raw['j']<(mat_size-Max_bin_distance//binsize))]
                                else:
                                    sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']>upper_limit+100) & (sc_interactions_pd['j']>upper_limit+100)]
                            else:
                                    sc_interactions_pd=pd.DataFrame()
                                    
                            if start_index==0:
                                sc_interactions_pd_a=sc_interactions_pd
                                #sc_interactions_pd_raw_a=sc_interactions_pd_raw
                            else:
                                sc_interactions_pd['i']=sc_interactions_pd['i']+start_index
                                sc_interactions_pd['j']=sc_interactions_pd['j']+start_index
                                #sc_interactions_pd_raw['i']=sc_interactions_pd_raw['i']+start_index
                                #sc_interactions_pd_raw['j']=sc_interactions_pd_raw['j']+start_index
                                sc_interactions_pd_a=pd.concat([sc_interactions_pd_a,sc_interactions_pd])
                                #sc_interactions_pd_raw_a=pd.concat([sc_interactions_pd_raw_a,sc_interactions_pd_raw])
                        except:
                            sc_interactions_pd_a=pd.DataFrame([])
                            pass
                    chr_time.append(time.perf_counter()-t)
                    chr_loop_number.append(len(sc_interactions_pd_a))
                    
                    print(f'coast:{time.perf_counter()-t:.4f}s')
                    t = time.perf_counter()
                all_time[str(reso)+'_'+chrom]=chr_time
                all_loop_number[str(reso)+'_'+chrom]=chr_loop_number
    
    
    
    
    if Figure==True:
        
        #all_time; all_loop_number
        k=0
        A_PD=pd.DataFrame()
        for reso in resolutions:
            for chrom in proc_chroms:
                temp_pd=pd.DataFrame([all_time[str(reso)+'_'+chrom],all_loop_number[str(reso)+'_'+chrom]])
                temp_pd=temp_pd.T
                temp_pd['reso']=reso
                temp_pd['chrom']=chrom
                if k==0:
                    A_PD=temp_pd
                else:
                    A_PD=pd.concat([A_PD,temp_pd],axis=0)
                k+=1
                
                
        
        from plot.plot_units import colors_diff,colors_nature,color_3
        A_PD.columns=['time','loopnumber','reso','chrom']
        #svae pd
        time_analysis_dir='/mnt/nas/Datasets/2020_Dip-C/time_analyasis'
        A_PD.to_csv(time_analysis_dir+'/singlekernal.csv',sep='\t',header=True)
        
        fig,ax=plt.subplots(figsize=(14,4))
        sns.barplot(A_PD,x='chrom',y='time',hue='reso',palette=colors_nature,alpha=0.9)
        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylim(0)
        plt.show()   
        
        
        #compared singlekernal and multiplekernal
        for reso in resolutions:
            for chrom in proc_chroms:
                temp_pd1=A_PD[(A_PD['reso']==reso) & (A_PD['chrom']==chrom)]
                print(reso,chrom,np.mean(temp_pd1['time']))
                
        
        
        
        #read_BPD
        B_PD=pd.read_csv(time_analysis_dir+'/multiplekernal.csv',sep='\t',header=0,index_col=0)
        
        Com_PD=pd.DataFrame()
        for k,reso in enumerate(resolutions):
            
            temp_pd1=A_PD[A_PD['reso']==reso]
            temp_pd2=B_PD[B_PD['reso']==reso]
            
            temp1=pd.DataFrame([np.mean(temp_pd1['time']),reso,'singlekernal'])
            
            temp2=pd.DataFrame([(temp_pd2['time'].values[0]/len(temp_pd1)),reso,'multikernal'])
            
            
            print(temp1.loc[0]/temp2.loc[0])
            Temp=pd.concat([temp1.T,temp2.T],axis=0)
            
            
            
            if k==0:
                Com_PD=Temp
            else:
                Com_PD=pd.concat([Com_PD,Temp],axis=0)
        
        
        Com_PD.columns=['time','reso','proc_num']
        fig,ax=plt.subplots(figsize=(4,4))
        sns.barplot(Com_PD,x='reso',y='time',hue='proc_num',palette=color_3)
        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.set_ylim(0)
        plt.show()   
        
        
            
            


import networkx as nx
from node2vec import Node2Vec
from sklearn.metrics.pairwise import cosine_similarity
def scloops(indir,outdir,chrom_lens, binsize, dist, neighborhood_limit_lower, \
            neighborhood_limit_upper,alpha=0.5, rank = 0, n_proc = 1, max_mem = 2, logger = None,Figure=False):        
    
    logger.set_rank(rank)
    try:
        os.makedirs(outdir)
    except:
        pass
    #table for loop distance
    MAXR=100
    DIS=np.zeros((2*MAXR+1,2*MAXR+1))
    for i in range(MAXR+1):
        for j in range(MAXR+1):
            dis=np.sqrt(i*i+j*j)
            DIS[MAXR+i,MAXR+j]=dis;DIS[MAXR+i,MAXR-j]=dis;DIS[MAXR-i,MAXR+j]=dis;DIS[MAXR-i,MAXR-j]=dis
    
    upper_limit=neighborhood_limit_upper;
    lower_limit=neighborhood_limit_lower; 
    Max_bin_distance=dist
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    completed_filenames = os.listdir(indir)
    
    #chrom='chr1'
    if binsize<=50000:
        mat_size=np.max([int((dist//binsize)*2+100),1200])
    else:
        mat_size=np.max([int((dist//binsize)*2+100),300]) 
        
    #statics_time()
    #t = time.time()
    for chrom in proc_chroms:
        using_suffix=".".join([chrom, "sc", "npz"])
        sc_filenames=[os.path.join(indir,name) for name in completed_filenames if name.split('.',1)[1]==using_suffix]
        using_suffix2=".".join([chrom, "aggr", "npz"])
        #filenames=[name.rstrip(using_suffix)+'.'+using_suffix2 for name in sc_filenames] 
        filenames=[os.path.join(indir,re.sub("\.*"+using_suffix,'.'+using_suffix2,os.path.basename(name)) ) for name in sc_filenames]
        setnames = [os.path.basename(fname)[:-(len(using_suffix2)+1)] for fname in filenames]
        logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
        #k=1000;filename=filenames[k]  #1,2, 10,100,1000
        for k, filename in enumerate(filenames):
            #address=os.path.join(outdir1, ".".join([setnames[k],"interactions",chrom, "csv"]))
            sc_address=os.path.join(outdir, ".".join([setnames[k],'sc',"interactions",chrom, "csv"]))
            #raw_address=os.path.join(outdir, ".".join([setnames[k],'sc_raw',"interactions",chrom, "csv"]))
            if os.path.exists(sc_address):
                #continue
                sc_interactions_pd_a=pd.read_csv(sc_address,sep = ",",index_col = None,header=0) 
                #sc_interactions_pd_raw=pd.read_csv(raw_address,sep = ",",index_col = None,header=0) 
            else:
                print( filename)
                #sub_mat=np.zeros((NUM,NUM))
                csr_mat = sparse.load_npz(filename)
                print("the chromosome: %s the th%d cell:, the mat size %d" % (chrom,k,csr_mat.shape[0]) )
                sub_mat = csr_mat.toarray()
                #sub_mat=sqrt_norm(sub_mat+sub_mat.T)
                sub_mat = sqrt_norm(sub_mat+sub_mat.T)
                #sub_mat=np.log2(1+np.abs(sub_mat))*np.sign(sub_mat)
                sub_mat = oe(sub_mat,expected=None)
                #show_mat(sub_mat,count=1)
                #read sc mat
                csr_mat = sparse.load_npz(sc_filenames[k])
                sc_sub_mat = csr_mat.toarray()
                sc_sub_mat = sc_sub_mat + sc_sub_mat.T
                sc_sub_mat=sc_sub_mat+np.eye(sc_sub_mat.shape[0])
                
                #Computing bin's embedding in single cell interaction networks
                g = nx.from_numpy_array(sc_sub_mat,create_using = nx.Graph())
                #g = nx.from_numpy_array(sub_mat,create_using = nx.Graph())
                model = Node2Vec(g, dimensions=32,walk_length=3,num_walks=50,workers=4,p=1,q=1) 
                #model = Node2Vec(g, dimensions=32,walk_length=10,num_walks=100,workers=4,p=1,q=1,weight_key='weight') 
                # Use temp_folder for big graphs
                model = model.fit(window=3)  # Any keywords acceptable by gensim.Word2Vec can be passed, `dimensions` and `workers` are automatically passed (from the Node2Vec constructor)
                #model = model.fit(window=5)
                #model.wv.most_similar('200.0')
                embeddings = {}
                for word in g.nodes():
                    embeddings[word] = model.wv[str(word)] 
                similarity=cosine_similarity(pd.DataFrame(embeddings).T)
                # from scipy.special import softmax
                # similarity=softmax(similarity)
                # thr=np.percentile(similarity,80)
                # similarity[similarity<thr]=0
                if Figure==True:
                    from plot.plot_units import show_mat
                    #start=4000;end=4250
                    start=1000;end=1200
                    sub_matrix=similarity[start:end,start:end]
                    sub_matrix[np.isnan(sub_matrix)] = 0
                    from plot.plot_units import generate_cmap,red,blue,white1,red1,blue1
                    mycmap=generate_cmap([blue1,white1,red1])
                    show_mat(sub_matrix,mycmap,count=0)
                    
                    #show_mat(np.corrcoef(sub_matrix),mycmap,count=0)
                    
                    # show_mat(sub_mat[start:end,:][:,start:end],count=1)
                    # show_mat(sc_sub_mat[start:end,:][:,start:end],count=1)
                #print('deleting edgelist', setname, chrom)
                #sys.stdout.flush()
                del g
                gc.collect()
                
                #sub_mat=sqrt_norm(sub_mat+sub_mat.T)
                #sc_sub_mat = sqrt_norm(sc_sub_mat+sc_sub_mat.T)
                #sub_mat=np.log2(1+np.abs(sub_mat))*np.sign(sub_mat)
                #show_mat(sc_sub_mat)
                #sc_sub_mat = oe(sc_sub_mat,expected=None)
                add_n=sub_mat.shape[0]-sc_sub_mat.shape[0]
                if add_n>0:
                    sc_sub_mat=np.pad(sc_sub_mat,((0,add_n),(0,add_n)),'constant',constant_values = (0,0))
                    similarity=np.pad(similarity,((0,add_n),(0,add_n)),'constant',constant_values = (0,0))
                
                # if Figure==True:
                #     #show sc cells
                #     start=800;end=1000
                #     sub_matrix=sc_sub_mat[start:end,start:end]
                #     sub_matrix[np.isnan(sub_matrix)] = 0

                #     from plot.plot_units import show_mat,generate_cmap,color,white1
                #     red1=color((207,0,25));blue1=color((0,42,111));
                    
                #     #sub_matrix[sub_matrix>0]=1
                #     show_mat(sub_matrix,count=1) 
                    
                #     sc_fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],'sc',"heatmap",chrom, "pdf"]))
                #     try:
                #         os.makedirs(os.path.join(outdir,'fig'))
                #     except:
                #         pass 
                #     #import matplotlib.pyplot as plt
                #     plt.savefig(sc_fig_address)
                #     plt.close()
                    
                #     sub_matrix=sub_mat[start:end,start:end]
                #     sub_matrix[np.isnan(sub_matrix)] = 0
                #     show_mat(sub_matrix,generate_cmap([white1,red1]),count=1) 
                    
                    
                #     aggr_fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],'aggr',"heatmap",chrom, "pdf"]))
                #     try:
                #         os.makedirs(os.path.join(outdir,'fig'))
                #     except:
                #         pass 
                #     #import matplotlib.pyplot as plt
                #     plt.savefig(aggr_fig_address)
                #     plt.close()
                #show_mat(sub_mat)
                #show_mat(sub_mat,count=1)
                #extracting significant interactions of single cells
                if sub_mat.shape[0]<=mat_size:
                    #sub_mat=sub_mat[:,600:1800][600:1800,:];sc_sub_mat=sc_sub_mat[:,600:1800][600:1800,:]
                    sc_interactions_pd_a=compute_significance_interaction(sub_mat,sc_sub_mat,similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                        upper_limit,\
                                                                        lower_limit,alpha) 
                else:
                    submatrices=get_submatric(sub_mat,sc_sub_mat,similarity,upper_limit,Max_bin_distance,binsize,mat_size)
                    for i,(sub_matric,sc_sub_matric,sub_similarity,start_index,end_flag) in enumerate(submatrices):
                        print(i,sub_matric.shape,end_flag)
                        # if i==0:
                        #     break
                        try:
                            if sc_sub_matric.sum()>np.diagonal(sc_sub_matric).sum():
                                
                                sc_interactions_pd=compute_significance_interaction(sub_matric,sc_sub_matric,sub_similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                                upper_limit,\
                                                                                lower_limit,alpha) 
                                sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']>upper_limit) & (sc_interactions_pd['j']>upper_limit)]
                                #sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']>upper_limit) & (sc_interactions_pd_raw['j']>upper_limit)]   
                                
                                if end_flag==0:
                                    sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd['j']<(mat_size-Max_bin_distance//binsize))]
                                    #sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd_raw['j']<(mat_size-Max_bin_distance//binsize))]
                                else:
                                    sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']>upper_limit+100) & (sc_interactions_pd['j']>upper_limit+100)]
                            else:
                                    sc_interactions_pd=pd.DataFrame()
                                    
                            if start_index==0:
                                sc_interactions_pd_a=sc_interactions_pd
                                #sc_interactions_pd_raw_a=sc_interactions_pd_raw
                            else:
                                sc_interactions_pd['i']=sc_interactions_pd['i']+start_index
                                sc_interactions_pd['j']=sc_interactions_pd['j']+start_index
                                #sc_interactions_pd_raw['i']=sc_interactions_pd_raw['i']+start_index
                                #sc_interactions_pd_raw['j']=sc_interactions_pd_raw['j']+start_index
                                sc_interactions_pd_a=pd.concat([sc_interactions_pd_a,sc_interactions_pd])
                                #sc_interactions_pd_raw_a=pd.concat([sc_interactions_pd_raw_a,sc_interactions_pd_raw])
                        except:
                            sc_interactions_pd_a=pd.DataFrame([])
                            pass
                
                
                #sc_interactions_pd=sc_interactions_pd_a 
                #sc_interactions_pd_raw=sc_interactions_pd_raw_a
                if len(sc_interactions_pd_a)!=0:
                    sc_interactions_pd=sc_interactions_pd_a[['i','j','strength']]
                    #sc_interactions_pd['value']=[sub_mat[i,j] for i,j in zip(sc_interactions_pd['i'].tolist(),sc_interactions_pd['j'].tolist())]
                    sc_interactions_pd.to_csv(sc_address,sep = ",",index = False,header=True)  
                # else:
                #     #sc_interactions_pd.to_csv(sc_address,sep = ",",index = False,header=True) 
                #     pd.DataFrame(['i','j','strength']).T.to_csv(sc_address,sep = ",",index = False,header=False)
                
                if Figure==True:
                    from plot.plot_units import show_mat
                    start=1000;end=1400
                    #start=300;end=500
                    sub_matrix=similarity[start:end,start:end]
                    sub_matrix[np.isnan(sub_matrix)] = 0
                    from plot.plot_units import generate_cmap,red,blue,white1,red1,blue1
                    mycmap=generate_cmap([blue1,white1,red1])
                    
                    show_mat(sub_matrix,mycmap,count=0)
                    
    
                    # show_mat(sub_mat[start:end,:][:,start:end],count=1)
                    # show_mat(sc_sub_mat[start:end,:][:,start:end],count=1)
                    
                    
                    #sub_matrix=sub_mat[start:end,:][:,start:end]
                    #show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=1)  
                    #sub_matrix=sc_sub_mat[start:end,:][:,start:end]
                    #show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=1)
                    
                    
                    
                    #show sc cells
                    print(len(sc_interactions_pd))
                    #start=4550;end=4750
                    #start=1000;end=1200
                    #sc_interactions_pd['i']=sc_interactions_pd['i']
                    #sc_interactions_pd['j']=sc_interactions_pd['j']
                    temp=sc_interactions_pd.loc[(sc_interactions_pd['i']>=start)&(sc_interactions_pd['j']>=start)&(sc_interactions_pd['i']<end)\
                                        &(sc_interactions_pd['j']<end)]
                        
                        
                        
                    
                    from plot.plot_units import show_mat_loops
                    # show_mat_loops(sub_matrix,start,end,upper_limit,temp,mycmap,count=1) 
                    
                    
                    # from plot.plot_units import generate_cmap,red,blue,white1,red1,blue1,black,red2,red3
                    # mycmap=generate_cmap([white1,red,red2,red3,black])    
                   
                    
                    # sub_matrix=All_strengthes_aggr[start-upper_limit:end-upper_limit,start-upper_limit:end-upper_limit]
                    # sub_matrix[np.isnan(sub_matrix)] = 0
                    
                    # show_mat(sub_matrix+sub_matrix.T,mycmap,count=1)
                    
                    # from plot.plot_units import show_mat_loops
                    # show_mat_loops(sub_matrix+sub_matrix.T,start,end,upper_limit,temp,mycmap,count=1)
                    
                    
                        
                     
                        
                    sub_matrix=sc_sub_mat[start:end,start:end]
                    sub_matrix[np.isnan(sub_matrix)] = 0
                    
                    show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=1)  
                    
                    
                   
                    sc_fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],'sc',"loops",chrom, "pdf"]))
                    try:
                        os.makedirs(os.path.join(outdir,'fig'))
                    except:
                        pass 
                    #import matplotlib.pyplot as plt
                    plt.savefig(sc_fig_address)
                    plt.close()
                   
                
                    #show aggr cells
                    print(len(sc_interactions_pd))
                    #sc_interactions_pd['i']=sc_interactions_pd['i']
                    #sc_interactions_pd['j']=sc_interactions_pd['j']
                    # temp=aggr_interactions_pd.loc[(aggr_interactions_pd['i']>=start)&(aggr_interactions_pd['j']>=start)&(aggr_interactions_pd['i']<end)\
                    #                     &(aggr_interactions_pd['j']<end)]
                    sub_matrix=sub_mat[start:end,start:end]
                    sub_matrix[np.isnan(sub_matrix)] = 0
                    from plot.plot_units import show_mat_loops
                    show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=1) 
                    
                    
                    aggr_fig_address=os.path.join(outdir,'fig', ".".join([setnames[k],'aggr',"loops",chrom, "pdf"]))
                    try:
                        os.makedirs(os.path.join(outdir,'fig'))
                    except:
                        pass 
                    #import matplotlib.pyplot as plt
                    plt.savefig(aggr_fig_address)
                    plt.close()
                    
                    
                    #sc_interactions_pd['i']=sc_interactions_pd['i']
                    #sc_interactions_pd['j']=sc_interactions_pd['j']
                    # temp=aggr_interactions_pd.loc[(aggr_interactions_pd['i']>=start)&(aggr_interactions_pd['j']>=start)&(aggr_interactions_pd['i']<end)\
                    #                     &(aggr_interactions_pd['j']<end)]
                    sub_matrix=similarity[start:end,start:end]
                    sub_matrix[np.isnan(sub_matrix)] = 0
                    from plot.plot_units import show_mat_loops
                    show_mat_loops(sub_matrix,start,end,upper_limit,temp,count=0)
                    
                    
                    #########plot decision graph
                    temp_raw=sc_interactions_pd_raw.loc[(sc_interactions_pd_raw['i']>=start)&(sc_interactions_pd_raw['j']>=start)&(sc_interactions_pd_raw['i']<end)\
                                        &(sc_interactions_pd_raw['j']<end)]
                    temp_raw.loc[temp_raw['cluster']==-1,'cluster']=0
                    
                    from plot.plot_units import red1,blue1,purple1    
                    colors=[blue1,red1];markers=['.','o']
                    ax=plt.figure(figsize=(4,4)).add_subplot(111)
                    for c in temp_raw['cluster'].unique():
                        color=colors[c]
                        
                        #for sc in [0,1]:
                        marker=markers[c]
                        indexes= (temp_raw['cluster']==c)
                        xs=temp_raw.loc[indexes,'strength']
                        ys=temp_raw.loc[indexes,'hdmd']
                        if c==0:
                            ax.scatter(xs,ys,s=10,c=color,marker=marker)
                        else:
                            ax.scatter(xs,ys,s=20,c=color,marker=marker)
                    
                    ax.set_xlabel('')
                    ax.set_ylabel('')
                    
                    ax.grid(False)
                    plt.gcf()
                    ax.set_facecolor('w')
                    
                    # ax.set_xlim([0.01,1])
                    # ax.set_ylim([0.07,1])
                    
                    ax = plt.gca()
                    ax.spines['right'].set_color('none')
                    ax.spines['top'].set_color('none')
                    ax.xaxis.set_ticks_position('bottom')
                    ax.yaxis.set_ticks_position('left')
                    
                    
                    #########plot rank graph
                    temp_raw['N_rank']=temp_raw['rank']/np.nanmax(temp_raw['rank'])
                    temp_raw['N_strength*hdmd']=temp_raw['strength*hdmd']/np.nanmax(temp_raw['strength*hdmd'])
                    
                    
                    ax=plt.figure(figsize=(4,4)).add_subplot(111)
                    for c in temp_raw['cluster'].unique():
                        color=colors[c]
                        
                        marker=markers[c]
                        indexes= (temp_raw['cluster']==c) 
                        xs=temp_raw.loc[indexes,'N_rank']
                        ys=temp_raw.loc[indexes,'N_strength*hdmd']
                        if c==0:
                            ax.scatter(xs,ys,s=10,c=color,marker=marker)
                        else:
                            ax.scatter(xs,ys,s=20,c=color,marker=marker)
                                
                    ax.set_xlabel('')
                    ax.set_ylabel('')
                    
                    # ax.set_xlim([0,0.06])
                    # ax.set_ylim([0,0.06])
                    
                    ax.grid(False)
                    plt.gcf()
                    ax.set_facecolor('w')
                    ax = plt.gca()
                    ax.spines['right'].set_color('none')
                    ax.spines['top'].set_color('none')
                    ax.xaxis.set_ticks_position('bottom')
                    ax.yaxis.set_ticks_position('left')
                    
                    
                    

               
                
                #sc_interactions_pd_raw=sc_interactions_pd_raw[['i','j','strength']]
                # #sc_interactions_pd['value']=[sub_mat[i,j] for i,j in zip(sc_interactions_pd['i'].tolist(),sc_interactions_pd['j'].tolist())]
                #sc_interactions_pd_raw.to_csv(raw_address,sep = ",",index = False,header=True)  
                
        #     rows=sc_interactions_pd['i'].values
        #     cols=sc_interactions_pd['j'].values
        #     matrix = sparse.csr_matrix(([1 for i in range(len(sc_interactions_pd))], (rows,cols)),shape=(mat_size,mat_size)).toarray()
        #     loop_num_mat+=matrix
        # allmatrix_sp=sparse.csr_matrix(loop_num_mat) 
        # sparse.save_npz(os.path.join(outdir, ".".join(["loop_num",chrom,"csv"])),allmatrix_sp)
        #loop_num_mat = coo_matrix((_data, (_row, _col)), shape=(4, 4), dtype=np.int)
        
        # from visual_show import show_mat_loops
        
        # print(len(sc_interactions_pd))
        # start=500;end=1000
        # #sc_interactions_pd['i']=sc_interactions_pd['i']
        # #sc_interactions_pd['j']=sc_interactions_pd['j']
        # temp=sc_interactions_pd.loc[(sc_interactions_pd['i']>=start)&(sc_interactions_pd['j']>=start)&(sc_interactions_pd['i']<end)\
        #                     &(sc_interactions_pd['j']<end)]
        # sub_matrix=sub_mat[start:end,start:end]
        
        # sub_matrix[np.isnan(sub_matrix)] = 0
        # show_mat_loops(sub_matrix.T,start,end,upper_limit,temp)

# def show_loop_results(indir,outdir, contact_statis_pd ,chrom_lens,\
#                       rank = 0, n_proc = 1,window=20):
    
#     import random
#     import matplotlib.pyplot as plt
#     from  visual_show import show_mat_loops,show_mat_loops_multigraph
#     proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
#     completed_filenames = os.listdir(indir)
#     filter_cells=contact_statis_pd.loc[contact_statis_pd[1]>250000,0].tolist()
    
    
#     for chrom in proc_chroms:
#         #chrom=proc_chroms[1]
#         using_suffix=".".join([chrom,  "rl", "npz"])
#         filenames=[os.path.join(indir,name) for name in completed_filenames if name.split('.',1)[1]==using_suffix] 
#         setnames = [os.path.basename(fname)[:-(len(using_suffix)+1)] for fname in filenames]
        
#         fliter_filenames=[filenames[setnames.index(cell)] for k,cell in enumerate(filter_cells) if cell in setnames]
        
        
        
#         #fliter_filenames=[ filename for k,filename in enumerate(filenames) if setnames[k] in filter_cells]
#         #filter_cells=[ setname for k,setname in enumerate(setnames) if setnames[k] in filter_cells]
        
#         l=386
#         #l=800
#         cell=filter_cells[l]
#         address=os.path.join(outdir, ".".join([cell,"interactions",chrom, "csv"]))
#         if os.path.exists(address):
#             sc_interactions_pd=pd.read_csv(address,sep = ",",index_col = None,header=0) 
#         filename=os.path.join(indir, ".".join([cell,chrom, "rl.npz"]))
#         csr_mat = sparse.load_npz(filename)
#         sub_mat = csr_mat.toarray()
    
#         print(len(sc_interactions_pd))
#         start=200;end=600;upper_limit=2e6
#         #sc_interactions_pd['i']=sc_interactions_pd['i']
#         #sc_interactions_pd['j']=sc_interactions_pd['j']
#         temp=sc_interactions_pd.loc[(sc_interactions_pd['i']>=start)&(sc_interactions_pd['j']>=start)&(sc_interactions_pd['i']<end)\
#                             &(sc_interactions_pd['j']<end)]
#         sub_matrix=sub_mat[start:end,start:end]
#         sub_matrix[np.isnan(sub_matrix)] = 0
#         show_mat_loops(sub_matrix.T,start,end,upper_limit,temp)


def bulkloops(indir,outdir,chrom_lens, binsize, dist, neighborhood_limit_lower, \
            neighborhood_limit_upper,alpha=1, rank = 0, n_proc = 1, max_mem = 2, logger = None,Figure=False):
    
    logger.set_rank(rank)
    try:
        os.makedirs(outdir)
    except:
        pass
    #table for loop distance
    MAXR=100
    DIS=np.zeros((2*MAXR+1,2*MAXR+1))
    for i in range(MAXR+1):
        for j in range(MAXR+1):
            dis=np.sqrt(i*i+j*j)
            DIS[MAXR+i,MAXR+j]=dis;DIS[MAXR+i,MAXR-j]=dis;DIS[MAXR-i,MAXR+j]=dis;DIS[MAXR-i,MAXR-j]=dis
    
    upper_limit=neighborhood_limit_upper;
    lower_limit=neighborhood_limit_lower; 
    Max_bin_distance=dist
    
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)

    #chrom='chr2'
    mat_size=1200
    for chrom in proc_chroms:
    
        output_add=os.path.join(indir,".".join(["sc",chrom,"csv"]))
        candidated_loops=pd.read_csv(output_add,header=0)
        mat_size = int(np.ceil(chrom_lens[chrom]/binsize))
        cell_mat = sp.sparse.csr_matrix((candidated_loops['weight'],(candidated_loops['i'],candidated_loops['j'])), \
                                        shape=(mat_size,mat_size)) 
        sub_mat=cell_mat.toarray()
        
        
        g = nx.Graph()
        g.add_weighted_edges_from(sub_mat)
        
        model = Node2Vec(g, dimensions=32,walk_length=3,num_walks=50,workers=4,p=1,q=1)  # Use temp_folder for big graphs
        model = model.fit(window=3)  # Any keywords acceptable by gensim.Word2Vec can be passed, `dimensions` and `workers` are automatically passed (from the Node2Vec constructor)

        #model.wv.most_similar('200.0')
        embeddings = {}
        for word in g.nodes():
            embeddings[word] = model.wv[str(word)] 
        similarity=cosine_similarity(pd.DataFrame(embeddings).T)
        
        
        if Figure==True:
            #read candidated_loops
            start=4550;end=4750
            sub_matrix=sub_mat[start:end,start:end]
            from plot.plot_units import show_mat
            show_mat(sub_matrix+sub_matrix.T,count=1)
        

        if sub_mat.shape[0]<=mat_size:
            #sub_mat=sub_mat[:,600:1800][600:1800,:];sc_sub_mat=sc_sub_mat[:,600:1800][600:1800,:]
            sc_interactions_pd_a,sc_interactions_pd_raw_a=compute_significance_interaction_bulk(sub_mat,similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                upper_limit,\
                                                                lower_limit,alpha) 
        else:
            submatrices=get_submatric_bulk(sub_mat,similarity,upper_limit,Max_bin_distance,binsize,mat_size)
            for i,(sub_matric,sub_similarity,start_index,end_flag) in enumerate(submatrices):
                print(i)
                # if i==6:
                #     break
            
            
                sc_interactions_pd,sc_interactions_pd_raw=compute_significance_interaction_bulk(sub_matric,sub_similarity,DIS,MAXR,Max_bin_distance,binsize,\
                                                                    upper_limit,\
                                                                    lower_limit,alpha) 
                
                sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']>upper_limit) & (sc_interactions_pd['j']>upper_limit)]
                sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']>upper_limit) & (sc_interactions_pd_raw['j']>upper_limit)]   
                
                if end_flag==0:
                    sc_interactions_pd=sc_interactions_pd[(sc_interactions_pd['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd['j']<(mat_size-Max_bin_distance//binsize))]
                    sc_interactions_pd_raw=sc_interactions_pd_raw[(sc_interactions_pd_raw['i']<(mat_size-Max_bin_distance//binsize)) & (sc_interactions_pd_raw['j']<(mat_size-Max_bin_distance//binsize))]
                
                if start_index==0:
                    sc_interactions_pd_a=sc_interactions_pd
                    sc_interactions_pd_raw_a=sc_interactions_pd_raw
                else:
                    sc_interactions_pd['i']=sc_interactions_pd['i']+start_index
                    sc_interactions_pd['j']=sc_interactions_pd['j']+start_index
                    sc_interactions_pd_raw['i']=sc_interactions_pd_raw['i']+start_index
                    sc_interactions_pd_raw['j']=sc_interactions_pd_raw['j']+start_index
                    
                    sc_interactions_pd_a=pd.concat([sc_interactions_pd_a,sc_interactions_pd])
                    sc_interactions_pd_raw_a=pd.concat([sc_interactions_pd_raw_a,sc_interactions_pd_raw])
                
                
        sc_interactions_pd=sc_interactions_pd_a 
        sc_interactions_pd_raw=sc_interactions_pd_raw_a

            
            
    
    
            
            
            
            
            
            


            
            
            
            
        
            
    
    
    
    

        
       
       
       
    
    
    
        
        
        
        
                    
            
        
        
        
    

    



