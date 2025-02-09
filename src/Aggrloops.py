# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 16:48:25 2022

@author: user
"""

import os
import pandas as pd
import numpy as np
import scipy as sp
from plot.plot_units import show_mat


def get_proc_chroms(chrom_lens, rank, n_proc):
    chrom_list = [(k, chrom_lens[k]) for k in list(chrom_lens.keys())]
    chrom_list.sort(key = lambda x: x[1])
    chrom_list.reverse()
    chrom_names = [i[0] for i in chrom_list]
    #chrom_names = list(chrom_lens.keys())
    #chrom_names.sort()
    
    indices = list(range(rank, len(chrom_names), n_proc))
    proc_chroms = [chrom_names[i] for i in indices]
    return proc_chroms

def aggreloops(indir,outdir ,chrom_lens,binsize,rank, n_proc,logger,Identify_loops=False,Figure=False):
    try:
        os.makedirs(outdir)
    except:
        pass
    proc_chroms = get_proc_chroms(chrom_lens, rank, n_proc)
    
    #proc_chroms = get_proc_chroms(chrom_lens, rank, 1)[0:6]
    completed_filenames = os.listdir(indir)
    # # aggr aggrloops
    # for chrom in proc_chroms:
    #     #chrom=proc_chroms[0]
    #     using_suffix=".".join(['aggr','interactions',chrom, "csv"])
    #     filenames=[os.path.join(indir,name) for name in completed_filenames if name.split('.',1)[1]==using_suffix] 
    #     #setnames = [os.path.basename(fname)[:-(len(using_suffix)+1)] for fname in filenames]
    #     logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
    #     #merge loops
    #     candidated_loops=pd.DataFrame()
    #     # all_loop_num=[]
    #     # old_loop_num=0
    #     for k,filename in enumerate(filenames):
    #         #filename=filenames[120]
    #         temp_loops=pd.read_csv(filename,index_col = None,header=0)
    #         temp_loops=temp_loops[['i','j']]
    #         print(k,len(temp_loops))
    #         candidated_loops=pd.concat([candidated_loops,temp_loops])
    #         #candidated_loops.drop_duplicates(inplace=True)
    #         # new_loop_num=len(candidated_loops)
    #         # all_loop_num.append(new_loop_num)
    #         # print("the ratio of new updated loops: ",(new_loop_num-old_loop_num)*100/len(temp_loops))
    #         # old_loop_num=new_loop_num
        
    #     candidated_loops=candidated_loops.sort_values(by=['i','j'])
    #     candidated_loops['weight']=1
    #     candidated_loops=candidated_loops.groupby(['i','j']).sum()
    #     candidated_loops=candidated_loops.reset_index()
    #     #candidated_loops=candidated_loops[['i','j','weight']]
        
    #     output_add=os.path.join(outdir,".".join(["aggr",chrom,"csv"]))
    #     candidated_loops.to_csv(output_add,sep = ",",index = False,header=True)
        
    #     # import seaborn as sns 
    #     # sns.scatterplot(range(len(all_loop_num)),np.log2(all_loop_num))
        
    #     #show candidated_loops
    #     #candidated_loops=candidated_loops[candidated_loops['weight']>5]
    #     mat_size = int(np.ceil(chrom_lens[chrom]/binsize))
    #     cell_mat = sp.sparse.csr_matrix((candidated_loops['weight'],(candidated_loops['i'],candidated_loops['j'])), \
    #                                     shape=(mat_size,mat_size)) 
    #     # cell_mat = sp.sparse.csr_matrix(([1]*len(candidated_loops),(candidated_loops['i'],candidated_loops['j'])), \
    #     #                                 shape=(mat_size,mat_size)) 
    #     cell_mat=cell_mat.toarray()
    #     from plot.plot_units import show_mat
    #     show_mat(cell_mat,count=1)
    ##### aggr sc loops
    for chrom in proc_chroms:
        #chrom='chr2'
        output_add=os.path.join(outdir,".".join(["sc",chrom,"csv"]))
        if os.path.exists(output_add)==False:
            using_suffix=".".join(['sc','interactions',chrom, "csv"])
            filenames=[os.path.join(indir,name) for name in completed_filenames if name.split('.',1)[1]==using_suffix] 
            #setnames = [os.path.basename(fname)[:-(len(using_suffix)+1)] for fname in filenames]
            logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
            #merge loops
            candidated_loops=pd.DataFrame()
            # all_loop_num=[]
            # old_loop_num=0
            for k,filename in enumerate(filenames):
                #filename=filenames[1700]
                temp_loops=pd.read_csv(filename,index_col = None,header=0)
                temp_loops=temp_loops[['i','j']]
                print(k,len(temp_loops))
                candidated_loops=pd.concat([candidated_loops,temp_loops])
                #candidated_loops.drop_duplicates(inplace=True)
                # new_loop_num=len(candidated_loops)
                # all_loop_num.append(new_loop_num)
                # print("the ratio of new updated loops: ",(new_loop_num-old_loop_num)*100/len(temp_loops))
                # old_loop_num=new_loop_num
            candidated_loops=candidated_loops.sort_values(by=['i','j'])
            candidated_loops['weight']=1
            candidated_loops=candidated_loops.groupby(['i','j']).sum()
            candidated_loops=candidated_loops.reset_index()
            #candidated_loops=candidated_loops[['i','j','weight']]
            candidated_loops.to_csv(output_add,sep = ",",index = False,header=True)
            
        if Figure==True:
            output_add=os.path.join(outdir,".".join(["sc",chrom,"csv"]))
            candidated_loops=pd.read_csv(output_add,header=0)
            mat_size = int(np.ceil(chrom_lens[chrom]/binsize))
            cell_mat = sp.sparse.csr_matrix((candidated_loops['weight'],(candidated_loops['i'],candidated_loops['j'])), \
                                            shape=(mat_size,mat_size)) 
            sub_mat=cell_mat.toarray()#read candidated_loops
            start=4000;end=4550
            sub_matrix=sub_mat[start:end,start:end]
            
            show_mat(sub_matrix+sub_matrix.T,count=1)
    
    
    #### aggr sc weight
    for chrom in proc_chroms:
        #chrom='chr2'
        output_add=os.path.join(outdir,".".join(["sc_raw",chrom,"csv"]))
        if os.path.exists(output_add):
            using_suffix=".".join(['sc_raw','interactions',chrom, "csv"])
            
            filenames=[os.path.join(indir,name) for name in completed_filenames if name.split('.',1)[1]==using_suffix] 
            #setnames = [os.path.basename(fname)[:-(len(using_suffix)+1)] for fname in filenames]
            logger.write(f'\tprocessor {rank}: computing for chromosome {chrom}', verbose_level = 1, allow_all_ranks = True)
            #merge loops
            candidated_loops=pd.DataFrame()
            # all_loop_num=[]
            # old_loop_num=0
            for k,filename in enumerate(filenames):
                #filename=filenames[1700]
                temp_loops=pd.read_csv(filename,index_col = None,header=0)
                temp_loops=temp_loops[['i','j','strength']]
                print(k,len(temp_loops))
                candidated_loops=pd.concat([candidated_loops,temp_loops])
                #candidated_loops.drop_duplicates(inplace=True)
                # new_loop_num=len(candidated_loops)
                # all_loop_num.append(new_loop_num)
                # print("the ratio of new updated loops: ",(new_loop_num-old_loop_num)*100/len(temp_loops))
                # old_loop_num=new_loop_num
            candidated_loops=candidated_loops.sort_values(by=['i','j'])
            #candidated_loops['weight']=1
            candidated_loops=candidated_loops.groupby(['i','j']).sum()
            candidated_loops=candidated_loops.reset_index()
            #candidated_loops=candidated_loops[['i','j','weight']]
            candidated_loops.to_csv(output_add,sep = ",",index = False,header=True)
            
        if Figure==True:
            
            candidated_loops=pd.read_csv(output_add,header=0)
            
            
            mat_size = int(np.ceil(chrom_lens[chrom]/binsize))
            cell_mat = sp.sparse.csr_matrix((candidated_loops['strength'],(candidated_loops['i'],candidated_loops['j'])), \
                                            shape=(mat_size,mat_size)) 
            sub_mat=cell_mat.toarray()#read candidated_loops
            start=4550;end=4750
            sub_matrix=sub_mat[start:end,start:end]
            show_mat(sub_matrix+sub_matrix.T,count=1)


    


    




        
        
    
        
           
            
         
            
        
        
        
        
            
            

    
    
    
    
                