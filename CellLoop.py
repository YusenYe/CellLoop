# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 14:32:26 2022

@author: user
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
#########################Initial parameters##############################################
###please write the directory where CellLoop package is located#####
main_dir='/home/dell/Desktop/CellLoop_test/CellLoop'
os.chdir(main_dir)
######select the data type and bin resolution you want to analyze######
DATASET='Dip-C'
#DATASET='HiRES'
#DATASET='GAGE-seq'
BINSIZE=100e3
if BINSIZE<=10e3:
    MAXDIST=5000000
    knn_cell_num=50
else:
    MAXDIST=10000000
    knn_cell_num=50
#default parameter to cut-off for removing short-range reads
LOW_CUTOFF=1e3   
if DATASET=='Dip-C':
    GENOME='mm10'
    #dataset dir
    INDIR="/mnt/nas/Datasets/2020_Dip-C/GSE162511"##
    FILE_SUFFIX=".contacts.pairs.txt.gz"
    MINCONCNUM=250000
    CHR_COLUMNS=[1,3]
    POS_COLUMNS=[2,4]
    feature='higashi_embed'
    OUTDIR=INDIR+'_'+str(int(BINSIZE//1000))+'kb'+'_knn='+str(knn_cell_num)
    
elif DATASET=='HiRES':
    GENOME='mm10'
    INDIR="/mnt/nas/Datasets/2023_Science/GSE223917"
    FILE_SUFFIX=".pairs.gz"
    MINCONCNUM=250000
    CHR_COLUMNS=[1,3]
    POS_COLUMNS=[2,4]
    feature='RNA_embed'
    OUTDIR=INDIR+'_'+str(int(BINSIZE//1000))+'kb'+'_knn='+str(knn_cell_num)
        
elif DATASET=='GAGE-seq':
    GENOME='mm10'
    INDIR="/mnt/data/Datasets/GAGE-seq/SingleCells"
    FILE_SUFFIX=".pairs"
    MINCONCNUM=200000
    CHR_COLUMNS=[0,2]
    POS_COLUMNS=[1,3]
    feature='RNA_embed'
    OUTDIR=INDIR+'_'+str(int(BINSIZE//1000))+'kb'+'_knn='+str(knn_cell_num)


if GENOME=='hg19':
    CHR_LEN=main_dir+'/ext/hg19.chrom.sizes'
    FILTER_FILE=main_dir+'/ext/hg19_filter_regions.txt'
elif GENOME=='mm10':
    CHR_LEN=main_dir+'/ext/mm10.chrom.sizes'
    FILTER_FILE=main_dir+'/ext/mm10_filter_regions.txt'
elif GENOME=='hg38':
    CHR_LEN='main_dir+/ext/hg38.chrom.sizes.txt'
    FILTER_FILE=main_dir+'/ext/hg38_filter_regions.txt'
###################################################################################
import argparse
import src.logger
#import time
from src.bin_reads import bin_sets
import pandas as pd 
#from src.RepresentLearning_95 import get_rl_for_all_95
#from src.sc_compartments import read_sc_compartment
#from src.sc_tadboundaries import sctad_boundary
from src.sc_interactions import scloops
#from visual_show import visual_show
#from src.commonorspecific_boundaries import comorspebound
#from src.commonorspecific_loops import comorspeloops
from src.sample_aggr_cells import get_aggr_cells
#from src.Computing_contactprob import con_prob_dis
from src.Aggrloops import aggreloops
from src.clustercells_loops import clustercells
#from src.aggremaps_visualanalyzing import aggremaps,aggreloops_cluster,aggr_bulkmaps,aggr_bulkloops
#from src.statistic_analysis_loops import statistic_analysis_loops
import multiprocessing
import scanpy as sc

#import pandas as pd
#########################################################################################
def main():
    #__spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    parser = create_parser()
    args = parser.parse_args()
    chrom_dict = parse_chrom_lengths(args.chrom, args.chr_lens, args.genome, args.max_chrom_number)
    parallel_mode, rank, n_proc, parallel_properties = determine_parallelization_options(args.parallel, args.threaded, args.num_proc)
    
    if rank == 0 and not os.path.exists(args.outdir):
        os.makedirs(args.outdir)
    if parallel_mode == "parallel":
        parallel_properties['comm'].Barrier()
    threaded = True if parallel_mode == "threaded" else False
    logger = src.logger.Logger(f'{args.outdir}/scloops.log', rank = rank, verbose_threshold = args.verbose, threaded = threaded)
    
    ###############################################################################################
    #step 1; binning
    bin_dir = os.path.join(args.outdir, "binned")
    if 'bin' in args.steps:
        logger.write('starting the binning step')
        logger.flush()
        if parallel_mode == 'nonparallel':
            #indir=args.indir; suffix=args.suffix; binsize = args.binsize; outdir = bin_dir; 
            #chr_columns = args.chr_columns; pos_columns = args.pos_columns; 
            #low_cutoff = args.low_cutoff; mincontactnum=args.mincontactnum;dataset=args.dataset;n_proc = n_proc;rank = rank;logger = logger
            bin_sets(args.indir, args.suffix, binsize = args.binsize, outdir = bin_dir, \
                      chr_columns = args.chr_columns, pos_columns = args.pos_columns, \
                      low_cutoff = args.low_cutoff, mincontactnum=args.mincontactnum,dataset=args.dataset,n_proc = n_proc, rank = rank, logger = logger)
        elif parallel_mode == 'threaded':
            params = [(args.indir, args.suffix, args.binsize, bin_dir, args.chr_columns, args.pos_columns,\
                      args.low_cutoff,args.mincontactnum,args.dataset,n_proc,i,logger) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(bin_sets, params)
        logger.write("binning completed")
        logger.flush()
        
    ################### you can get using filelist ##########################
    filenames=os.listdir(bin_dir)
    filenames.sort()
    #for filelist
    f=open(args.indir+'/filelist.txt','w')
    for line in filenames:
        print(line)
        f.write(line.rstrip('.bedpe')+'\n')
    f.close()
    
    
    ####################################################################################################
    #step 1.1 #choose feature compute knn graph of single cells
    ####################################################################################################
    #indir=args.indir 
    if feature=='higashi_embed':
        feature_dir=main_dir+'/Dataset/'+DATASET+'/CellFeatures/'+feature
        filenames=os.listdir(bin_dir)
        filenames.sort()
        #for filelist
        f=open(feature_dir+'/filelist.txt','w')
        for line in filenames:
            print(line)
            f.write(line.rstrip('.bedpe')+'\n')
        f.close()
    
    elif feature=='RNA_embed':
        feature_dir=main_dir+'/Dataset/'+DATASET+'/CellFeatures/'+feature
        adata_rna=sc.read(feature_dir+'/adata_rna.h5ad')
        filelist_dir=args.indir+'/filelist.txt'
        filenames=pd.read_csv(filelist_dir, index_col=None,header=None,sep='\t')
        # filenames=os.listdir(bin_dir)
        # filenames.sort()
        adata_rna.obs_names=adata_rna.obs['cell']
        #adata_rna=adata_rna[[cellname in filenames for cellname in adata_rna.obs_names],:]
        filenames=filenames[[cellname in adata_rna.obs_names for cellname in filenames[0].values]]
        adata_rna=adata_rna[filenames[0].values.tolist(),:]
        # sc.pp.neighbors(adata_rna, n_neighbors=15, n_pcs=20)
        # sc.tl.umap(adata_rna)
        # sc.tl.leiden(adata_rna)
        # sc.pl.umap(adata_rna,color=['n_genes','leiden'], palette='tab20',size=10,wspace=1,hspace=0.1)
        adata_rna.write(feature_dir+'/adata_rna.h5ad')
        #for filelist
        f=open(feature_dir+'/filelist.txt','w')
        for line in filenames[0]:
            print(line)
            f.write(line.rstrip('.bedpe')+'\n')
        f.close()
    
    
    #############################################################################################
    #### step 2. Generating the enhanced contact map of the single cells#########################
    # sampling and getting aggregated cells
    #choose feature compute knn graph of single cells
    SampleAggrCell_dir=os.path.join(args.outdir,'sampleaggrcell_dir')
    #k=20
    if 'sampleaggrcell' in args.steps:
        logger.write('starting the sampling aggragated cells step')
        logger.flush()
        if parallel_mode == 'nonparallel':
            #bin_dir=bin_dir;outdir=SampleAggrCell_dir;feature_dir=feature_dir;feature=feature;chrom_lens=chrom_dict;
            #n_proc=n_proc;rank=rank;binsize=args.binsize;k=knn_cell_num;
            get_aggr_cells(bin_dir=bin_dir,outdir=SampleAggrCell_dir,feature_dir=feature_dir,feature=feature,chrom_lens=chrom_dict,\
                           n_proc=n_proc, rank=rank,binsize=args.binsize,k=knn_cell_num)
                # bin_sets(args.indir, args.suffix, binsize = args.binsize, outdir = bin_dir, \
            #           chr_columns = args.chr_columns, pos_columns = args.pos_columns, \
            #           low_cutoff = args.low_cutoff, n_proc = n_proc, rank = rank, logger = logger)
        elif parallel_mode == 'threaded':
            params = [(bin_dir,SampleAggrCell_dir,feature_dir,feature,chrom_dict,n_proc,i,\
                      args.binsize,knn_cell_num) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(get_aggr_cells, params)
        logger.write("binning completed")
        logger.flush()
        
    
    ### Step 3 Calling chromatin loops of single cells################################################################
    scloops_dir=os.path.join(args.outdir,'scloops_new')
    # temp_keys={'chr2','chr13'}
    # chrom_dict={key:chrom_dict[key] for key in chrom_dict.keys() & temp_keys }
    alpha=1
    scloops_outdir=os.path.join(scloops_dir+'_'+feature,'alpha='+str(alpha))
    # scloops_dir+'_'+feature+'_alpha='+str(alpha)
    if 'scloops' in args.steps:
        indir=SampleAggrCell_dir+'_'+feature
        logger.write('starting the chromatin loops calling')
        logger.flush()
        if parallel_mode=="nonparallel":
            #chrom_lens = chrom_dict;outdir = scloops_outdir;indir = indir
            #binsize = args.binsize;dist = args.dist
            #neighborhood_limit_lower = args.local_lower_limit; neighborhood_limit_upper = args.local_upper_limit
            #rank = rank; n_proc = n_proc; max_mem = args.max_memory; logger = logger
            scloops(indir = indir, outdir = scloops_outdir, chrom_lens = chrom_dict, \
                            binsize = args.binsize, dist = args.dist, \
                            neighborhood_limit_lower = args.local_lower_limit, \
                            neighborhood_limit_upper = args.local_upper_limit,alpha=alpha,\
                            rank = rank, n_proc = n_proc, max_mem = args.max_memory, logger = logger)          
        elif parallel_mode == 'threaded':
            params = [(indir, scloops_outdir, chrom_dict, \
                            args.binsize,args.dist,args.local_lower_limit,args.local_upper_limit,alpha,\
                            i,n_proc, args.max_memory,logger) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(scloops, params)
        logger.write("calling completed")
        logger.flush()

    #############################################################################################
    ##### Step 4. Generating chromatin loop frequency map(LFmap)##############################################################
    aggreloops_dir=os.path.join(args.outdir,'aggreloops','alpha='+str(alpha))
    if 'aggreloops' in args.steps:
        logger.write('starting aggregated loops')
        logger.flush()
        if parallel_mode=="nonparallel":
            #indir = scloops_outdir; outdir = aggreloops_dir;chrom_lens = chrom_dict;
            #binsize = args.binsize;dist = args.dist
            #rank = rank; n_proc = n_proc; max_mem = args.max_memory; logger = logger
            aggreloops(indir = scloops_outdir,outdir = aggreloops_dir,\
                       chrom_lens = chrom_dict, binsize = args.binsize,
                       rank = rank, n_proc = n_proc, logger = logger)          
        elif parallel_mode == 'threaded':
            params = [(scloops_outdir, aggreloops_dir, chrom_dict, \
                            args.binsize,\
                            i,n_proc,logger) for i in range(n_proc)]
            with multiprocessing.Pool(n_proc) as pool:
                pool.starmap(aggreloops, params)
        logger.write("aggregating completed")
        logger.flush()
            
    
                  
    # ############a##################################################################################
    # ##### step 5. try cluster cells using loop features ##########################################
    # clustercells_dir=os.path.join(args.outdir,'clustercells','alpha='+str(alpha))
    # if 'nclustercells' in args.steps:
    #     if parallel_mode=="nonparallel":
    #         #indir = aggreloops_dir; scloops_dir=scloops_outdir;outdir =  clustercells_dir;chrom_lens = chrom_dict;
    #         #binsize = args.binsize;dist = args.dist
    #         #rank = rank; n_proc = n_proc; logger = logger
    #         clustercells(indir = aggreloops_dir,scloops_dir=scloops_outdir,outdir = clustercells_dir, chrom_lens = chrom_dict, binsize = args.binsize,
    #                       rank = rank, n_proc = n_proc, logger = logger)          
        



def parse_chrom_lengths(chrom, chrom_lens_filename, genome, max_chrom_number):
    if not max_chrom_number or max_chrom_number == -1:
        if not chrom or chrom == "None":
            chrom_count = 22 if genome.startswith('hg') else 19 if genome.startswith("mm") else None
            if not chrom_count:
                raise("Genome name is not recognized. Use --max-chrom-number")
            chrom = ['chr' + str(i) for i in range(1, chrom_count + 1)]
        else:
            chrom = [c.strip() for c in chrom.split()]
    else:
        chrom = ['chr' + str(i) for i in range(1, max_chrom_number + 1)]
    with open(chrom_lens_filename) as infile:
        lines = infile.readlines()
    chrom_lens = {line.split()[0]: int(line.split()[1]) for line in lines if line.split()[0] in chrom}
    return chrom_lens

def determine_parallelization_options(parallel, threaded, n_proc):
    if parallel and threaded:
        raise "Only one of 'parallel' or 'threaded' flags can be set. \
                If using a job scheduling system on a cluster with multiple machines, use parallel.\
                If using a single machine with multiple CPUs, use threaded"
    else:
        if parallel:
            #print('it;s parallel')
            from mpi4py import MPI
            comm = MPI.COMM_WORLD
            n_proc = comm.Get_size()
            rank = comm.Get_rank()
            #print('n proc is', n_proc)
            #print('rank is', rank)
            mode = 'parallel'
            properties = {'comm': comm}
        elif threaded:
            import multiprocessing
            if n_proc < 1:
                raise Exception('if threaded flag is set, n should be a positive integer')
            n_proc = n_proc
            mode = 'threaded'
            rank = 0
            properties = {}
        else:
            mode = 'nonparallel'
            n_proc = 1
            rank = 0
            properties = {}
    return mode, rank, n_proc, properties

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--indir', action = 'store', required =False, \
                        help = 'input directory',default=INDIR)
    parser.add_argument('-s', '--suffix', required = False, \
                        help = 'suffix of the input files', default=FILE_SUFFIX)
    parser.add_argument('-o', '--outdir', action = 'store', \
                        required = False, help = 'output directory', default=OUTDIR)
    parser.add_argument('-c', '--chr-columns', action = 'store', nargs = 2, \
                        type = int, help = 'two integer column numbers for chromosomes', required = False, default = CHR_COLUMNS)
    parser.add_argument('-p', '--pos-columns', action = 'store', nargs = 2, \
                        type = int, help = 'two integer column numbers for read positions', required = False, default = POS_COLUMNS)
    parser.add_argument('-l', '--chr-lens', action = 'store', \
                        help = 'path to the chromosome lengths file', required = False,\
                           default=CHR_LEN)
    parser.add_argument('-g', '--genome', action = 'store', help = 'genome name; hgxx or mmxx', \
                        required = False, default = "hg19")
    
    parser.add_argument('--chrom', action = 'store', help = 'chromosome to process', \
                        required = False, default = None)
    parser.add_argument('--dist', type = int, help = 'distance from diagonal to consider', \
                        default =MAXDIST, required = False)
    
    parser.add_argument('--mincontactnum', type = int, help = 'minimum contact number for single cells', \
                        default =MINCONCNUM, required = False)
    parser.add_argument('--dataset', action = 'store', help = 'choose dataset', \
                        default = DATASET, required = False)
        
    parser.add_argument('--binsize', type = int, help = 'bin size used for binning the reads', \
                        required = False, default = BINSIZE)
    parser.add_argument('--low-cutoff', type = int, help = 'cut-off for removing short-range reads', \
                        default = LOW_CUTOFF, required = False)
    
    parser.add_argument('--parallel', action = 'store_true', default = False, \
                        help = 'if set, will attempt to run in parallel mode', required = False)
    parser.add_argument('--threaded', action = 'store_true', default =True, \
                        help = 'if set, will attempt to use multiprocessing on single machine', required = False)
    parser.add_argument('-n', '--num-proc', help = 'number of processes used in threaded mode',
                        required = False, default = 10, type = int)
    
    parser.add_argument('--local-lower-limit', default = 2, type = int, required = False, \
                        help = 'number of bins around center (in each direction) to exlude from neighborhood')
    parser.add_argument('--local-upper-limit', default = 5, type = int, required = False, \
                        help = 'number of bins around center (in each direction) forming the neighborhood')
   
        
    parser.add_argument('--filter-file', default=FILTER_FILE, required = False, \
                        help = "bed file of regions to be filtered. Regions should be binned")
    parser.add_argument('--max-memory', default = 2, type = float, required = False, \
                        help = 'memory available in GB, that will be used in constructing dense matrices')
    
    parser.add_argument('--verbose', type = int, required = False, default = 0,
                        help = 'integer between 0 and 3 (inclusive), 0 for the least amount to print to log file') 
    
    parser.add_argument('--steps', nargs = "*", default = ['bin','sampleaggrcell','scloops','aggreloops'] , \
                        required = False, help = 'steps to run. Default is all steps.')
    parser.add_argument('--max-chrom-number', action = "store", required = False, type = int, \
                        help = "biggest chromosome number to consider in genome, for example 22 for hg", default = -1)
    parser.add_argument('--prefix', action = 'store', required = False, default = None,
                        help = "a prefix that will be added to the name of the final output files")
    return parser

if __name__ == "__main__":
    main()
    #passa
    








