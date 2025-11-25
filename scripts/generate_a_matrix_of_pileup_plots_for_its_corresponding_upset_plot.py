#!/usr/bin/env python3

import pandas as pd
import numpy as np
import cooltools
import cooler
import bioframe
import pickle
import os.path
import sys


def generate_pileup_matrices(dots_file, cooler_meta, libraries_names, dot_subsets_names, resolution, flank):

    dots=pd.read_csv(dots_file, sep='\t')
    clrs=pd.read_csv(cooler_meta, sep='\t')
    
    hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
    hg38_cens = bioframe.fetch_centromeres('hg38')
    hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
     
    # create an empty array to save pileup matrices
    mtx_arr = np.zeros((len(libraries_names),len(dot_subsets_names)), dtype=object)
    for i in range(len(libraries_names)): # rows in upset plot
        clr = cooler.Cooler(clrs.path[clrs.name == libraries_names[i]].reset_index(drop=True)[0]+'::resolutions/'+resolution)
        exp = cooltools.expected_cis(
            clr,
            view_df=hg38_arms,
            nproc=4,)
        for j in range(len(dot_subsets_names)): # columns in upset plot
            print(f'dot_subsets_name: {dot_subsets_names[j]}; library name: {libraries_names[i]}')
            subset_dots = dots.loc[dots['sample_name'] == dot_subsets_names[j],['chrom1','start1','end1','chrom2','start2','end2']]
            pileup_mtx=cooltools.pileup(clr, subset_dots, view_df=hg38_arms, expected_df=exp, flank=flank)
            mtx_arr[i,j]=np.nanmean(pileup_mtx, axis=2)
    
    return mtx_arr

dots_file=sys.argv[1]
cooler_meta=sys.argv[2]
libraries_names=sys.argv[3].split(",")
dot_subsets_names=sys.argv[4].split(",")
resolution=sys.argv[5]
flank=int(sys.argv[6])
output=sys.argv[7]


pileup_mtx=generate_pileup_matrices(dots_file=dots_file,
                                    cooler_meta=cooler_meta,
                                    libraries_names=libraries_names,
                                    dot_subsets_names=dot_subsets_names,
                                    resolution=resolution,
                                    flank=flank)

with open(output, 'wb') as f:
    pickle.dump(pileup_mtx, f)



