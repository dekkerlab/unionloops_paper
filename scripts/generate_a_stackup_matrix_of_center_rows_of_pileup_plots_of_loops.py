#!/usr/bin/env python3

import pandas as pd
import numpy as np
import cooltools
import cooler
import bioframe
import pickle
import os.path
import sys


def generate_a_stackup_matrix_of_center_rows_of_pileup_matrices(dots_file, cooler_file, resolution, flank):

    dots=pd.read_csv(dots_file, sep='\t')[['chrom1','start1','end1','chrom2','start2','end2']]
    
    hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
    hg38_cens = bioframe.fetch_centromeres('hg38')
    hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
    
    clr = cooler.Cooler(cooler_file+'::resolutions/'+resolution)
    exp = cooltools.expected_cis(
        clr,
        view_df=hg38_arms,
        nproc=4,)
    pileup_mtx=cooltools.pileup(clr, dots, view_df=hg38_arms, expected_df=exp, flank=flank)
    # only select the center row of the pileup plot of each loop
    stackup_mtx=pileup_mtx[pileup_mtx.shape[0] // 2, :, :]
    
    return stackup_mtx.T

dots_file=sys.argv[1]
cooler_file=sys.argv[2]
resolution=sys.argv[3]
flank=int(sys.argv[4])
output=sys.argv[5]


stackup_mtx=generate_a_stackup_matrix_of_center_rows_of_pileup_matrices(dots_file=dots_file,
                                                                        cooler_file=cooler_file,
                                                                        resolution=resolution,
                                                                        flank=flank)

with open(output, 'wb') as f:
    pickle.dump(stackup_mtx, f)



