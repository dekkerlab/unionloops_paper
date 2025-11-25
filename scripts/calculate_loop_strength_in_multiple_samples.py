#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import cooler
import bioframe
import cooltools

input_cooler_paths=pd.read_csv(sys.argv[1], sep='\t')
dots=pd.read_csv(sys.argv[2], sep='\t').loc[:,['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2']]
libraries_names=sys.argv[3].split(",")
resolution=sys.argv[4]
flank=int(sys.argv[5])
output=sys.argv[6]


hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

def quantify_loops(mtx):
    sq_size = mtx.shape[0]
    midpoint = int(np.floor(sq_size/2))
    mid_9pixels_mean = np.nanmean(mtx[midpoint-1:midpoint+2,midpoint-1:midpoint+2])
    neighboring_size = int(np.ceil(0.3*sq_size) // 2 * 2 + 1) # get the closest odd number of 30% sq_size
    upper_left_mean = np.nanmean(mtx[:neighboring_size,:neighboring_size])
    upper_right_mean = np.nanmean(mtx[:neighboring_size, (sq_size-neighboring_size):])
    lower_right_mean = np.nanmean(mtx[(sq_size-neighboring_size):, (sq_size-neighboring_size):])
    return mid_9pixels_mean/np.nanmean(np.array([upper_left_mean,upper_right_mean,lower_right_mean]))

def quantify_individual_loops_in_stack(stack):
    loop_strength=np.array([])
    for loop_idx in np.arange(stack.shape[-1]):
        loop_strength=np.append(loop_strength, quantify_loops(stack[:,:,loop_idx]))
    return loop_strength

loop_strengths=pd.DataFrame()
for cond in libraries_names:
    clr = cooler.Cooler(input_cooler_paths.loc[input_cooler_paths.name == cond,'path'].reset_index(drop=True)[0]+'::resolutions/'+resolution)
    expected = cooltools.expected_cis(
        clr,
        view_df=hg38_arms,
        nproc=4,)
    pileup_mtx=cooltools.pileup(clr, dots, view_df=hg38_arms, expected_df=expected, flank=flank)
    loop_strengths[cond] = quantify_individual_loops_in_stack(pileup_mtx)

loop_strengths_formatted=pd.concat([dots, loop_strengths], axis=1)
# loop_strengths_formatted.replace([np.inf, -np.inf], np.nan, inplace=True) # drop +/-inf
# loop_strengths_formatted.dropna(inplace=True)

loop_strengths_formatted.to_csv(output, sep='\t', index=False, header=True)

