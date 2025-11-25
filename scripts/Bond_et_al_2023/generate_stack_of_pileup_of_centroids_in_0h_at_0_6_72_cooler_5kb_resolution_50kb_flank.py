#!/usr/bin/env python3

import pandas as pd
import numpy as np
import cooltools
import cooler
import bioframe
import pickle


centroids=pd.read_csv('/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/unionloops_5kb_old/centroids/HIC_K562_0000min_biorep1_to_4__hg38.hg38.mapq_30.5kb.dots.no.HiCCUPS.filtering.tsv', sep='\t')[['chrom1','start1','end1','chrom2','start2','end2']]

hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

# conditions=['K562_0000min','K562_0360min','K562_4320min']
conditions=['K562_0000min','K562_0360min','K562_4320min','MEGA']

cooler_meta='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/cooler_paths_of_Bond_et_al_2023.tsv'
clrs=pd.read_csv(cooler_meta, sep='\t')


stack={}
for cond in conditions:
    clr=clrs.loc[clrs.name == cond, 'path'].reset_index(drop=True)[0]
    clr = cooler.Cooler(clr+'::resolutions/5000')
    exp = cooltools.expected_cis(
        clr,
        view_df=hg38_arms,
        nproc=4,)
    stack[cond]=cooltools.pileup(clr, centroids, view_df=hg38_arms, expected_df=exp, flank=50000)

with open('/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/unionloops_5kb_old/centroids/stack_of_pileup_of_centroids_in_0h_at_0_6_72_cooler_5kb_resolution_50kb_flank.pickle', 'wb') as f:
    pickle.dump(stack, f)
