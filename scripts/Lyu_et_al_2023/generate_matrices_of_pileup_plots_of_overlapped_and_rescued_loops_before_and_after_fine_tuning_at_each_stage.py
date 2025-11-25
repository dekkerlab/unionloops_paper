#!/usr/bin/env python3

import pandas as pd
import numpy as np
import cooltools
import cooler
import bioframe
import pickle

cluster_file_new='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/unionloops_5kb_new/results/clusters/clusters_of_enriched_pixels.resolution.5kb.tsv'
cluster_new=pd.read_csv(cluster_file_new, sep='\t')[['chrom1','start1','end1','chrom2','start2','end2','condition','region','c_label','c_size']]

dots_dir='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/unionloops_5kb_old/'
dots_old=pd.read_csv(dots_dir+'union_dot_list_of_Lyu_et_al_2023_5kb_old_method.txt', sep='\t')
# remove MEGA only loops
dots_old=dots_old.loc[dots_old.sample_name != 'MEGA']

dots_file='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/unionloops_5kb_new/results/Lyu_et_al_2023_union_loop_list_5kb.tsv'
dots_new=pd.read_csv(dots_file, sep='\t')[['chrom1','start1','end1','chrom2','start2','end2','sample_name']]
# remove MEGA only loops
dots_new=dots_new.loc[dots_new.sample_name != 'MEGA']

conditions=['ESC_R1R2','DE_R1R2','PGT_R1R2','PP_R1R2','BETA_R1R2']

dots_file_old_no_filtering='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/unionloops_5kb_old/centroids/union_dot_list_of_Lyu_et_al_2023_5kb_old_method_no_HiCCUPS_filtering.txt'

dots_old_in_cluster_new=pd.merge(dots_old, cluster_new, on=['chrom1','start1','end1','chrom2','start2','end2'])
dots_new_in_cluster_new=pd.merge(dots_new, cluster_new, on=['chrom1','start1','end1','chrom2','start2','end2'])


overlapped_loops={}
for cond in conditions:
    overlapped_loops[cond]={}
    # left: old; right: new
    df=pd.merge(dots_old_in_cluster_new.loc[dots_old_in_cluster_new.sample_name.str.contains(cond)],
                dots_new_in_cluster_new.loc[dots_new_in_cluster_new.sample_name.str.contains(cond)],
                on=['region','c_label','c_size'],
                how="outer",
                indicator=True)
    
    # filtered loops (old)
    dots_old_no_filtering=pd.read_csv(dots_file_old_no_filtering, sep='\t')
    dots_old_no_filtering=dots_old_no_filtering.loc[dots_old_no_filtering.sample_name.str.contains(cond), ['chrom1','start1','end1','chrom2','start2','end2']]
    dots_old_no_filtering_in_cluster_new=pd.merge(dots_old_no_filtering, cluster_new, on=['chrom1','start1','end1','chrom2','start2','end2'])

    columns_to_convert=['start1','end1','start2','end2']
    
    df_old_only=df.loc[df._merge == 'left_only', ['chrom1_x', 'start1_x', 'end1_x', 'chrom2_x', 'start2_x', 'end2_x', 'region', 'c_label', 'c_size']].set_axis(['chrom1','start1','end1','chrom2','start2','end2','region','c_label','c_size'], axis=1).drop_duplicates().reset_index(drop=True)
    df_old_only[columns_to_convert] = df_old_only[columns_to_convert].astype(int)
    overlapped_loops[cond]['old_only']=df_old_only
    
    df_new_only=df.loc[df._merge == 'right_only', ['chrom1_y', 'start1_y', 'end1_y', 'chrom2_y', 'start2_y', 'end2_y', 'region', 'c_label', 'c_size']].set_axis(['chrom1','start1','end1','chrom2','start2','end2', 'region','c_label','c_size'], axis=1).drop_duplicates().reset_index(drop=True)
    df_new_only[columns_to_convert] = df_new_only[columns_to_convert].astype(int)
    overlapped_loops[cond]['new_only_in_new']=df_new_only
    
    # find filtered loops (old) in new only clusters
    df_dots_old_no_filtering_in_new_only_clusters=pd.merge(dots_old_no_filtering_in_cluster_new, df_new_only, on=['region','c_label','c_size'])[['chrom1_x', 'start1_x', 'end1_x', 'chrom2_x', 'start2_x', 'end2_x', 'region', 'c_label', 'c_size']].set_axis(['chrom1','start1','end1','chrom2','start2','end2','region','c_label','c_size'], axis=1).drop_duplicates().reset_index(drop=True)
    df_dots_old_no_filtering_in_new_only_clusters[columns_to_convert] = df_dots_old_no_filtering_in_new_only_clusters[columns_to_convert].astype(int)
    overlapped_loops[cond]['new_only_in_old_no_filtering']=df_dots_old_no_filtering_in_new_only_clusters
    
    # multiple loops belong to the same cluster
    df_old_in_both=df.loc[df._merge == 'both', ['chrom1_x', 'start1_x', 'end1_x', 'chrom2_x', 'start2_x', 'end2_x', 'region', 'c_label', 'c_size']].set_axis(['chrom1','start1','end1','chrom2','start2','end2','region','c_label','c_size'], axis=1).drop_duplicates().reset_index(drop=True)
    df_old_in_both[columns_to_convert] = df_old_in_both[columns_to_convert].astype(int)
    overlapped_loops[cond]['old_in_both']=df_old_in_both
    
    df_new_in_both=df.loc[df._merge == 'both', ['chrom1_y', 'start1_y', 'end1_y', 'chrom2_y', 'start2_y', 'end2_y', 'region', 'c_label', 'c_size']].set_axis(['chrom1','start1','end1','chrom2','start2','end2','region','c_label','c_size'], axis=1).drop_duplicates().reset_index(drop=True)
    df_new_in_both[columns_to_convert] = df_new_in_both[columns_to_convert].astype(int)
    overlapped_loops[cond]['new_in_both']=df_new_in_both
        
methods=['old_only', 'old_in_both', 'new_in_both', 'new_only_in_old_no_filtering', 'new_only_in_new']
cooler_meta='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/cooler_paths_of_Lyu_et_al_2023.tsv'
clrs=pd.read_csv(cooler_meta, sep='\t')
resolution=5000
flank=50000
output='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/unionloops_5kb_new/pileup_matrices_of_overlapped_and_rescued_loops_before_and_after_fine_tuning_at_each_stage.pickle'


hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

matrices={}
for cond in conditions:
    for method in methods:
        clr = cooler.Cooler(clrs.path[clrs.name == cond].reset_index(drop=True)[0]+'::resolutions/'+str(resolution))
        exp = cooltools.expected_cis(
            clr,
            view_df=hg38_arms,
            nproc=4,)
        dots=overlapped_loops[cond][method][['chrom1','start1','end1','chrom2','start2','end2']]
        pileup_mtx=cooltools.pileup(clr, dots, view_df=hg38_arms, expected_df=exp, flank=flank)
        matrices[f'{cond}_{method}']=np.nanmean(pileup_mtx, axis=2)
        
with open(output, 'wb') as f:
    pickle.dump(matrices, f)





