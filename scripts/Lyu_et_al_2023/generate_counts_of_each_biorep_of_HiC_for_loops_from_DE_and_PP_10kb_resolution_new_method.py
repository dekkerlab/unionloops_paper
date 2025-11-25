#!/usr/bin/env python3

import cooler
import numpy as np
import pandas as pd

# Annotate with bin1_id, bin2_id, and count for all pixels 
# bin_id are same for all conditions with same resolution
def annotate_counts(enriched_pixels, clr_dict):
    # bin1_id, bin2_id
    df=enriched_pixels.copy()
    clr_for_bin_id=cooler.Cooler(list(clr_dict.values())[0])
    df['bin1_id']=df.apply(lambda x: clr_for_bin_id.offset((x.chrom1, x.start1, x.end1)), axis=1)
    df['bin2_id']=df.apply(lambda x: clr_for_bin_id.offset((x.chrom2, x.start2, x.end2)), axis=1)
    
    for condition, clr_file in clr_dict.items():
        clr=cooler.Cooler(clr_file)
        df[condition+'_count']=df.apply(lambda x: clr.matrix(balance=False)[x.bin1_id, x.bin2_id].item(), axis=1)
    return df

dots_file='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/unionloops_10kb_new/results/Lyu_et_al_2023_union_loop_list_10kb.tsv'
dots=pd.read_csv(dots_file, sep='\t')

start_point='DE_R1R2'
end_point='PP_R1R2'
dots=dots.loc[(dots.sample_name.str.contains(start_point)) | (dots.sample_name.str.contains(end_point))].reset_index(drop=True)
dots['status']='static'
dots.loc[(dots.sample_name.str.contains(start_point)) & ~(dots.sample_name.str.contains(end_point)), 'status']='lost'
dots.loc[~(dots.sample_name.str.contains(start_point)) & (dots.sample_name.str.contains(end_point)), 'status']='gained'

resolution=10000
coolers={'DE_R1': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/coolers/results/coolers_library/DE_R1__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'DE_R2': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/coolers/results/coolers_library/DE_R2__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'PP_R1': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/coolers/results/coolers_library/PP_R1__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'PP_R2': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/coolers/results/coolers_library/PP_R2__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}'
        }

count_df=annotate_counts(dots, coolers).drop(['sample_name', 'ESC_R1R2', 'DE_R1R2', 'PGT_R1R2', 'PP_R1R2', 'BETA_R1R2', 'MEGA', 'bin1_id', 'bin2_id'], axis=1)

count_df.to_csv('/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Lyu_et_al_2023/In_situ_HiC/unionloops_10kb_new/counts_of_each_biorep_of_HiC_for_loops_from_DE_and_PP_10kb_resolution_new_method.tsv', header=True, index=False, sep='\t')







