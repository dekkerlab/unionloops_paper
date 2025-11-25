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

dots_file='/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/unionloops_5kb_new/results/Bond_et_al_2023_union_loop_list_5kb.tsv'
dots=pd.read_csv(dots_file, sep='\t')

start_point='K562_0000min'
end_point='K562_4320min'
dots=dots.loc[(dots.sample_name.str.contains(start_point)) | (dots.sample_name.str.contains(end_point))].reset_index(drop=True)
dots['status']='static'
dots.loc[(dots.sample_name.str.contains(start_point)) & ~(dots.sample_name.str.contains(end_point)), 'status']='lost'
dots.loc[~(dots.sample_name.str.contains(start_point)) & (dots.sample_name.str.contains(end_point)), 'status']='gained'

resolution=5000
coolers={'K562_0000min_R1': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_0000min_biorep1__hg38/results/coolers_library/HIC_K562_0000min_biorep1__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'K562_0000min_R2': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_0000min_biorep2__hg38/results/coolers_library/HIC_K562_0000min_biorep2__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'K562_0000min_R3': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_0000min_biorep3__hg38/results/coolers_library/HIC_K562_0000min_biorep3__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'K562_0000min_R4': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_0000min_biorep4__hg38/results/coolers_library/HIC_K562_0000min_biorep4__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'K562_4320min_R1': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_4320min_biorep1__hg38/results/coolers_library/HIC_K562_4320min_biorep1__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'K562_4320min_R2': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_4320min_biorep2__hg38/results/coolers_library/HIC_K562_4320min_biorep2__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'K562_4320min_R3': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_4320min_biorep3__hg38/results/coolers_library/HIC_K562_4320min_biorep3__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}',
         'K562_4320min_R4': '/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/coolers/HIC_K562_4320min_biorep4__hg38/results/coolers_library/HIC_K562_4320min_biorep4__hg38.hg38.mapq_30.1000.mcool'+f'::resolutions/{resolution}'
        }

count_df=annotate_counts(dots, coolers).drop(['sample_name', 'K562_0000min', 'K562_0360min', 'K562_4320min', 'MEGA', 'bin1_id', 'bin2_id'], axis=1)

count_df.to_csv('/home/jiangyuan.liu5-umw/projects/4DNII/public_data_verification_for_union_list_loop_caller/public_datasets/Bond_et_al_2023/In_situ_HiC/unionloops_5kb_new/counts_of_each_biorep_of_HiC_for_loops_from_0_and_72h_5kb_resolution_new_method.tsv', header=True, index=False, sep='\t')







