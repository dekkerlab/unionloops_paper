#!/usr/bin/env python3

import pandas as pd
import sys
import cooler

input_cooler_paths=pd.read_csv(sys.argv[1], sep='\t')
resolutions=sys.argv[2].split(",")
resolutions=[int(resolution) for resolution in resolutions]
ignore_diags=int(sys.argv[3])
output=sys.argv[4]


def generate_count_summary(conditions, clr_filenames, resolution, ignore_diags=2):
    total=[]
    total_trans=[]
    total_cis=[]
    total_cis_without_diags=[]
    for cond in conditions:
        clr=cooler.Cooler(clr_filenames[cond]+f'::resolutions/{resolution}')
        all_pixels = clr.pixels()[:]
        total.append(all_pixels['count'].sum())
        cis=[]
        cis_without_diags=[]
        for chrom in clr.chroms()[:].name:
            bins = clr.bins().fetch(chrom)[:]
            pixels = clr.pixels().fetch(chrom)[:]
            bins = bins[['chrom','start','end']].reset_index()
            merged_pixels = pixels.merge(bins, left_on='bin1_id', right_on='index')
            merged_pixels = merged_pixels.merge(bins, left_on='bin2_id', right_on='index', suffixes=('1', '2'))
            merged_pixels_without_diags=merged_pixels.loc[merged_pixels.start2 - merged_pixels.start1 >= ignore_diags*resolution].copy()
            cis.append(merged_pixels['count'].sum())
            cis_without_diags.append(merged_pixels_without_diags['count'].sum())
        total_cis.append(sum(cis))
        total_cis_without_diags.append(sum(cis_without_diags))
        total_trans.append(all_pixels['count'].sum() - sum(cis))

    count_summary = pd.DataFrame(
        {'condition': conditions,
         'total': total,
         'total_trans': total_trans,
         'total_cis': total_cis,
         f'total_cis_without_{ignore_diags}_diags': total_cis_without_diags
        })
    count_summary['total_trans_percent']=count_summary.total_trans/count_summary.total
    count_summary['total_cis_percent']=count_summary.total_cis/count_summary.total
    count_summary[f'total_cis_without_{ignore_diags}_diags_percent']=count_summary[f'total_cis_without_{ignore_diags}_diags']/count_summary['total']
    count_summary['resolution']=resolution

    
    return count_summary
    
output_list=[]
for resolution in resolutions:
    df=generate_count_summary(input_cooler_paths.name.tolist(), 
                              dict(zip(input_cooler_paths.name, input_cooler_paths.path)), 
                              resolution,
                              ignore_diags)
    output_list.append(df)

output_df=pd.concat(output_list, axis=0, ignore_index=True)
output_df.to_csv(output, sep='\t', index=False, header=True)


   
    
    
