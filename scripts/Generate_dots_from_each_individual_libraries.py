#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import cooler
import bioframe
import cooltools
import pickle
from cooltools.api.dotfinder import *

input_cooler_paths=pd.read_csv(sys.argv[1], sep='\t')
output_folder=sys.argv[2]

# Check whether the specified path exists or not
isExist = os.path.exists(output_folder)
if not isExist:
   # Create a new directory because it does not exist
   os.makedirs(output_folder)

hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

# Generate enriched pixels at 5kb resolution for each matrix
dots_dict={}
for i in range(len(input_cooler_paths)):
    clr=cooler.Cooler(input_cooler_paths.path[i]+'::resolutions/5000')
    expected = cooltools.expected_cis(
        clr,
        view_df=hg38_arms,
        nproc=4,)
    
    dots_dict[input_cooler_paths.name[i]] = cooltools.dots(
        clr,
        expected=expected,
        view_df=hg38_arms,
        lambda_bin_fdr=0.1,
        max_loci_separation=10_000_000,
        tile_size=5_000_000,
        max_nans_tolerated=4,
        clustering_radius=20_000,
        nproc=4,)


dots_filenames=[]
for i in dots_dict.keys():
    filename=os.path.join(output_folder,os.path.basename(input_cooler_paths.loc[input_cooler_paths.name == i,'path'].copy().reset_index(drop=True)[0])).replace('.1000.mcool', '.5kb.dots.tsv')
    dots_dict[i].to_csv(filename, sep='\t', index=False, header=True)
    dots_filenames.append(filename)

meta_dots_files = {
    'name' : list(dots_dict.keys()),
    'path' : dots_filenames
}


df = pd.DataFrame(meta_dots_files)                          
df.to_csv(os.path.join(output_folder,'dots_meta_5kb.tsv'), sep='\t', index=False, header=True)

