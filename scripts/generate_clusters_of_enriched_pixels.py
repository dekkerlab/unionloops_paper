#!/usr/bin/env python3

import pandas as pd
import numpy as np
import sys
import os
import bioframe
import cooler
import cooltools
from cooltools.api.dotfinder import *


# Modified clustering_step() function in dotfinder.py
# Only return clusters instead of centroids
def generate_clusters(
    scored_df,
    dots_clustering_radius,
    assigned_regions_name="region",
    obs_raw_name=observed_count_name,
):
    """
    Group together adjacent significant pixels into clusters after
    the lambda-binning multiple hypothesis testing by iterating over
    assigned regions and calling `clust_2D_pixels`.

    Parameters
    ----------
    scored_df : pandas.DataFrame
        DataFrame with enriched pixels that are ready to be
        clustered and are annotated with their genomic  coordinates.
    dots_clustering_radius : int
        Birch-clustering threshold.
    assigned_regions_name : str | None
        Name of the column in scored_df to use for grouping pixels
        before clustering. When None, full chromosome clustering is done.
    obs_raw_name : str
        name of the column with raw observed pixel counts
    Returns
    -------
    centroids : pandas.DataFrame
        Pixels from 'scored_df' annotated with clustering information.

    Notes
    -----
    'dots_clustering_radius' in Birch clustering algorithm corresponds to a
    double the clustering radius in the "greedy"-clustering used in HiCCUPS

    """
    # make sure provided pixels are annotated with genomic corrdinates and raw counts column is present:
    if not {"chrom1", "chrom2", "start1", "start2", obs_raw_name}.issubset(scored_df):
        raise ValueError("Scored pixels provided for clustering are not annotated")

    scored_df = scored_df.copy()
    if (
        not assigned_regions_name in scored_df.columns
    ):  # If input scores are not annotated by regions:
        logging.warning(
            f"No regions assigned to the scored pixels before clustering, using chromosomes"
        )
        scored_df[assigned_regions_name] = np.where(
            scored_df["chrom1"] == scored_df["chrom2"], scored_df["chrom1"], np.nan
        )

    # cluster within each regions separately and accumulate the result:
    pixel_clust_list = []
    scored_pixels_by_region = scored_df.groupby(assigned_regions_name, observed=True)
    for region, _df in scored_pixels_by_region:
        logging.info(f"clustering enriched pixels in region: {region}")
        # Using genomic corrdinated for clustering, not bin_id
        pixel_clust = clust_2D_pixels(
            _df,
            threshold_cluster=dots_clustering_radius,
            bin1_id_name="start1",
            bin2_id_name="start2",
        )
        pixel_clust_list.append(pixel_clust)
    logging.info("Clustering is complete")

    # concatenate clustering results ...
    # indexing information persists here ...
    if not pixel_clust_list:
        logging.warning("No clusters found for any regions! Output will be empty")
        empty_output = pd.DataFrame(
            [],
            columns=list(scored_df.columns)
            + [
                assigned_regions_name + "1",
                assigned_regions_name + "2",
                "c_label",
                "c_size",
                "cstart1",
                "cstart2",
            ],
        )
        return empty_output  # Empty dataframe with the same columns as anticipated
    else:
        pixel_clust_df = pd.concat(
            pixel_clust_list, ignore_index=False
        )  # Concatenate the clustering results for different regions

    # now merge pixel_clust_df and scored_df DataFrame ...
    # TODO make a more robust merge here
    df = pd.merge(
        scored_df, pixel_clust_df, how="left", left_index=True, right_index=True
    )
    # TODO check if next str-cast is neccessary
    df[assigned_regions_name + "1"] = df[assigned_regions_name].astype(str)
    df[assigned_regions_name + "2"] = df[assigned_regions_name].astype(str)

    return df



# Modified dots() function in dotfinder.py
# https://github.com/open2c/cooltools/blob/48bec4008d765b9e0ee0a72782a04d4bca73e930/cooltools/api/dotfinder.py

def generate_clusters_of_enriched_pixels(
    clr,
    expected,
    expected_value_col="balanced.avg",
    clr_weight_name="weight",
    view_df=None,
    kernels=None,
    max_loci_separation=10_000_000,
    max_nans_tolerated=1,  # test if this has desired behavior
    n_lambda_bins=40,  # update this eventually
    lambda_bin_fdr=0.1,
    clustering_radius=20_000,
    cluster_filtering=None,
    tile_size=5_000_000,
    nproc=1,
):
    """
    Call dots on a cooler {clr}, using {expected} defined in regions specified
    in {view_df}.

    All convolution kernels specified in {kernels} will be all applied to the {clr},
    and statistical testing will be performed separately for each kernel. A convolutional
    kernel is a small squared matrix (e.g. 7x7) of zeros and ones
    that defines a "mask" to extract local expected around each pixel. Since the
    enrichment is calculated relative to the central pixel, kernel width should
    be an odd number >=3.

    Parameters
    ----------
    clr : cooler.Cooler
        A cooler with balanced Hi-C data.
    expected : DataFrame in expected format
        Diagonal summary statistics for each chromosome, and name of the column
        with the values of expected to use.
    expected_value_col : str
        Name of the column in expected that holds the values of expected
    clr_weight_name : str
        Name of the column in the clr.bins to use as balancing weights.
        Using raw unbalanced data is not supported for dot-calling.
    view_df : viewframe
        Viewframe with genomic regions, at the moment the view has to match the
        view used for generating expected. If None, generate from the cooler.
    kernels : { str:np.ndarray } | None
        A dictionary of convolution kernels to be used for calculating locally adjusted
        expected. If None the default kernels from HiCCUPS are going to be recommended
        based on the resolution of the cooler.
    max_loci_separation : int
        Miaximum loci separation for dot-calling, i.e., do not call dots for
        loci that are further than max_loci_separation basepair apart. default 10Mb.
    max_nans_tolerated : int
        Maximum number of NaNs tolerated in a footprint of every used kernel
        Adjust with caution, as large max_nans_tolerated, might lead to artifacts in
        pixels scoring.
    n_lambda_bins : int
        Number of log-spaced bins, where FDR-testing will be performed independently.
        TODO: generate lambda-bins on the fly based on the dynamic range of the data (i.e. maximum pixel count)
    lambda_bin_fdr : float
        False discovery rate (FDR) for multiple hypothesis testing BH-FDR procedure, applied per lambda bin.
    clustering_radius : None | int
        Cluster enriched pixels with a given radius. "Brightest" pixels in each group
        will be reported as the final dot-calls. If None, no clustering is performed.
    cluster_filtering : bool
        whether to apply additional filtering to centroids after clustering, using cluster_filtering_hiccups()
    tile_size : int
        Tile size for the Hi-C heatmap tiling. Typically on order of several mega-bases, and <= max_loci_separation.
        Controls tradeoff between memory consumption and speed of execution.
    nproc : int
        Number of processes to use for multiprocessing.

    Returns
    -------
    dots : pandas.DataFrame
        BEDPE-style dataFrame with genomic coordinates of called dots and additional annotations.

    Notes
    -----
    'clustering_radius' in Birch clustering algorithm corresponds to a
    double the clustering radius in the "greedy"-clustering used in HiCCUPS
    (to be tested).

    TODO describe sequence of processing steps

    """

    #### Generate viewframes ####
    if view_df is None:
        view_df = make_cooler_view(clr)
    else:
        try:
            _ = is_compatible_viewframe(
                view_df,
                clr,
                check_sorting=True,
                raise_errors=True,
            )
        except Exception as e:
            raise ValueError("view_df is not a valid viewframe or incompatible") from e

    # check balancing status
    if clr_weight_name:
        # check if cooler is balanced
        try:
            _ = is_cooler_balanced(clr, clr_weight_name, raise_errors=True)
        except Exception as e:
            raise ValueError(
                f"provided cooler is not balanced or {clr_weight_name} is missing"
            ) from e
    else:
        raise ValueError("calling dots on raw data is not supported.")

    # add checks to make sure cis-expected is symmetric
    # make sure provided expected is compatible
    try:
        _ = is_valid_expected(
            expected,
            "cis",
            view_df,
            verify_cooler=clr,
            expected_value_cols=[
                expected_value_col,
            ],
            raise_errors=True,
        )
    except Exception as e:
        raise ValueError("provided expected is not compatible") from e
    expected = expected.set_index(["region1", "region2", "dist"]).sort_index()

    # Prepare some parameters.
    binsize = clr.binsize
    loci_separation_bins = bp_to_bins(max_loci_separation, binsize)
    tile_size_bins = bp_to_bins(tile_size, binsize)

    # verify provided kernels or recommend them (HiCCUPS)...
    if kernels and is_compatible_kernels(kernels, binsize, max_nans_tolerated):
        warnings.warn(
            "Compatibility checks for 'kernels' are not fully implemented yet, use at your own risk"
        )
    else:
        # recommend them (default hiccups ones for now)
        kernels = recommend_kernels(binsize)
    # deduce kernel_width - overall footprint
    kernel_width = max(len(k) for k in kernels.values())  # 2*w+1
    kernel_half_width = int((kernel_width - 1) / 2)  # former w parameter

    # try to guess required lambda bins using "max" value of pixel counts
    # statistical: lambda-binning edges ...
    if not 40 <= n_lambda_bins <= 50:
        raise ValueError(f"Incompatible n_lambda_bins={n_lambda_bins}")
    BASE = 2 ** (1 / 3)  # very arbitrary - parameterize !
    ledges = np.concatenate(
        (
            [-np.inf],
            np.logspace(
                0,
                n_lambda_bins - 1,
                num=n_lambda_bins,
                base=BASE,
                dtype=np.float64,
            ),
            [np.inf],
        )
    )

    # list of tile coordinate ranges
    tiles = list(
        generate_tiles_diag_band(
            clr, view_df, kernel_half_width, tile_size_bins, loci_separation_bins
        )
    )

    # 1. Calculate genome-wide histograms of scores.
    time_start = time.perf_counter()
    gw_hist = scoring_and_histogramming_step(
        clr,
        expected,
        expected_value_col=expected_value_col,
        clr_weight_name=clr_weight_name,
        tiles=tiles,
        kernels=kernels,
        ledges=ledges,
        max_nans_tolerated=max_nans_tolerated,
        loci_separation_bins=loci_separation_bins,
        nproc=nproc,
    )
    elapsed_time = time.perf_counter() - time_start
    logging.info(f"Done building histograms in {elapsed_time:.3f} sec ...")

    # 2. Determine the FDR thresholds.
    threshold_df, qvalues = determine_thresholds(gw_hist, lambda_bin_fdr)
    logging.info("Determined thresholds for every lambda-bin ...")

    # 3. Filter using FDR thresholds calculated in the histogramming step
    time_start = time.perf_counter()
    filtered_pixels = scoring_and_extraction_step(
        clr,
        expected,
        expected_value_col=expected_value_col,
        clr_weight_name=clr_weight_name,
        tiles=tiles,
        kernels=kernels,
        ledges=ledges,
        thresholds=threshold_df,
        max_nans_tolerated=max_nans_tolerated,
        loci_separation_bins=loci_separation_bins,
        nproc=nproc,
        bin1_id_name="bin1_id",
        bin2_id_name="bin2_id",
    )
    elapsed_time = time.perf_counter() - time_start
    logging.info(f"Done extracting enriched pixels in {elapsed_time:.3f} sec ...")

    # 4. Post-processing
    logging.info(f"Begin post-processing of {len(filtered_pixels)} filtered pixels")
    logging.info("preparing to extract needed q-values ...")

    # annotate enriched pixels
    filtered_pixels_qvals = annotate_pixels_with_qvalues(filtered_pixels, qvalues)
    filtered_pixels_annotated = cooler.annotate(
        filtered_pixels_qvals, clr.bins()[["chrom", "start", "end"]], replace=True
    )
    if not clustering_radius:
        # TODO: make sure returned DataFrame has the same columns as "postprocessed_calls"
        # columns to return before-clustering
        output_cols = []
        output_cols += bedpe_required_cols
        output_cols += [
            observed_count_name,
        ]
        output_cols += [f"la_exp.{k}.value" for k in kernels]
        output_cols += [f"la_exp.{k}.qval" for k in kernels]
        return filtered_pixels_annotated[output_cols]

    # 4a. clustering
    # Clustering is done independently for every region, therefore regions must be assigned:
    filtered_pixels_annotated = assign_regions(filtered_pixels_annotated, view_df)
    clusters = generate_clusters(
        filtered_pixels_annotated,
        clustering_radius,
    ).reset_index(drop=True)

    return clusters




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
    
    dots_dict[input_cooler_paths.name[i]] = generate_clusters_of_enriched_pixels(
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
    filename=os.path.join(output_folder,os.path.basename(input_cooler_paths.loc[input_cooler_paths.name == i,'path'].copy().reset_index(drop=True)[0])).replace('.1000.mcool', '.5kb.clusters.tsv')
    dots_dict[i].to_csv(filename, sep='\t', index=False, header=True)
    dots_filenames.append(filename)

meta_dots_files = {
    'name' : list(dots_dict.keys()),
    'path' : dots_filenames
}


df = pd.DataFrame(meta_dots_files)                          
df.to_csv(os.path.join(output_folder,'clusters_meta_5kb.tsv'), sep='\t', index=False, header=True)
