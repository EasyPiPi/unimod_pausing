#### config file ####
configfile: "config.yml"

#### import packages ####
import os
import pandas as pd
import numpy as np
import re
import functools
import operator

# helper for using multiindex to select rows
idx = pd.IndexSlice

##### metadata #####
sel_col = ['assay', 'cell_line', 'reference', 'group', 'read_type', 'replicate']
# metadata for samples
def make_metadata(file_path, selected_column = sel_col):
    metadata = pd.read_csv(file_path, dtype=str)
    metadata.set_index(selected_column, drop = False, inplace = True)
    metadata.sort_index(inplace = True)
    return metadata

#### metadata ####
metadata_aoi = make_metadata("metadata/metadata_aoi.csv")
croprodata = pd.read_csv("metadata/copro_sample.csv", dtype=str)

# simulations for LRT
lrt_params = pd.read_csv("metadata/simulation_params_lrt.csv", dtype=str)

#### store wildcard for easier access ####
sample_wildcard = "{assay}-{cell_line}-{reference}-{group}-{read_type}-{replicate}"
expand_sample_wildcard = "{df.assay}-{df.cell_line}-{df.reference}-{df.group}-{df.read_type}-{df.replicate}"
combine_wildcard = "{assay}-{cell_line}-{reference}-{group}-{read_type}"
expand_combine_wildcard = "{df.assay}-{df.cell_line}-{df.reference}-{df.group}-{df.read_type}"

#### rules ####
rule all:
    input:
        #### Prepocessing ####
        ## PRO-seq ##
        "indicator/proseq/all.done",
        # ## coPRO-seq ##
        # expand(os.path.join("ext_data/copro/hg38", "{sample}" + ".bw"), sample = croprodata.file_name),
        # # ## find active TSSs ##
        # # "outputs/read_dt/human_transcript_granges.rds",
        # #### Rate estimates ####
        # ## Simulation ##
        # "indicator/simulation/pause_release.done",
        # "indicator/simulation/steric_hindrance.done",
        # ## Experiment - replicate ##
        # expand(os.path.join("outputs/within_replicate", expand_sample_wildcard, "pause_release", "rate.RDS"), df = metadata.itertuples()),
        # ## Experiment - sample ##
        # expand(os.path.join("outputs/within_sample", expand_combine_wildcard, "pause_release", "rate.RDS"), df = metadata.itertuples()),
        # expand(os.path.join("outputs/within_sample", expand_combine_wildcard, "pause_release", "rate.RDS"), df = metadata_aoi.itertuples()),
        # expand(os.path.join("outputs/within_sample", expand_combine_wildcard, "steric_hindrance", "rate.csv"), df = metadata.itertuples()),
        # # compare transcription rate and pause sites across samples #
        # "indicator/compare_rates_across_samples/run.done",
        # #### LRT ####
        # ## Simulation ##
        # expand(os.path.join("outputs/simulation/tables/lrt_matched_cov", "{param_id}_{lambda_exp}.RDS"), param_id = lrt_params.param_id, lambda_exp = ["lrt_high", "lrt_median", "lrt_low"]),
        # ## Experiment ##
        # expand(os.path.join("outputs/between_samples", "{df.assay}-{df.cell_line}-{df.reference}-{df.read_type}-{normalization}-{replicates}", "omega.csv"), df = metadata.itertuples(), normalization = ["identity"], replicates = ["all"]),
        # expand(os.path.join("indicator/visualize_two_samples", "{df.assay}-{df.cell_line}-{df.reference}-{df.read_type}-{normalization}-{replicates}" + ".done"), df = metadata.itertuples(), normalization = ["identity"], replicates = ["all"]),
        # ## GSEA analyses ##
        # # expand(os.path.join("indicator/do_gsea", "{df.assay}-{df.cell_line}-{df.reference}-{df.read_type}-{normalization}-{replicates}" + ".complete"), df = metadata.itertuples(), normalization = ["identity"], replicates = ["all"]),

##### load rules #####
include: "rules/proseq.smk"
# include: "rules/unimod.smk"
# include: "rules/simulation.smk"
