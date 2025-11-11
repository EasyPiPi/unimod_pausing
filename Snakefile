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
metadata_procap = pd.read_csv("metadata/copro_sample.csv", dtype=str)
metadata_comparison = pd.read_csv("metadata/metadata_comparison.csv", dtype=str)

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
        #### Experimental Data ####
        ## Preprocessing ##
        # PRO-seq 
        "indicator/proseq/all.done",
        # PRO-cap liftover 
        expand(os.path.join("ext_data/copro/hg38", "{sample}" + ".bw"), sample = metadata_procap.file_name),
        ## find active TSSs ##
        "outputs/read_dt/human_transcript_granges.rds",
        ## Rate estimates ##
        expand(os.path.join("outputs/within_sample", expand_combine_wildcard, "pause_release", "rate.RDS"), df = metadata_aoi.itertuples()),
        ## LRT ##
        expand(os.path.join("outputs/between_samples", "{df.group_1}" + "_vs_" + "{df.group_2}", "omega.csv"), df = metadata_comparison.itertuples()),
        expand(os.path.join("indicator/visualize_two_samples", "{df.group_1}" + "_vs_" + "{df.group_2}" + ".done"), df = metadata_comparison.itertuples()),
        # #### Simulation ####
        # "indicator/simulation/pause_release.done",
        # "indicator/simulation/steric_hindrance.done",
        # ## LRT ##
        # expand(os.path.join("outputs/simulation/tables/lrt_matched_cov", "{param_id}_{lambda_exp}.RDS"), param_id = lrt_params.param_id, lambda_exp = ["lrt_high", "lrt_median", "lrt_low"]),


##### load rules #####
include: "rules/proseq.smk"
include: "rules/unimod.smk"
# include: "rules/simulation.smk"
