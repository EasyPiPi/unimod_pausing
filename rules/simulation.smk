#### run simulation ####
# now run simulation on cluster, may need a better implementation afterwards

#### summarize results ####
rule subsample_simulation_for_pause_release:
    input:
        rds = os.path.join("outputs/simulation/data", "{param_id}.RDS")
    params:
        helper = "scripts/unimod/helper_function_em_pause_release.R",
        sample_cell = 5000, # number of cells for subsampling,
        sample_n = 50, # number of times to sample
        lambda_exp = "median",
        matched_len = 20000,
        sel_sample = "{param_id}",
        count_rnap = True # whether to count RNAP number before pause site or not
    threads:1
    log:
        os.path.join("logs/simulation/rate/pause_release", "{param_id}.log")
    output:
        rate_tbl = os.path.join("outputs/simulation/tables/rate/pause_release", "{param_id}.RDS"),
        rnap_tbl = os.path.join("outputs/simulation/tables/rate/pause_release", "{param_id}.csv")
    script:
        "../scripts/simulation/subsample_simulation_pause_release.R"

rule subsample_simulation_for_steric_hindrance:
    input:
        rds = os.path.join("outputs/simulation/data", "{param_id}.RDS")
    params:
        helper = "scripts/unimod/helper_function_em_steric_hindrance.R",
        sample_cell = 5000, # number of cells for subsampling,
        sample_n = 50, # number of times to sample
        lambda_exp = "median",
        matched_len = 20000,
        sel_sample = "{param_id}"
    threads:1
    log:
        os.path.join("logs/simulation/rate/steric_hindrance", "{param_id}.log")
    output:
        rate_tbl = os.path.join("outputs/simulation/tables/rate/steric_hindrance", "{param_id}.RDS"),
        rnap_tbl = os.path.join("outputs/simulation/tables/rate/steric_hindrance", "{param_id}.csv")
    script:
        "../scripts/simulation/subsample_simulation_steric_hindrance.R"

rule summarize_simulation_for_pause_release:
    input:
        expand(os.path.join("outputs/simulation/tables/rate/pause_release", "{param_id}.RDS"), param_id = metadata_pr_params.param_id),
        meta_pr = "metadata/simulation_params_pause_release.csv",
        meta_st = "metadata/simulation_params_steric_hindrance.csv"
    params:
        table = "outputs/simulation/tables/rate/pause_release",
        figure = "outputs/simulation/figures/pause_release_spacing50_varied_zeta"
    threads:1
    log:
        os.path.join("logs/simulation/rate/pause_release", "all.log")
    output:
        touch("indicator/simulation/pause_release.done")
    script:
        "../scripts/simulation/summarize_simulation_subsamples_pause_release.R"

rule summarize_simulation_for_steric_hindrance:
    input:
        expand(os.path.join("outputs/simulation/tables/rate/steric_hindrance", "{param_id}.RDS"), param_id = metadata_st_params.param_id),
        meta_st = "metadata/simulation_params_steric_hindrance.csv"
    params:
        table = "outputs/simulation/tables/rate/steric_hindrance",
        figure = "outputs/simulation/figures/steric_hindrance_varied_zeta"
    threads:1
    log:
        os.path.join("logs/simulation/rate/steric_hindrance", "all.log")
    output:
        touch("indicator/simulation/steric_hindrance.done")
    script:
        "../scripts/simulation/summarize_simulation_subsamples_steric_hindrance.R"

rule subsample_simulation_for_lrts_with_matched_coverage:
    input:
        rds = os.path.join("outputs/simulation/data_lrt", "{param_id}_pos.RDS")
    params:
        helper = "scripts/unimod/helper_function_em_pause_release.R",
        sample_cell = 5000, # number of cells for subsampling,
        sample_n = 100, # number of times to sample
        lambda_exp = "{lambda_exp}", # scaling factor to match simulation to coverage in experimental data
        matched_len = 20000,
        sel_sample = "{param_id}",
        count_rnap = False
    threads:1
    wildcard_constraints:
       param_id="[^_]+"
    log:
        os.path.join("logs/simulation/lrt_matched_cov", "{param_id}_{lambda_exp}.log")
    output:
        rate_tbl = os.path.join("outputs/simulation/tables/lrt_matched_cov", "{param_id}_{lambda_exp}.RDS"),
        rnap_tbl = os.path.join("outputs/simulation/tables/lrt_matched_cov", "{param_id}_{lambda_exp}.csv")
    script:
        "../scripts/simulation/subsample_simulation_pause_release.R"
