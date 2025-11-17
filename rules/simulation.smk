#### generate simulated data ####
# Simulated data are generated using SimPol 

#### summarize results ####
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
