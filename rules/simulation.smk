#### generate simulated data ####
# Simulated data are generated using SimPol 

#### subsample read counts ####
rule subsample_simulation_for_lrt_pause_escape:
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
        os.path.join("logs/simulation/lrt_pause_escape", "{param_id}_{lambda_exp}.log")
    output:
        rate_tbl = os.path.join("outputs/simulation/tables/lrt_pause_escape", "{param_id}_{lambda_exp}.RDS"),
        rnap_tbl = os.path.join("outputs/simulation/tables/lrt_pause_escape", "{param_id}_{lambda_exp}.csv"),
        bw = os.path.join("outputs/simulation/tables/lrt_pause_escape", "{param_id}_{lambda_exp}.bw")
    script:
        "../scripts/simulation/pause_escape/subsample_simulation_pause_escape.R"

rule subsample_simulation_for_lrt_pause_distribution:
    input:
        rds = os.path.join("outputs/simulation/data_fk", "{param_id}", "positions/position_matrix_400000.csv")
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
        os.path.join("logs/simulation/lrt_pause_distribution", "{param_id}_{lambda_exp}.log")
    output:
        rate_tbl = os.path.join("outputs/simulation/tables/lrt_pause_distribution", "{param_id}_{lambda_exp}.RDS"),
        bw = os.path.join("outputs/simulation/tables/lrt_pause_distribution", "{param_id}_{lambda_exp}.bw")
    script:
        "../scripts/simulation/pause_distribution/subsample_simulation_pause_distribution.R"

#### perform LRT ####
rule lrt_pause_escape:
    input:
        expand(os.path.join("outputs/simulation/tables/lrt_pause_escape", "{param_id}_{lambda_exp}.RDS"), param_id = metadata_lrt_params.param_id, lambda_exp = ["lrt_high", "lrt_median", "lrt_low"])
    params:
        helper = "scripts/unimod/helper_function_em_two_condition.R",
        table_dir = "outputs/simulation/tables/lrt_pause_escape",
        figure_dir = "outputs/simulation/figures/lrt_pause_escape"
    threads:1
    log:
        os.path.join("logs/simulation/lrt_visualization", "pause_escape.log")
    output:
        done = touch(os.path.join("indicator/simulation/lrt_visualization", "pause_escape.done"))
    script:
        "../scripts/simulation/pause_escape/lrt_rates.R"

rule lrt_pause_distribution:
    input:
        expand(os.path.join("outputs/simulation/tables/lrt_pause_distribution", "{param_id}_{lambda_exp}.RDS"), param_id = metadata_dist_params.param_id, lambda_exp = ["lrt_high", "lrt_median", "lrt_low"])
    params:
        helper = "scripts/unimod/helper_function_em_two_condition.R",
        table_dir = "outputs/simulation/tables/lrt_pause_distribution",
        figure_dir = os.path.join("outputs/simulation/figures/lrt_pause_distribution", "{lambda_exp}"),
        lambda_exp = "{lambda_exp}"
    threads:1
    log:
        os.path.join("logs/simulation/lrt_visualization", "{lambda_exp}", "pause_distribution.log")
    output:
        done = touch(os.path.join("indicator/simulation/lrt_visualization", "{lambda_exp}", "pause_distribution.done"))
    script:
        "../scripts/simulation/pause_distribution/lrt_fk.R"

#### visualize results ####
rule visualize_simulation_pause_escape:
    input:
        os.path.join("indicator/simulation/lrt_visualization", "pause_escape.done")
    params:
        table_dir = "outputs/simulation/tables/lrt_pause_escape",
        figure_dir = "outputs/simulation/figures/lrt_pause_escape",
        lambda_exp = "{lambda_exp}"
    threads:1
    log:
        os.path.join("logs/simulation/lrt_visualization", "{lambda_exp}", "visualize_pause_escape.log")
    output:
        done = touch(os.path.join("indicator/simulation/lrt_visualization", "{lambda_exp}", "visualize_pause_escape.done"))
    script:
        "../scripts/simulation/visualize_pause_escape_bw.R"

rule visualize_simulation_pause_distribution:
    input:
        os.path.join(os.path.join("indicator/simulation/lrt_visualization", "{lambda_exp}", "pause_distribution.done"))
    params:
        table_dir = "outputs/simulation/tables/lrt_pause_distribution",
        figure_dir = os.path.join("outputs/simulation/figures/lrt_pause_distribution", "{lambda_exp}"),
        lambda_exp = "{lambda_exp}"
    threads:1
    log:
        os.path.join("logs/simulation/lrt_visualization", "{lambda_exp}", "visualize_pause_distribution.log")
    output:
        done = touch(os.path.join("indicator/simulation/lrt_visualization", "{lambda_exp}", "visualize_pause_distribution.done"))
    script:
        "../scripts/simulation/visualize_pause_distribution_bw.R"