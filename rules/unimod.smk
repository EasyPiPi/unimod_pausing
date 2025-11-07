rule generate_counting_region:
    input:
        gtf = rules.get_human_transcripts.output.gtf,
        tsn = "outputs/read_dt/max_tsn_per_gene_{cell_line}.rds"
    params:
        tss_length = config["tss_length"], # parameter k
        tts_length = config["tts_length"], # parameter m
        gb_min_length = config["gb_min_length"], 
        gb_max_length = config["gb_max_length"], 
        dist_to_tss = config["gb_max_length"]
    threads:1
    log:
        os.path.join("logs/generate_counting_region", "{cell_line}.log")
    output:
        grng = "outputs/read_dt/granges_for_read_counting_{cell_line}.RData",
        tss = "outputs/read_dt/tss_{cell_line}.bed"
    script:
        "../scripts/unimod/generate_counting_region.R"

rule analyze_one_sample_pause_release:
    input:
        grng = rules.generate_counting_region.output.grng,
        bwp1_p3 = os.path.join("outputs/bigwig/p3", combine_wildcard + "_plus.bw"),
        bwm1_p3 = os.path.join("outputs/bigwig/p3", combine_wildcard + "_minus.bw")
    params:
        helper = "scripts/unimod/helper_function.R",
        em = "scripts/unimod/helper_function_em_pause_release.R",
        result_dir = os.path.join("outputs/within_sample", combine_wildcard, "pause_release")
    threads:1
    log:
        os.path.join("logs/analyze_one_sample", combine_wildcard + "_pause_release.log")
    output:
        rate_calibrated = os.path.join("outputs/within_sample", combine_wildcard, "pause_release", "rate_calibrated.csv"),
        rate_tbl = os.path.join("outputs/within_sample", combine_wildcard, "pause_release", "rate.RDS")
    script:
        "../scripts/unimod/analyze_one_sample_poisson_pause_release.R"

rule analyze_one_sample_steric_hindrance:
    input:
        grng = rules.generate_counting_region.output.grng,
        bwp1_p3 = os.path.join("outputs/bigwig/p3", combine_wildcard + "_plus.bw"),
        bwm1_p3 = os.path.join("outputs/bigwig/p3", combine_wildcard + "_minus.bw"),
        scale = "outputs/between_samples/table/scale_factor.csv"
    params:
        helper = "scripts/unimod/helper_function.R",
        em = "scripts/unimod/helper_function_em_steric_hindrance.R",
        sample_id = combine_wildcard,
        result_dir = os.path.join("outputs/within_sample", combine_wildcard, "steric_hindrance")
    threads:1
    log:
        os.path.join("logs/analyze_one_sample", combine_wildcard + "_steric_hindrance.log")
    output:
        rate_tbl = os.path.join("outputs/within_sample", combine_wildcard, "steric_hindrance", "rate.csv")
    script:
        "../scripts/unimod/analyze_one_sample_poisson_steric_hindrance.R"

rule compare_rates_across_samples:
    input:
        expand(os.path.join("outputs/within_sample", expand_combine_wildcard, "pause_release", "rate.RDS"), df = metadata.itertuples())
    threads:1
    log:
        os.path.join("logs/compare_rates_across_samples", "run.log")
    output:
        touch("indicator/compare_rates_across_samples/run.done")
    script:
        "../scripts/unimod/compare_rates_across_samples.R"

rule analyze_two_samples:
    input:
        rc1 = os.path.join("outputs/within_sample", "{assay}-{cell_line}-{reference}-control-{read_type}", "pause_release", "rate.RDS"),
        rc2 = os.path.join("outputs/within_sample", "{assay}-{cell_line}-{reference}-treated-{read_type}", "pause_release", "rate.RDS"),
        rate1 = os.path.join("outputs/within_sample", "{assay}-{cell_line}-{reference}-control-{read_type}", "steric_hindrance", "rate.csv"),
        rate2 = os.path.join("outputs/within_sample", "{assay}-{cell_line}-{reference}-treated-{read_type}", "steric_hindrance", "rate.csv"),
        spike_in = "metadata/scaling_factor.csv"
    params:
        quantile_normalization = "{normalization}",
        replicates = "{replicates}", # Whether using all loci or only gene bodies in the analyzed set to calculate lambda
        result_dir = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}"),
        helper_tc = "scripts/unimod/helper_function_em_two_condition.R",
        helper_pr = "scripts/unimod/helper_function_em_pause_release.R"
    threads:1
    log:
        os.path.join("logs/analyze_two_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}" + ".log")
    output:
        omega = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}", "omega.csv"),
        beta = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}", "beta.csv"),
        alpha = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}", "alpha.csv")
    script:
        "../scripts/unimod/analyze_two_samples.R"

rule visualize_two_samples:
    input:
        grng = rules.generate_counting_region.output.grng,
        gtf = rules.get_human_transcripts.output.gtf,
        spike_in = "metadata/scaling_factor.csv",
        bwp1_p3 = os.path.join("outputs/bigwig/p3", "{assay}-{cell_line}-{reference}-control-{read_type}" + "_plus.bw"),
        bwm1_p3 = os.path.join("outputs/bigwig/p3", "{assay}-{cell_line}-{reference}-control-{read_type}" + "_minus.bw"),
        bwp2_p3 = os.path.join("outputs/bigwig/p3", "{assay}-{cell_line}-{reference}-treated-{read_type}" + "_plus.bw"),
        bwm2_p3 = os.path.join("outputs/bigwig/p3", "{assay}-{cell_line}-{reference}-treated-{read_type}" + "_minus.bw"),
        beta = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}", "beta.csv"),
        omega = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}", "omega.csv"),
        alpha = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}", "alpha.csv"),
        tf_target = "outputs/between_samples/table/gsea_targets.csv"
    params:
        result_dir = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}")
    threads:1
    log:
        os.path.join("logs/visualize_two_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}" + ".log")
    output:
        touch(os.path.join("indicator/visualize_two_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}" + ".done"))
    script:
        "../scripts/unimod/visualize_two_samples.R"

# rule do_gsea:
#     input:
#         omega = rules.analyze_two_samples.output.omega,
#         beta = rules.analyze_two_samples.output.beta
#     params:
#         prop = 0.1, # proportion of genes to be considered in a GSEA enrichment analysis,
#         result_dir = os.path.join("outputs/between_samples", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}", "go")
#     threads:2
#     log:
#         os.path.join("logs/do_gsea", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}" + ".log")
#     output:
#         touch(os.path.join("indicator/do_gsea", "{assay}-{cell_line}-{reference}-{read_type}-{normalization}-{replicates}" + ".complete"))
#     script:
#         "../scripts/go/do_gsea_within_study.R"
