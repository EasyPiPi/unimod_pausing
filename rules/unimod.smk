rule generate_counting_region:
    input:
        gtf = rules.get_human_transcripts.output.gtf,
        tsn = "outputs/read_dt/max_tsn_per_gene_{cell_line}.rds"
    params:
        tss_length = config["tss_length"], # parameter k
        tts_length = config["tts_length"], # parameter m
        gb_min_length = config["gb_min_length"], 
        gb_max_length = config["gb_max_length"], 
        dist_to_tss = config["dist_to_tss"]
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

rule analyze_two_samples:
    input:
        rc1 = os.path.join("outputs/within_sample", "PROseq-DLD1-aoi-{group_1}-SE", "pause_release", "rate.RDS"),
        rc2 = os.path.join("outputs/within_sample", "PROseq-DLD1-aoi-{group_2}-SE", "pause_release", "rate.RDS"),
        spike_in = "metadata/scaling_factor.csv"
    params:
        helper_tc = "scripts/unimod/helper_function_em_two_condition.R",
        helper_pr = "scripts/unimod/helper_function_em_pause_release.R",
        result_dir = os.path.join("outputs/between_samples", "{group_1}" + "_vs_" + "{group_2}")
    threads:1
    log:
        os.path.join("logs/analyze_two_samples", "{group_1}" + "_vs_" + "{group_2}" + ".log")
    output:
        omega = os.path.join("outputs/between_samples", "{group_1}" + "_vs_" + "{group_2}", "omega.csv"),
        beta = os.path.join("outputs/between_samples", "{group_1}" + "_vs_" + "{group_2}", "beta.csv")
    script:
        "../scripts/unimod/analyze_two_samples_DLD1.R"

rule visualize_two_samples:
    input:
        grng = "outputs/read_dt/granges_for_read_counting_DLD1.RData",
        gtf = rules.get_human_transcripts.output.gtf,
        spike_in = "metadata/scaling_factor.csv",
        bwp1_p3 = os.path.join("outputs/bigwig/p3", "PROseq-DLD1-aoi-{group_1}-SE" + "_plus.bw"),
        bwm1_p3 = os.path.join("outputs/bigwig/p3", "PROseq-DLD1-aoi-{group_1}-SE" + "_minus.bw"),
        bwp2_p3 = os.path.join("outputs/bigwig/p3", "PROseq-DLD1-aoi-{group_2}-SE" + "_plus.bw"),
        bwm2_p3 = os.path.join("outputs/bigwig/p3", "PROseq-DLD1-aoi-{group_2}-SE" + "_minus.bw"),
        beta = os.path.join("outputs/between_samples", "{group_1}" + "_vs_" + "{group_2}", "beta.csv"),
        omega = os.path.join("outputs/between_samples", "{group_1}" + "_vs_" + "{group_2}", "omega.csv")
    params:
        result_dir = os.path.join("outputs/between_samples", "{group_1}" + "_vs_" + "{group_2}")
    threads:1
    log:
        os.path.join("logs/visualize_two_samples", "{group_1}" + "_vs_" + "{group_2}" + "visualization.log")
    output:
        touch(os.path.join("indicator/visualize_two_samples", "{group_1}" + "_vs_" + "{group_2}" + ".done"))
    script:
        "../scripts/unimod/visualize_two_samples.R"
