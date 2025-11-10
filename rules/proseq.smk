#### build bwa index ####
rule compute_chrominfo:
    input:
    # download from https://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/
        genome = os.path.join("ext_data/genome/dm6", "Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz")
    output:
        genome = os.path.join("ext_data/genome/dm6", "Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"),
        # https://www.biostars.org/p/272373/
        chrominfo = os.path.join("ext_data/genome/dm6", "chrominfo.txt")
    shell:
        """
        gzip -dk {input.genome}
        samtools faidx {output.genome}
        cut -f1,2 {output.genome}.fai > {output.chrominfo}
        """

rule build_bwa_index:
    input:
        genome = os.path.join(config["ROOT_DIR"], "ext_data/genome/dm6", "Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz")
    output:
        complete = touch("indicator/bwa/index/dm6.complete")
    params:
        dir = "ext_data/genome/dm6/bwa/index"
    shell:
    # consider using realpath to find the relative path therefore we can copy symlink there
    # e.g., test_dir=$(realpath --relative-to=. ~/projects/comparaReg/external_resources/baboon/genome.fasta)
    # or simply use ln -rfs ~/tmp/abs.txt abs
        """
        mkdir -p {params.dir}
        cd {params.dir}
        ln -s {input.genome} genome.fasta
        bwa index genome.fasta
        """

# for single end libraries
def get_fastq_se(wildcards):
    # need to be changed when more datasets are included
    if (wildcards.assay == "PROseq") & (wildcards.reference in ["aoi"]):
        subdf = metadata_aoi.loc[[tuple(wildcards)], :]
        r1 = os.path.join(config["GEO_PROSEQ_DIR"], subdf.file_name.values[0] + ".fastq.gz")
    return r1

rule make_symlink_se:
    input:
        r1 = get_fastq_se
    output:
        r1 = temp(os.path.join("tmp/raw_data/p3", sample_wildcard + ".fastq.gz"))
    wildcard_constraints:
        # need to be changed when more datasets are included
        reference = "aoi|hoyt"
    shell:
        """
        ln -s {input.r1} {output.r1}
        """

rule process_proseq_se:
    input:
        r1 = rules.make_symlink_se.output.r1,
        chrominfo = config["HUMAN_CHROMINFO"],
        index = config["HUMAN_BWA_INDEX"]
    output:
        complete = touch(os.path.join("indicator/proseq/p3", sample_wildcard + ".complete"))
    params:
        others = "-SE -P",
        prefix = os.path.join("tmp/raw_data/p3", sample_wildcard),
        tmp_dir = os.path.join("tmp/proseq/p3", sample_wildcard),
        output = os.path.join("outputs/proseq/p3", sample_wildcard)
    wildcard_constraints:
        assay = "PROseq",
        read_type = "SE"
    threads:4
    shell:
        """
        mkdir -p {params.tmp_dir}
        bash scripts/proseq/proseq2.0.bsh {params.others} -i {input.index} -c {input.chrominfo} -I {params.prefix} -T {params.tmp_dir} -O {params.output} --thread={threads}
        rm -rf {params.tmp_dir}
        """

rule process_proseq_se_dme:
    input:
        r1 = rules.make_symlink_se.output.r1,
        complete = rules.build_bwa_index.output.complete,
        chrominfo = rules.compute_chrominfo.output.chrominfo,
        index = 'ext_data/genome/dm6/bwa/index/genome.fasta'
    output:
        complete = touch(os.path.join("indicator/proseq/dme/p3", sample_wildcard + ".complete"))
    params:
        others = "-SE -P",
        prefix = os.path.join("tmp/raw_data/p3", sample_wildcard),
        tmp_dir = os.path.join("tmp/proseq/dme/p3", sample_wildcard),
        output = os.path.join("outputs/proseq/dme/p3", sample_wildcard)
    wildcard_constraints:
        assay = "PROseq",
        read_type = "SE"
    threads:4
    shell:
        """
        mkdir -p {params.tmp_dir}
        bash scripts/proseq/proseq2.0.bsh {params.others} -i {input.index} -c {input.chrominfo} -I {params.prefix} -T {params.tmp_dir} -O {params.output} --thread={threads}
        rm -rf {params.tmp_dir}
        """

def get_bam_files(wc):
    proseq_dir_pattern = "-".join([wc.assay, wc.cell_line, wc.reference, wc.group, wc.read_type, wc.replicate])
    files = os.listdir(os.path.join("outputs/proseq/p3", proseq_dir_pattern))
    # print(files)
    r = re.compile(".*.bam")
    bam_file = list(filter(r.match, files))[0]
    # print(bam_file)
    return os.path.join("outputs/proseq/p3", proseq_dir_pattern, bam_file)

rule extact_5prime_end_proseq_se:
    input:
        complete = rules.process_proseq_se.output.complete
    output:
        bed = os.path.join("tmp/proseq/p5", sample_wildcard, sample_wildcard + ".bed.gz")
    params:
        bam = os.path.join("outputs/proseq/p3", sample_wildcard, sample_wildcard + "_QC.sort.bam"),
        tmp = os.path.join("tmp/proseq/p5", sample_wildcard),
        prefix = sample_wildcard
    wildcard_constraints:
        assay = "PROseq",
        read_type = "SE"
    shell:
        """
        mkdir -p {params.tmp}
        bedtools bamtobed -i {params.bam} 2> {params.tmp}/kill.warnings | awk 'BEGIN{{OFS="\\t"}} ($5 > 0) {{print $0}}' | awk 'BEGIN{{OFS="\\t"}} ($6 == "+") {{print $1,$3-1,$3,$4,$5,"-"}}; ($6 == "-") {{print $1,$2,$2+1,$4,$5,"+"}}' | gzip > {output.bed}
        """

rule normalize_5prime_end_data:
    input:
        chrominfo = config["HUMAN_CHROMINFO"],
        bed = os.path.join("tmp/proseq/p5", sample_wildcard, sample_wildcard + ".bed.gz")
    output:
        bwp = os.path.join("outputs/proseq/p5", sample_wildcard, sample_wildcard + "_plus.bw"),
        bwm = os.path.join("outputs/proseq/p5", sample_wildcard, sample_wildcard + "_minus.bw"),
        bwp_rpm = os.path.join("outputs/proseq/p5", sample_wildcard, sample_wildcard + "_plus.rpm.bw"),
        bwm_rpm = os.path.join("outputs/proseq/p5", sample_wildcard, sample_wildcard + "_minus.rpm.bw")
    params:
        tmp = os.path.join("tmp/proseq/p5", sample_wildcard),
        prefix = sample_wildcard
    shell:
        """
        zcat {input.bed} | grep "rRNA\|chrM" -v | grep "_" -v | sort-bed - | gzip > {params.tmp}/{params.prefix}.nr.rs.bed.gz
        ## Convert to bedGraph ... Can't gzip these, unfortunately.
        bedtools genomecov -bg -i {params.tmp}/{params.prefix}.nr.rs.bed.gz -g {input.chrominfo} -strand + > {params.tmp}/{params.prefix}_plus.bedGraph
        bedtools genomecov -bg -i {params.tmp}/{params.prefix}.nr.rs.bed.gz -g {input.chrominfo} -strand - > {params.tmp}/{params.prefix}_minus.noinv.bedGraph
        ## Invert read counts on the minus strand.
        cat {params.tmp}/{params.prefix}_minus.noinv.bedGraph | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,-1*$4}}' > {params.tmp}/{params.prefix}_minus.bedGraph
        ## normalized by RPM
        readCount=`zcat {input.bed} | grep "" -c`
        cat {params.tmp}/{params.prefix}_plus.bedGraph | awk 'BEGIN{{OFS="\\t"}} {{print $1,$2,$3,$4*1000*1000/'$readCount'/1}}' > {params.tmp}/{params.prefix}_plus.rpm.bedGraph
        cat {params.tmp}/{params.prefix}_minus.bedGraph | awk 'BEGIN{{OFS="\t"}} {{print $1,$2,$3,$4*1000*1000/'$readCount'/1}}' > {params.tmp}/{params.prefix}_minus.rpm.bedGraph
        ## Convert to bigwig (RPM)
        bedGraphToBigWig {params.tmp}/{params.prefix}_plus.rpm.bedGraph {input.chrominfo} {output.bwp_rpm}
        bedGraphToBigWig {params.tmp}/{params.prefix}_minus.rpm.bedGraph {input.chrominfo} {output.bwm_rpm}
        ## Convert to bigwig
        bedGraphToBigWig {params.tmp}/{params.prefix}_plus.bedGraph {input.chrominfo} {output.bwp}
        bedGraphToBigWig {params.tmp}/{params.prefix}_minus.bedGraph {input.chrominfo} {output.bwm}
        ## cleanup
        rm -rf {params.tmp}
        """

# retrieve bigwig files, used by both rule merge_bigwig and rule analyze_one_sample_pause_release_replicate
def get_bigwig_files(wc, strand):
    try:
        proseq_dir = os.path.join("outputs/proseq", wc.prime)
    except:
        proseq_dir = os.path.join("outputs/proseq", "p3")
    try:
        proseq_dir_pattern = "-".join([wc.assay, wc.cell_line, wc.reference, wc.group, wc.read_type, wc.replicate])
    except:
        proseq_dir_pattern = "-".join([wc.assay, wc.cell_line, wc.reference, wc.group, wc.read_type])
    # select dirs based on wildcards
    r = re.compile(proseq_dir_pattern)
    sel_dir = list(filter(r.match, os.listdir(proseq_dir)))
    sel_dir = [os.path.join(proseq_dir, x) for x in sel_dir]
    # list files in selected dirs
    sel_file = [os.listdir(x) for x in sel_dir]
    sel_file = functools.reduce(operator.iconcat, sel_file, []) # flatten list
    # only keep bigwig files
    r = re.compile(".*" + strand + ".bw")
    bw_file = list(filter(r.match, sel_file))
    bw_file = [os.path.join(x, y) for x, y in zip(sel_dir, bw_file)]
    return bw_file

# combine PRO-seq samples
rule merge_bigwig:
    input:
        expand(os.path.join("outputs/proseq/p5", expand_sample_wildcard, expand_sample_wildcard + "_plus.bw"), df = metadata_aoi.itertuples()),
        expand(os.path.join("indicator/proseq/p3", expand_sample_wildcard + ".complete"), df = metadata_aoi.itertuples()),
        chrominfo = config["HUMAN_CHROMINFO"],
        bw_plus = lambda wc: get_bigwig_files(wc, strand = "plus"),
        bw_minus = lambda wc: get_bigwig_files(wc, strand = "minus")           
    params:
        threshold = 1e9
    output:
        bedgraph_plus = temp(os.path.join('tmp/proseq/{prime}', combine_wildcard + '_plus_sorted.bedGraph')),
        bedgraph_minus = temp(os.path.join('tmp/proseq/{prime}', combine_wildcard + '_minus_sorted.bedGraph')),
        bw_plus = os.path.join("outputs/bigwig/{prime}", combine_wildcard + "_plus.bw"),
        bw_minus = os.path.join("outputs/bigwig/{prime}", combine_wildcard + "_minus.bw")
    run:
        # specify how to handle studies without replicates
        # if wildcards.reference not in ["chivu", "aoi"]:
        shell(
            """
            touch {output.bedgraph_plus} {output.bedgraph_minus}
            ln -sfr {input.bw_plus} {output.bw_plus}
            ln -sfr {input.bw_minus} {output.bw_minus}
            """
        )

# expand to generate read counts for both 5 prime and 3 prime end
rule proseq_done:
    input:
        "indicator/bwa/index/dm6.complete",
        # expand(os.path.join("indicator/proseq/dme/p3", expand_sample_wildcard + ".complete"), df = metadata_aoi.itertuples()),
        expand(os.path.join("outputs/bigwig/{prime}", expand_combine_wildcard + "_plus.bw"), df = metadata_aoi.itertuples(), prime = ["p3", "p5"])
    output:
        touch("indicator/proseq/all.done")

# crossmap PRO-cap data
rule crossmap_copro:
    input:
        bw = os.path.join("ext_data/copro/hg19", "{sample}" + ".bw"),
        chain = os.path.join("ext_data/chain", "hg19ToHg38.over.chain.gz")
    output:
        bw = os.path.join("ext_data/copro/hg38", "{sample}" + ".bw")
    params:
        bw = os.path.join("ext_data/copro/hg38", "{sample}")
    shell:
        "CrossMap bigwig {input.chain} {input.bw} {params.bw}"

# find active TSSs based on PRO-cap data and plot read counts
rule get_human_transcripts:
    input:
        non_olp_gn = config["HUMAN_NONOLP_GENE"],
        gtf = config["HUMAN_GTF"]
    log:
        os.path.join("logs/get_human_transcripts", "run.log")
    output:
        gtf = "outputs/read_dt/human_transcript_granges.rds"
    script:
        "../scripts/proseq/get_human_transcripts.R"

rule find_active_tss_DLD1:
    input:
        bwp_rep1 = os.path.join("ext_data/copro/hg38", "GSM4296337_PRO-cap-01-NELFC-AID-untreated-rep1.plus.bw"),
        bwm_rep1 = os.path.join("ext_data/copro/hg38", "GSM4296337_PRO-cap-01-NELFC-AID-untreated-rep1.minus.bw"),
        bwp_rep2 = os.path.join("ext_data/copro/hg38", "GSM4296339_PRO-cap-03-NELFC-AID-untreated-rep2.plus.bw"),
        bwm_rep2 = os.path.join("ext_data/copro/hg38", "GSM4296339_PRO-cap-03-NELFC-AID-untreated-rep2.minus.bw"),
        gtf = rules.get_human_transcripts.output.gtf
    params:
        rc_cutoff = config["rc_cutoff"]
    threads:12
    log:
        os.path.join("logs/find_active_tss_DLD1", "run.log")
    output:
        max_tsn_gn = "outputs/read_dt/max_tsn_per_gene_DLD1.rds"
    script:
        "../scripts/proseq/find_active_tss_DLD1.R"
