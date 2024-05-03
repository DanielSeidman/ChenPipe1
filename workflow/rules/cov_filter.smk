rule compute_d4:
    input:
        unpack(get_bams)
    output:
        "results/{ref_name}/callable_sites/{sample}.mosdepth.global.dist.txt",
        temp("results/{ref_name}/callable_sites/{sample}.per-base.d4"),
        summary="results/{ref_name}/callable_sites/{sample}.mosdepth.summary.txt"
    conda:
        "../envs/cov_filter.yml"
    log:
        "logs/{ref_name}/compute_d4/{sample}.txt"
    benchmark:
        "benchmarks/{ref_name}/compute_d4/{sample}.txt"
    params:
        prefix = os.path.join(DEFAULT_STORAGE_PREFIX, "results/{ref_name}/callable_sites/{sample}")
    shell:
        "mosdepth --d4 -t {threads} {params.prefix} {input.bam} &> {log}"

rule merge_d4:
    input:
        unpack(get_input_for_coverage)
    output:
        "results/{ref_name}/callable_sites/all_samples.d4"
    conda:
        "../envs/cov_filter.yml"
    log:
        "logs/{ref_name}/merge_d4/log.txt"
    benchmark:
        "benchmarks/{ref_name}/merge_d4/benchmark.txt"
    shell:
        "d4tools merge {input.d4files} {output} &> {log}"

rule collect_covstats:
    input:
        unpack(get_input_covstats)
    output:
        "results/{ref_name}/summary_stats/all_cov_sumstats.txt"  
    run:
        covStats = collectCovStats(input.covStatFiles)
        with open(output[0], "w") as f:
            print("chrom\tmean_cov\tstdev_cov", file=f)
            for chrom in covStats:
                print(chrom, covStats[chrom]['mean'], covStats[chrom]['stdev'], sep="\t", file=f)

rule create_cov_bed:
    input:
        stats = "results/{ref_name}/summary_stats/all_cov_sumstats.txt",
        d4 = "results/{ref_name}/callable_sites/all_samples.d4"
    output:
        covbed = "results/{ref_name}/callable_sites/{prefix}_callable_sites_cov.bed"
    benchmark:
        "benchmarks/{ref_name}/covbed/{prefix}_benchmark.txt"
    params:
        cov_threshold_stdev = config["cov_threshold_stdev"],
        cov_threshold_lower = config["cov_threshold_lower"],
        cov_threshold_upper = config["cov_threshold_upper"],
        cov_threshold_rel = config["cov_threshold_rel"]
    conda:
        "../envs/cov_filter.yml"
    script:
        "../scripts/create_coverage_bed.py"

rule callable_bed:
    input:
        cov = "results/{ref_name}/callable_sites/{prefix}_callable_sites_cov.bed",
        map = "results/{ref_name}/callable_sites/{prefix}_callable_sites_map.bed"
    output:
        callable_sites = "results/{ref_name}/{prefix}_callable_sites.bed",
        tmp_cov = temp("results/{ref_name}/callable_sites/{prefix}_temp_cov.bed")
    conda:
        "../envs/cov_filter.yml"
    benchmark:
        "benchmarks/{ref_name}/callable_bed/{prefix}_benchmark.txt"
    params:
        merge = config['cov_merge']
    shell:
        """
        bedtools sort -i {input.cov} | bedtools merge -d {params.merge} -i - > {output.tmp_cov}
        bedtools intersect -a {output.tmp_cov} -b {input.map} | bedtools sort -i - | bedtools merge -i - > {output.callable_sites}
        """
