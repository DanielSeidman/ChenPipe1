rule picard_intervals:
    input:
        ref = "config/{ref_name}.fasta",
        fai = "config/{ref_name}.fasta.fai",
        dictf = "config/{ref_name}.dict"
    output:
        intervals = temp("results/{ref_name}/intervals/picard_interval_list.list")
    params:
        minNmer = int(config['minNmer'])
    conda:
        '../envs/bam2vcf.yml'
    log:
        "logs/{ref_name}/picard_intervals/log.txt"
    benchmark:
        "benchmarks/{ref_name}/picard_intervals/benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['process_ref']['threads'] * 4000
    shell:
        "picard ScatterIntervalsByNs REFERENCE={input.ref} OUTPUT={output.intervals} MAX_TO_MERGE={params.minNmer} OUTPUT_TYPE=ACGT &> {log}\n"

rule format_interval_list:
    input:
        intervals = "results/{ref_name}/intervals/picard_interval_list.list"
    output:
        intervals = "results/{ref_name}/intervals/master_interval_list.list"
    run:
        with open(output.intervals, "w") as out:
            with open(input.intervals, "r") as inp:
                for line in inp:
                    if not line.startswith("@"):
                        line = line.strip().split("\t")
                        chrom, start, end = line[0], line[1], line[2]
                        print(f"{chrom}:{start}-{end}", file=out)
    

checkpoint create_db_intervals:
    input:
        ref = "config/{ref_name}.fasta",
        fai = "config/{ref_name}.fasta.fai",
        dictf = "config/{ref_name}.dict",
        intervals = "results/intervals/master_interval_list.list"
    output:
        fof = "results/{ref_name}/intervals/db_intervals/intervals.txt",
        out_dir = directory("results/{ref_name}/intervals/db_intervals"),
    params:
        max_intervals = get_db_interval_count
    log:
        "logs/{ref_name}/db_intervals/log.txt"
    benchmark:
        "benchmarks/{ref_name}/db_intervals/benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['create_db_intervals']['threads'] * 4000
    conda:
        '../envs/bam2vcf.yml'
    shell:
        """
        gatk SplitIntervals -L {input.intervals} \
        -O {output.out_dir} -R {input.ref} -scatter {params} \
        -mode INTERVAL_SUBDIVISION \
        --interval-merging-rule OVERLAPPING_ONLY &> {log}
        ls -l {output.out_dir}/*scattered.interval_list > {output.fof}
        """

checkpoint create_gvcf_intervals:
    input:
        ref = "config/{ref_name}.fasta",
        fai = "config/{ref_name}.fasta.fai",
        dictf = "config/{ref_name}.dict",
        intervals = "results/{ref_name}/intervals/master_interval_list.list"
    output:
        fof = "results/{ref_name}/intervals/gvcf_intervals/intervals.txt",
        out_dir = directory("results/{ref_name}/intervals/gvcf_intervals"),
    params:
        max_intervals = config["num_gvcf_intervals"]
    log:
        "logs/{ref_name}/gvcf_intervals/log.txt"
    benchmark:
        "benchmarks/{ref_name}/gvcf_intervals/benchmark.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['create_gvcf_intervals']['threads'] * 4000
    conda:
        '../envs/bam2vcf.yml'
    shell:
        """
        gatk SplitIntervals -L {input.intervals} \
        -O {output.out_dir} -R {input.ref} -scatter {params} \
        -mode BALANCING_WITHOUT_INTERVAL_SUBDIVISION \
        --interval-merging-rule OVERLAPPING_ONLY  &> {log}
        ls -l {output.out_dir}/*scattered.interval_list > {output.fof}
        """