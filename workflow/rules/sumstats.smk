rule bam_sumstats:
    input:
        unpack(get_bams),
        ref = "results/{ref_name}/data/genome/{ref_name}.fna",
    output:
        cov = "results/{ref_name}/summary_stats/{sample}_coverage.txt",
        alnSum = "results/{ref_name}/summary_stats/{sample}_AlnSumMets.txt",
    conda:
        "../envs/fastq2bam.yml"
    shell:
        """
        samtools coverage --output {output.cov} {input.bam}
        samtools flagstat -O tsv {input.bam} > {output.alnSum}
        """

rule sentieon_bam_stats:
    input:
        unpack(get_bams),
        indexes = expand("results/{{ref_name}}/data/genome/{{ref_name}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
        ref = "results/{ref_name}/data/genome/{ref_name}.fna"
    params:
        lic = config['sentieon_lic']
    output:
        insert_file = "results/{ref_name}/summary_stats/{sample}_insert_metrics.txt",
        qd = "results/{ref_name}/summary_stats/{sample}_qd_metrics.txt",
        gc = "results/{ref_name}/summary_stats/{sample}_gc_metrics.txt",
        gc_summary = "results/{ref_name}/summary_stats/{sample}_gc_summary.txt",
        mq = "results/{ref_name}/summary_stats/{sample}_mq_metrics.txt"
    conda:
        "../envs/sentieon.yml"
    shell:
        """
        export SENTIEON_LICENSE={params.lic}
        sentieon driver -r {input.ref} \
        -t {threads} -i {input.bam} \
        --algo MeanQualityByCycle {output.mq} \
        --algo QualDistribution {output.qd} \
        --algo GCBias --summary {output.gc_summary} {output.gc} \
        --algo InsertSizeMetricAlgo {output.insert_file}
        """

rule collect_fastp_stats:
    input:
        collect_fastp_stats_input
    output:
        "results/{ref_name}/summary_stats/{sample}_fastp.out"
    shell:
        "cat {input} > {output}"

rule collect_sumstats:
    input:
        unpack(get_input_sumstats)
    output:
        "results/{ref_name}/summary_stats/{prefix}_bam_sumstats.txt"
    run:
        if not config['sentieon']:
            FractionReadsPassFilter, NumReadsPassFilter = collectFastpOutput(input.fastpFiles)
            aln_metrics = collectAlnSumMets(input.alnSumMetsFiles)
            SeqDepths, CoveredBases = collectCoverageMetrics(input.coverageFiles)
            printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumReadsPassFilter, output[0])
        else:
            FractionReadsPassFilter, NumReadsPassFilter = collectFastpOutput(input.fastpFiles)
            aln_metrics = collectAlnSumMets(input.alnSumMetsFiles)
            SeqDepths, CoveredBases = collectCoverageMetrics(input.coverageFiles)
            median_inserts, median_insert_std = collect_inserts(input.insert_files)
            printBamSumStats(SeqDepths, CoveredBases, aln_metrics, FractionReadsPassFilter, NumReadsPassFilter, output[0], median_inserts, median_insert_std)