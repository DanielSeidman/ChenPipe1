rule trim_galore_call:
    input:
        r1="results/filtered_fastqs/{sample}/{run}_1.fastq.gz",
        r2="results/filtered_fastqs/{sample}/{run}_2.fastq.gz",
    output:
        r1trim="results/filtered_fastqs/{sample}/{run}_trimgalore/{run}_1.fastq.gz",
        r2trim="results/filtered_fastqs/{sample}/{run}_trimgalore/{run}_2.fastq.gz",
    conda:
        "../envs/fastq2bam.yml"
    threads:
        resources['trim_galore_call']['threads']
    resources:
        mem_mb = lambda wildcards,attempt: attempt * resources['trim_galore_call']['mem']

    shell:
        "trim_galore --paired {input.r1} {input.r2} -o {run}_trimgalore/ "


rule bwa_map:
    params:
        rg=get_read_group,
        ref=config['reference']
    input:
        r1 = "results/filtered_fastqs/{sample}/{run}_trimgalore/{run}_1.fastq.gz",
        r2 = "results/filtered_fastqs/{sample}/{run}_trimgalore/{run}_2.fastq.gz",
        indexes = expand("{params.ref}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
    output: 
        bam = temp("results/bams/preMerge/{sample}/{run}.bam"),
        bai = temp("results/bams/preMerge/{sample}/{run}.bam.bai"),

    conda:
        "../envs/fastq2bam.yml"
    threads:
        resources['bwa_map']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['bwa_map']['mem']
    log:
        "logs/bwa_mem/{sample}/{run}.txt"
    benchmark:
        "benchmarks/bwa_mem/{sample}_{run}.txt"
    shell:
        "bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -o {output.bam} - && samtools index {output.bam} {output.bai}"

rule merge_bams:
    input:
        merge_bams_input
    output:
        bam = temp("results/bams/postMerge/{sample}.bam"),
        bai = temp("results/bams/postMerge/{sample}.bam.bai")
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/merge_bams/{sample}.txt"
    benchmark:
        "benchmarks/merge_bams/{sample}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['merge_bams']['mem']
    shell:
        "samtools merge {output.bam} {input} && samtools index {output.bam} > {log}"

rule dedup:
    input:
        unpack(dedup_input)
    output:
        dedupBam = "results/bams/{sample}_final.bam",
        dedupBai = "results/bams/{sample}_final.bam.bai",
    conda:
        "../envs/sambamba.yml"
    resources:
        threads = resources['dedup']['threads'],
        mem_mb = lambda wildcards, attempt: attempt * resources['dedup']['mem']
    log:
        "logs/sambamba_dedup/{sample}.txt"
    benchmark:
        "benchmarks/sambamba_dedup/{sample}.txt"
    shell:
        "sambamba markdup -t {threads} {input.bam} {output.dedupBam} 2> {log}"