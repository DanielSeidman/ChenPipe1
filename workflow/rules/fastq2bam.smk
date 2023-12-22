rule trim_galore_call:
    input:
        unpack(get_reads)
    output:
        r1trim="results/{ref_name}/trimgalore_fastqs/{sample}/{run}_1_trimmed.fastq.gz",
        r2trim="results/{ref_name}/trimgalore_fastqs/{sample}/{run}_2_trimmed.fastq.gz",
    conda:
        "../envs/fastq2bam.yml"
    threads:
        resources['trim_galore_call']['threads']
    resources:
        mem_mb = lambda wildcards,attempt: attempt * resources['trim_galore_call']['mem']
    log:
        "logs/{ref_name}/trim_galore/{sample}/{run}.txt"
    shell:
        "trim_galore --paired {input.r1} {input.r2} -o results/{wildcards.ref_name}/trimgalore_fastqs/{wildcards.sample}/ "


rule bwa_map:
    input:
        r1 = "results/{ref_name}/trimgalore_fastqs/{sample}/{run}_1_trimmed.fastq.gz",
        r2 = "results/{ref_name}/trimgalore_fastqs/{sample}/{run}_2_trimmed.fastq.gz",
        ref = "config/{ref_name}.fasta",
        indexes=expand("config/{{ref_name}}.fasta.{ext}",ext=["sa", "pac", "bwt", "ann", "amb", "fai"]),
        dictf="config/{ref_name}.dict",
    output: 
        bam = temp("results/{ref_name}/bams/preMerge/{sample}/{run}.bam"),
        bai = temp("results/{ref_name}/bams/preMerge/{sample}/{run}.bam.bai"),
    params:
        rg=get_read_group,
    conda:
        "../envs/fastq2bam.yml"
    threads:
        resources['bwa_map']['threads']
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['bwa_map']['mem']
    log:
        "logs/{ref_name}/bwa_mem/{sample}/{run}.txt"
    benchmark:
        "benchmarks/{ref_name}/bwa_mem/{sample}_{run}.txt"
    shell:
        "bwa mem -M -t {threads} -R {params.rg} {input.ref} {input.r1} {input.r2} 2> {log} | samtools sort -o {output.bam} - && samtools index {output.bam} {output.bai}"

rule merge_bams:
    input:
        merge_bams_input
    output:
        bam = temp("results/{ref_name}/bams/postMerge/{sample}.bam"),
        bai = temp("results/{ref_name}/bams/postMerge/{sample}.bam.bai")
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{ref_name}/merge_bams/{sample}.txt"
    benchmark:
        "benchmarks/{ref_name}/merge_bams/{sample}.txt"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['merge_bams']['mem']
    shell:
        "samtools merge {output.bam} {input} && samtools index {output.bam} > {log}"

rule dedup:
    input:
        unpack(dedup_input)
    output:
        dedupBam = "results/{ref_name}/bams/{sample}_final.bam",
        dedupBai = "results/{ref_name}/bams/{sample}_final.bam.bai",
    conda:
        "../envs/sambamba.yml"
    resources:
        threads = resources['dedup']['threads'],
        mem_mb = lambda wildcards, attempt: attempt * resources['dedup']['mem']
    log:
        "logs/{ref_name}/sambamba_dedup/{sample}.txt"
    benchmark:
        "benchmarks/{ref_name}/sambamba_dedup/{sample}.txt"
    shell:
        "sambamba markdup -t {threads} {input.bam} {output.dedupBam} 2> {log}"