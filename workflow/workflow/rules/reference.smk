rule index_reference:
    input:
        ref = "config/{ref_name}.fasta"
    output:
        indexes = expand("config/{{ref_name}}.fasta.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"]),
        fai = "config/{ref_name}.fasta.fai",
        dictf = "config/{ref_name}.dict"
    conda:
        "../envs/fastq2bam.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt * resources['index_ref']['threads'] * 4000
    log:
        "logs/{ref_name}/index_ref/log.txt"
    benchmark:
        "benchmarks/{ref_name}/index_ref/benchmark.txt"
    shell:
        """
        bwa index {input.ref} 2> {log}
        samtools faidx {input.ref} --output {output.fai} >> {log}
        samtools dict {input.ref} -o {output.dictf} >> {log} 2>&1
        """
