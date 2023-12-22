ruleorder: download_reference > index_reference

rule download_reference:
    input:
        ref = get_ref
    output:
        ref = "config/{ref_name}.fasta"
    params:
        dataset = "config/{ref_name}_dataset.zip",
        outdir = "config/{ref_name}"
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/download_ref/log.txt"
    benchmark:
        "benchmarks/download_ref/benchmark.txt"
    shell:
        """
        if [ -z "{input.ref}" ]  # check if this is empty
        then
            mkdir -p {params.outdir}
            datasets download genome accession --exclude-gff3 --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.refGenome} \
            && 7z x {params.dataset} -aoa -o{params.outdir} \
            && cat {params.outdir}/ncbi_dataset/data/{wildcards.refGenome}/*.fasta > {output.ref}
        else
            cp {input.ref} {output.ref}
        fi
        """
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
        mem_mb = lambda wildcards, attempt: attempt * resources['index_ref']['mem']
    log:
        "logs/index_ref/log.txt"
    benchmark:
        "benchmarks/index_ref/benchmark.txt"
    shell:
        """
        bwa index {input.ref} 2> {log}
        samtools faidx {input.ref} --output {output.fai} >> {log}
        samtools dict {input.ref} -o {output.dictf} >> {log} 2>&1
        """
