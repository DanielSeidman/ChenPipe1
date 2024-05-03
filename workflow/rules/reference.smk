ruleorder: download_reference > index_reference
# localrules: copy_reference, download_reference

# This does not work with SLURM as of 4/3/24. See here for more info:https://github.com/snakemake/snakemake-executor-plugin-slurm/issues/60
# rule copy_reference:
#     """Copies user-specified reference genome path to results dir to maintain ref_name wildcard"""
#     input:
#         ref = get_ref
#     output:
#         ref = "results/{ref_name}/data/genome/{ref_name}.fna"
#     log:
#         "logs/{ref_name}/copy_ref/log.txt"
#     shell:
#         #probably don't need to unzip but might as well.
#         """
#         gunzip -c {input.ref} 2> {log} > {output.ref} || cp {input.ref} {output.ref} &> {log}
#         """

rule download_reference:
    input:
        ref = get_ref
    output:
        ref = "results/{ref_name}/data/genome/{ref_name}.fna"
    params:
        dataset = "results/{ref_name}/data/genome/{ref_name}_dataset.zip",
        outdir = "results/{ref_name}/data/genome/{ref_name}"
    conda:
        "../envs/fastq2bam.yml"
    log:
        "logs/{ref_name}/download_ref/log.txt"
    benchmark:
        "benchmarks/{ref_name}/download_ref/benchmark.txt"
    shell:
        """
        if [ -z "{input.ref}" ]  # check if this is empty
        then
            mkdir -p {params.outdir}
            datasets download genome accession --exclude-gff3 --exclude-protein --exclude-rna --filename {params.dataset} {wildcards.ref_name} \
            && 7z x {params.dataset} -aoa -o{params.outdir} \
            && cat {params.outdir}/ncbi_dataset/data/{wildcards.ref_name}/*.fna > {output.ref}
        else
            gunzip -c {input.ref} 2> {log} > {output.ref} || cp {input.ref} {output.ref} &> {log}
        fi
        """

rule index_reference:
    input:
        ref = "results/{ref_name}/data/genome/{ref_name}.fna"
    output:
        indexes = expand("results/{{ref_name}}/data/genome/{{ref_name}}.fna.{ext}", ext=["sa", "pac", "bwt", "ann", "amb"]),
        fai = "results/{ref_name}/data/genome/{ref_name}.fna.fai",
        dictf = "results/{ref_name}/data/genome/{ref_name}.dict"
    conda:
        "../envs/fastq2bam.yml"
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
