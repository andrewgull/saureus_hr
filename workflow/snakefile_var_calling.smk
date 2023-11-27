# this pipeline maps mutant reads onto a corresponding parent's genome
from snakemake.io import touch, directory, expand


#configfile: "config_var_calling.yaml"
RUNS = ["1", "2"]



rule all:
    input:
        expand("results/final/{parent}_{mutant}_all.done", parent=config["parent"], mutant=config["mutant"])

rule QC:
    input:
        reads = expand("resources/mutants/{parent}/{mutant}/{mutant}_{run}.fq.gz", parent=config["parent"],
            mutant=config["mutant"], run=RUNS)
    output:
        directory("results/fastqc/{parent}_mutants")
    threads: 12
    message: "executing FastQC with {threads} threads on {wildcards.parent} mutant reads"
    log: fastqc="results/logs/{parent}_fastqc.log",
         multiqc = "results/logs/{parent}_multiqc.log"
    conda: "envs/qc.yaml"
    shell:
        "mkdir {output} && "
        "fastqc -t {threads} -o {output} {input} &> {log.fastqc} && "
        "multiqc -q {output} &> {log.multiqc}"

rule call_variants:
     input:
         reference="results/assemblies_joined/{parent}/assembly.fasta",
         R1="resources/mutants/{parent}/{mutant}/{mutant}_1.fq.gz",
         R2="resources/mutants/{parent}/{mutant}/{mutant}_2.fq.gz"
     output:
         directory("results/mutants_variants/{parent}/{mutant}")
     threads: 12
     #message: "executing SNIPPY with {threads} threads on {wildcards.mutant} read files"
     log: "results/logs/{parent}_{mutant}_snippy.log"
     conda: "envs/snippy.yaml"
     shell:
         "snippy --cpus {threads} --outdir {output} --ref {input.reference} --R1 {input.R1} --R2 {input.R2} &> {log}"

rule final:
    input:
      qc="results/fastqc/{parent}_mutants",
      vars="results/mutants_variants/{parent}/{mutant}"
    output: touch("results/final/{parent}_{mutant}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")
