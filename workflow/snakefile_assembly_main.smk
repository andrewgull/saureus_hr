from snakemake.io import touch, directory, temp, expand

# command to run the pipeline on 14 threads and 10Gb of RAM:
# snakemake --use-conda --cores 14 --resources mem_mb=10000
# snakemake --dag results/final/DA63360_all.done | dot -Tpng > dag_new.png

rule all:
    input:
        expand("results/final/{strain}_all.done", strain=config['strains'])

rule filter_nanopore:
    input: "resources/data_raw/{strain}/Nanopore/{strain}_all.fastq.gz"
    output: "results/data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    message: "executing filtlong on {wildcards.strain} long reads"
    log: "results/logs/{strain}_filtlong.log"
    conda: "filtlong-env"
    threads: 14
    params: min_len=config["min_nanopore_length"]
    shell: "filtlong --min_length {params.min_len} {input} 2> {log} | pigz -c -p {threads} > {output}"

# Make an assembly with Unicycler or FLye-Medaka-Polypolish
rule adaptive_hybrid_assembly:
    input:
        short_reads_1 = "results/data_filtered/{strain}/Illumina/{strain}_1.fq.gz",
        short_reads_2 = "results/data_filtered/{strain}/Illumina/{strain}_2.fq.gz",
        long_reads = "results/data_filtered/{strain}/Nanopore/{strain}_all.fastq.gz"
    output:
        assembly_dir = directory("results/assemblies/{strain}"),
        draft_dir = directory("results/drafts/{strain}"),
        polish_dir = directory("results/polished/{strain}")
    threads: 18
    message:
        "executing assembly script with {threads} threads on {wildcards.strain} reads"
    log:
        "results/logs/{strain}_assembly.log"
    conda: "uni+fmp-env"
    params: basecaller=config["basecaller"], genome_size=config["genome_size"], coverage=config["coverage"], genome_length=config["genome_length"], cov_threshold=config["cov_threshold"]
    script:
        "scripts/adaptive_hybrid_assembly.py"

rule QC_assembly:
    input: "results/assemblies/{strain}", "resources/busco_downloads"
    output:
        directory("results/qualcheck_assembly/{strain}")
    threads: 18
    message: "executing BUSCO and QUAST with {threads} threads on {wildcards.strain} assembly"
    conda: "busco-quast-env"
    params: tax_dataset=config["tax_dataset"]
    script:
        "scripts/QC_assembly.py"

rule bwa_map:
    input:
        assembly_dir = "results/assemblies/{strain}",
        short_read_1= "resources/parents/{strain}/Illumina/{strain}_1.fq.gz",
        short_read_2 = "resources/parents/{strain}/Illumina/{strain}_2.fq.gz"
    output:
        # an output file marked as temp is deleted after all rules that use it as an input are completed
        temp("results/mapping/{strain}/assembly.sam")
    threads: 18
    message: "executing BWA with {threads} threads on {wildcards.strain} assembly"
    log: mem = "results/logs/{strain}_bwa_mem.log",
        index = "results/logs/{strain}_bwa_index.log"
    conda: "bwa-env"
    shell:
        "bwa index {input.assembly_dir}/assembly.fasta &> {log.index} && "
        "bwa mem -t {threads} {input.assembly_dir}/assembly.fasta {input.short_read_1} {input.short_read_2} -o {output} &> {log.mem}"

rule samtools:
    input: "results/mapping/{strain}/assembly.sam"
    output: "results/mapping/{strain}/assembly.bam"
    threads: 18
    message: "executing SAMTOOLS: VIEW-SORT with {threads} threads on {wildcards.strain} mapping file"
    log: "results/logs/{strain}_samtools.log"
    conda: "samtools-env"
    shell:
        "samtools view -b {input} | samtools sort -o {output} -O BAM -@ {threads} &> {log} && samtools index {output}"

rule unmapped:
    input: "results/mapping/{strain}/assembly.bam"
    output:
        r1 = "results/mapping/{strain}/unmapped_1.fastq",
        r2 = "results/mapping/{strain}/unmapped_2.fastq"
    threads: 18
    message: "executing SAMTOOLS: VIEW-FASTQ with {threads} threads on {wildcards.strain} BAM file"
    log: "results/logs/{strain}_unmapped.log"
    conda: "samtools-env"
    shell:
        "samtools view -@ {threads} -u -f 12 -F 256 {input} | samtools fastq -1 {output.r1} -2 {output.r2} -@ {threads} &> {log}"

rule plasmid_assembly:
    input:
        r1 = "results/mapping/{strain}/unmapped_1.fastq",
        r2 = "results/mapping/{strain}/unmapped_2.fastq"
    output:
        directory("results/plasmids/{strain}")
    threads: 18
    message: "executing SPAdes in plasmid mode with {threads} threads on unmapped reads of {wildcards.strain}"
    log: "results/logs/{strain}_spades.log"
    conda: "spades-env"
    shell:
        # || true prevents the rule from failing when spades throws an error; this happens when unmapped files are too small
        "spades.py --plasmid -1 {input.r1} -2 {input.r2} -t {threads} -o {output} &> {log} || true"

rule assembly_summary:
    input: "results/assemblies/{strain}", "results/plasmids/{strain}"
    output: "results/assemblies_joined/{strain}/summary.tsv"
    threads: 1
    message: "summarizing unicycler and SPAdes assemblies of strain {wildcards.strain}"
    log: "results/logs/{strain}_assembly_summary.log"
    script: 
        "scripts/assembly_summary.py"

rule join_assemblies:
    input: "results/assemblies/{strain}", "results/plasmids/{strain}"
    output:
        "results/assemblies_joined/{strain}/assembly.fasta"
    threads: 1
    message: "joining Unicycler assembly and SPAdes plasmid assembly together, strain {wildcards.strain}"
    log: "results/logs/{strain}_joiner.log"
    script:
        "scripts/join_two_fastas.py"

rule genome_annotation:
    input: "results/assemblies_joined/{strain}/assembly.fasta"
    output:
        directory("results/annotations/{strain}/prokka")
    threads: 18
    message: "executing PROKKA with {threads} threads on full assembly of {wildcards.strain}"
    log: "results/logs/{strain}_prokka.log"
    conda: "prokka-env"
    params: centre=config["centre"], minlen=config["minlen"], 
            genus=config["genus"], species=config["species"]
    shell:
        # skip tRNAs search?
        "prokka --addgenes --addmrna --compliant --notrna --outdir {output} --prefix {wildcards.strain}_genomic --centre {params.centre} --genus {params.genus} "
        "--species {params.species} --strain {wildcards.strain} --kingdom Bacteria --cpus {threads} "
        "--mincontiglen {params.minlen} {input} &> {log}"

rule rename_gbk:
    input: "results/annotations/{strain}/prokka"
    output: "results/annotations/{strain}/prokka_renamed/{strain}_genomic.gbk"
    params:
        filename="{strain}_genomic.gbk"
    script:
        "scripts/rename_genomic_gbk.py"

rule trna_annotation:
    # tRNA genes only - tRNAScan-SE
    input: "results/assemblies_joined/{strain}/assembly.fasta"
    output:
        general = "results/annotations/{strain}/trna/trna_gen.txt",
        struct = "results/annotations/{strain}/trna/trna_struct.txt",
        iso = "results/annotations/{strain}/trna/trna_iso.txt",
        stats = "results/annotations/{strain}/trna/trna_stat.txt",
        bed = "results/annotations/{strain}/trna/trna_coords.bed",
        gff = "results/annotations/{strain}/trna/trna_feat.gff",
        fasta = "results/annotations/{strain}/trna/trna_seq.fasta"
    threads: 18
    message: "executing tRNAScan-SE with {threads} threads on full assembly of {wildcards.strain} strain"
    log: "results/logs/{strain}_trnascan.log"
    conda: "trnascan-env"
    shell:
        "tRNAscan-SE -B --forceow -o {output.general} -f {output.struct} -s {output.iso} -m {output.stats} -b {output.bed} "
        "-j {output.gff} -a {output.fasta} -l {log} --thread {threads} {input} &> {log}"

rule join_annotations:
    input:
        prokka="results/annotations/{strain}/prokka",
        trnascan="results/annotations/{strain}/trna/trna_seq.fasta"
        # mode="annotation"
    output:
        "results/annotations/{strain}/joined/annotation.fasta"  # it's not supposed to be used as input for RGI tool
    script:
        "scripts/join_two_fastas.py"

rule resistance_genes:
    input:
        "results/annotations/{strain}/prokka"
    output:
        "results/resistance_genes/{strain}/rgi_table.txt" # IT'S JUST A PREFIX!
    threads: 18
    message: "executing RGI with {threads} threads on predicted genes/proteins from {wildcards.strain}"
    log: "results/logs/{strain}_rgi.log"
    conda: "rgi-env"
    shell:
        "output=$(echo '{output}' | cut -d'.' -f 1) && "
        "rgi main --input_sequence {input}/{wildcards.strain}_genomic.faa --output_file $output  "
        "--input_type protein --local  --num_threads {threads} --include_loose --clean &> {log}"

rule final:
    input:
        qc_ass="results/qualcheck_assembly/{strain}",
        qc_ill_raw="results/qualcheck_reads/{strain}/Illumina/{strain}_summary.tsv",
        qc_nan_raw="results/qualcheck_reads/{strain}/Nanopore/{strain}_summary.tsv",
        rgi = "results/resistance_genes/{strain}/rgi_table.txt",
        prokka="results/annotations/{strain}/prokka",
        trnascan="results/annotations/{strain}/trna/trna_gen.txt",
        summary="results/assemblies_joined/{strain}/summary.tsv",
        renamed_gbk="results/annotations/{strain}/prokka_renamed/{strain}_genomic.gbk"
    output: touch("results/final/{strain}_all.done")
    shell: "echo 'DONE'"

onsuccess:
    print("Workflow finished, no errors")

