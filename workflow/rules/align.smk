rule get_sra_se:
	output:
		temp("data/sra/se/{accession}.fastq.gz"),
	log:
		"logs/get_sra/{accession}.log",
	wrapper:
		"0.77.0/bio/sra-tools/fasterq-dump"
                
rule get_sra_pe:
	output:
		temp("data/sra/pe/{accession}_1.fastq.gz"),
		temp("data/sra/pe/{accession}_2.fastq.gz"),
	log:
		"logs/get_sra/{accession}.log",
	wrapper: 
		"0.77.0/bio/sra-tools/fasterq-dump"
        
rule merge_fastqs:
    input:
        get_fq_merge,
    output:
        temp("data/merged/{sample}_{read}.fastq.gz"),
    log:
        "logs/merge-fastqs/{sample}_{read}.log",
    wildcard_constraints:
        read="single|1|2",
    shell:
        "zcat {input} > {output} 2> {log}"
        
rule bowtie2_align_se:
	input:
		sample=get_bowtie2_input,
		index=rules.bowtie2_index.output,
	output:
		temp("results/mapped/{sample}.bam")
	log:
		"logs/bowtie2/{sample}.log"
	params:
		index=lambda w, input: os.path.splitext(input.index[0])[0],  # prefix of reference genome index (built with bowtie2-build)
		extra=""  # optional parameters
	threads: 8  # Use at least two threads
	wrapper:
		"0.77.0/bio/bowtie2/align"
