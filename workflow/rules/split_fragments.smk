rule total_fragments:
	input:
		"results/aligned_reads/filtered/{sample}.bam",
		"results/aligned_reads/filtered/{sample}.bam.bai"
	output:
		"results/aligned_reads/split_fragments/{sample}_total.bam",
		"results/aligned_reads/split_fragments/{sample}_total.bam.bai"
	log:
		"logs/split_fragments/{sample}_total.log"
	shell:
		"""
		mv {input[0]} {output[0]} 2> {log}
		mv {input[1]} {output[1]} 2>> {log}
		"""
		
rule small_fragments:
	input:
		"results/aligned_reads/split_fragments/{sample}_total.bam",
		"results/aligned_reads/split_fragments/{sample}_total.bam.bai"
	output:
		"results/aligned_reads/split_fragments/{sample}_small.bam"
	log:
		"logs/split_fragments/{sample}_small.log"
	params:
		extra=" -f bam -F 'template_length < 100 and template_length > -100'"
	threads: 8
	wrapper:
		"v1.1.0/bio/sambamba/view"
		
rule large_fragments:
	input:
		"results/aligned_reads/split_fragments/{sample}_total.bam",
		"results/aligned_reads/split_fragments/{sample}_total.bam.bai"
	output:
		"results/aligned_reads/split_fragments/{sample}_large.bam"
	log:
		"logs/split_fragments/{sample}_large.log"
	params:
		extra=" -f bam -F 'not (template_length < 180 and template_length > -180)'"
	threads: 8
	wrapper:
		"v1.1.0/bio/sambamba/view"

rule samtools_index_small_fragments:
    input:
        "results/aligned_reads/split_fragments/{sample}_small.bam"
    output:
        "results/aligned_reads/split_fragments/{sample}_small.bam.bai"
    log:
        "logs/samtools_index/{sample}_small.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"

rule samtools_index_large_fragments:
    input:
        "results/aligned_reads/split_fragments/{sample}_large.bam"
    output:
        "results/aligned_reads/split_fragments/{sample}_large.bam.bai"
    log:
        "logs/samtools_index/{sample}_large.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"