if config["filter_chroms"]:
	rule filter_multireads:
		input:
			"results/aligned_reads/mapped/{sample}.bam"
		output:
			temp("results/aligned_reads/unireads/{sample}.bam")
		log:
			"logs/filter_multireads/{sample}.log"
		params:
			extra="-bh -q 30" # optional params string
		wrapper:
			"0.77.0/bio/samtools/view"

	rule filter_chroms:
		input:
			bam="results/aligned_reads/unireads/{sample}.bam",
			keep_chroms=rules.define_keep_chroms.output
		output:
			temp("results/aligned_reads/filtered/{sample}.bam")
		log:
			"logs/filter_chroms/{sample}.log"
		params:
			extra="-bh -L {}".format(input.keep_chroms) # optional params string
		wrapper:
			"0.77.0/bio/samtools/view"
else:
	rule filter_multireads:
		input:
			"results/aligned_reads/mapped/{sample}.bam"
		output:
			temp("results/aligned_reads/filtered/{sample}.bam")
		log:
			"logs/filter_multireads/{sample}.log"
		params:
			extra="-bh -q 30" # optional params string
		wrapper:
			"0.77.0/bio/samtools/view"
		
rule samtools_sort:
    input:
       "results/aligned_reads/filtered/{sample}.bam"
    output:
        "results/aligned_reads/sorted/{sample}.bam"
    log:
        "logs/samtools_sort/{sample}.log"
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.77.0/bio/samtools/sort"

rule samtools_index:
    input:
        "results/aligned_reads/sorted/{sample}.bam"
    output:
        "results/aligned_reads/sorted/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"

# rule align_stats:
# 
rule samtools_idxstats:
    input:
        bam="results/aligned_reads/sorted/{sample}.bam",
        idx="results/aligned_reads/sorted/{sample}.bam.bai"
    output:
        "results/aligned_reads/stats/{sample}.idxstats"
    log:
        "logs/samtools/idxstats/{sample}.log"
    wrapper:
        "0.77.0/bio/samtools/idxstats" 
# rule summarize_read_processing:
