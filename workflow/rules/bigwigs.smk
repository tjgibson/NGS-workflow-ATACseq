rule make_bigwigs_ind:
	input:
		bam = "results/aligned_reads/sorted/{sample}.bam",
		bai = "results/aligned_reads/sorted/{sample}.bam.bai"
	output:
		"results/bigwigs/coverage/individual/{sample}.bw"
	conda:
		"../envs/deeptools.yaml"
	params:
		extra="--binSize 10"
	threads: 8
	shell:
		"bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"

rule merge_bam:
	input:
		get_bam_merge
	output:
		temp("results/aligned_reads/merged/{sample_group}.bam")
	params:
		"" # optional additional parameters as string
	threads:  # Samtools takes additional threads through its option -@
		8     # This value - 1 will be sent to -@
	wrapper:
		"0.77.0/bio/samtools/merge"		

rule samtools_index_merged:
    input:
        "results/aligned_reads/merged/{sample}.bam"
    output:
        temp("results/aligned_reads/merged/{sample}.bam.bai")
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"
        
rule make_bigwigs_merged:
	input:
		bam = "results/aligned_reads/merged/{sample}.bam",
		bai = "results/aligned_readsmerged/{sample}.bam.bai"
	output:
		"results/bigwigs/coverage/merged/{sample}.bw"
	conda:
		"../envs/deeptools.yaml"
	params:
		extra="--binSize 10"
	threads: 8
	shell:
		"bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"
		
