wildcard_constraints:
        frag_size="small|large|total"
rule make_bigwigs_ind:
	input:
		bam = "results/aligned_reads/split_fragments/{sample}_{frag_size}.bam",
		bai = "results/aligned_reads/split_fragments/{sample}_{frag_size}.bam.bai"
	output:
		temp("results/bigwigs/coverage/individual/{sample}_{frag_size}.bw")
	conda:
		"../envs/deeptools.yaml"
	params:
		extra=config["params"]["bigwigs_ind"] 
	threads: 8
	shell:
		"bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"

rule merge_bam:
	input:
		get_bam_merge
	output:
		temp("results/aligned_reads/merged/{sample_group}_{frag_size}.bam")
	params:
		"" # optional additional parameters as string
	threads:  # Samtools takes additional threads through its option -@
		8     # This value - 1 will be sent to -@
	wrapper:
		"v1.1.0/bio/samtools/merge"		

rule samtools_index_merged:
    input:
        "results/aligned_reads/merged/{sample}_{frag_size}.bam"
    output:
        temp("results/aligned_reads/merged/{sample}_{frag_size}.bam.bai")
    log:
        "logs/samtools_index/{sample}_{frag_size}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "v1.1.0/bio/samtools/index"
        
rule make_bigwigs_merged:
	input:
		bam = "results/aligned_reads/merged/{sample}_{frag_size}.bam",
		bai = "results/aligned_reads/merged/{sample}_{frag_size}.bam.bai"
	output:
		temp("results/bigwigs/coverage/merged/{sample}_{frag_size}.bw")
	conda:
		"../envs/deeptools.yaml"
	params:
		extra=config["params"]["bigwigs_merged"] 
	threads: 8
	shell:
		"bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"


rule zscore_normalize_ind_bigwigs:
	input:
		"results/bigwigs/coverage/individual/{sample}_{frag_size}.bw"
	output:
		"results/bigwigs/zscore_normalized/individual/{sample}_{frag_size}.bw"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	script:
		"../scripts/zscore_normalize_bw.R"

rule zscore_normalize_merged_bigwigs:
	input:
		"results/bigwigs/coverage/merged/{sample}_{frag_size}.bw"
	output:
		"results/bigwigs/zscore_normalized/merged/{sample}_{frag_size}.bw"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	script:
		"../scripts/zscore_normalize_bw.R"