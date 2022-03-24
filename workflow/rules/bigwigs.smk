rule make_bigwigs_ind:
	input:
		bam = "results/aligned_reads/filtered/{sample}.bam",
		bai = "results/aligned_reads/filtered/{sample}.bam.bai"
	output:
		temp("results/bigwigs/coverage/individual/{sample}.bw")
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
		temp("results/aligned_reads/merged/{sample_group}.bam")
	params:
		"" # optional additional parameters as string
	threads:  # Samtools takes additional threads through its option -@
		8     # This value - 1 will be sent to -@
	wrapper:
		"v1.1.0/bio/samtools/merge"		

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
        "v1.1.0/bio/samtools/index"
        
rule make_bigwigs_merged:
	input:
		bam = "results/aligned_reads/merged/{sample}.bam",
		bai = "results/aligned_reads/merged/{sample}.bam.bai"
	output:
		temp("results/bigwigs/coverage/merged/{sample}.bw")
	conda:
		"../envs/deeptools.yaml"
	params:
		extra=config["params"]["bigwigs_merged"] 
	threads: 8
	shell:
		"bamCoverage --bam {input.bam} -o {output} -p {threads} {params.extra}"


rule zscore_normalize_ind_bigwigs:
	input:
		"results/bigwigs/coverage/individual/{sample}.bw"
	output:
		"results/bigwigs/zscore_normalized/individual/{sample}.bw"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	script:
		"../scripts/zscore_normalize_bw.R"

rule zscore_normalize_merged_bigwigs:
	input:
		"results/bigwigs/coverage/merged/{sample}.bw"
	output:
		"results/bigwigs/zscore_normalized/merged/{sample}.bw"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	script:
		"../scripts/zscore_normalize_bw.R"


rule compute_scaling_factors:
	input:
		mapping_stats=get_scaling_input
	output:
		"results/scaling_factors/individual_scaling_factors.tsv",
		"results/scaling_factors/merged_scaling_factors.tsv"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	script:
		"../scripts/compute_scaling_factors.R"


rule spikeIn_normalize_ind_bigwigs:
	input:
		bw=get_ind_spikeIn_input,
		scaling_factors="results/scaling_factors/individual_scaling_factors.tsv"
	output:
		"results/bigwigs/spikeIn_normalized/individual/{sample}.bw"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	script:
		"../scripts/spikeIn_normalize_bw.R"

rule spikeIn_normalize_merged_bigwigs:
	input:
		bw=get_merged_spikeIn_input,
		scaling_factors="results/scaling_factors/merged_scaling_factors.tsv"
	output:
		"results/bigwigs/spikeIn_normalized/merged/{sample}.bw"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	script:
		"../scripts/spikeIn_normalize_bw.R"
