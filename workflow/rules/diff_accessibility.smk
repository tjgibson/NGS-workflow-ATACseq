rule feature_counts:
	input:
		samples=get_featurecounts_input,  # list of sam or bam files
		annotation="results/peaks/final/{experiment}.saf",

	output:
		multiext(
			"results/count_tables/{experiment}",
			".featureCounts",
			".featureCounts.summary",
		),
	threads: 2
	params:
		r_path="",  # implicitly sets the --Rpath flag
		extra="-O -f -F 'SAF' -p",
	log:
		"logs/feature_counts/{experiment}.log",
	wrapper:
		"v1.3.2/bio/subread/featurecounts"

if config["run_diff_accessibility"]:
	rule DEseq2:
		input:
			"results/count_tables/{experiment}.featureCounts"
		output:
			"results/DEseq2/{experiment}.dds"
		params:
			samples=config["samples"],
			model=config["diff_accessibility"]["model"],
			count_threshold=config["diff_accessibility"]["count_threshold"],
		conda:
			"../envs/DEseq2.yaml",
		log:
			"logs/DEseq2/{experiment}.log",
		script:
			"../scripts/DEseq2.R"

	rule DEseq2_results:
		input:
			dds="results/DEseq2/{experiment}.dds",
			peaks="results/peaks/final/{experiment}.narrowPeak"
		output:
			"results/DEseq2/{experiment}_{contrast}_results.tsv"
		params:
			 contrast=get_contrast,
			 padj_cutoff=config["diff_accessibility"]["padj_cutoff"],
			 FC_cutoff=config["diff_accessibility"]["log2FC_cutoff"],
		conda:
			"../envs/DEseq2.yaml",
		log:
			"logs/DEseq2_results/{experiment}_{contrast}.log",
		script:
			"../scripts/DEseq2_results.R"