import pandas as pd



# read in table with sample metadata
samples = (
    pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

units = (
    pd.read_csv(config["units"], sep="\t", dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)

# function to check config files for inclusion of optional workflow steps
def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))

def get_fq_merge(wildcards):
	unit = units.loc[wildcards.sample]
	if all(pd.isna(unit["fq1"])):
		accession = unit["sra"]
		if all(unit["read_format"] == "SE"):
			return expand(
				"data/sra/se/{accession}.fastq.gz", accession=accession)
		else:
			return expand(
				"data/sra/pe/{accession}_{read}.fastq.gz", accession=accession, read = wildcards.read[-1])
	if all(unit["read_format"] == "SE"):
		return units.loc[wildcards.sample, "fq1"].tolist()
	fq = "fq{}".format(wildcards.read[-1])
	return units.loc[wildcards.sample, fq].tolist()


def get_NGmerge_input(wildcards):
	if not is_activated("mergeReads"):
		unit = units.loc[wildcards.sample]
		if all(pd.isna(unit["fq1"])):
			# SRA sample (always paired-end for now)
			accession = unit["sra"]
			if all(unit["read_format"] == "PE"):
				return expand("data/sra/pe/{accession}_{read}.fastq.gz", accession=accession, read=[1,2])
		fastqs = units.loc[(wildcards.sample, wildcards.sample), ["fq1", "fq2"]]
		if len(fastqs) == 2:
			return [fastqs.fq1, fastqs.fq2]
	unit = units.loc[wildcards.sample]
	return ["data/merged/{sample}_1.fastq.gz", "data/merged/{sample}_2.fastq.gz"]

def get_bam_merge(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	group = pd.unique(unit["sample_name"])
	return expand(
		"results/aligned_reads/split_fragments/{group}_{frag_size}.bam", group=group, frag_size = wildcards.frag_size)


def get_macs2_merged_input(wildcards):
	sample =  samples[samples["experiment"] == wildcards.experiment]
	in_samples = pd.unique(sample["sample_name"])
	return expand(
		"results/aligned_reads/split_fragments/{sample}_small.bam", sample=in_samples)

def get_macs2_merged_input_by_sample(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	sample_names = pd.unique(unit["sample_name"])
	return expand("results/aligned_reads/split_fragments/{sample}_small.bam",sample=sample_names)


def get_featurecounts_input(wildcards):
	sample =  samples[samples["experiment"] == wildcards.experiment]
	in_samples = pd.unique(sample["sample_name"])
	return expand(
		"results/aligned_reads/split_fragments/{sample}_small.bam", sample=in_samples)

def get_contrast(wildcards):
	return config["diff_accessibility"]["experiments"][wildcards.experiment]["contrasts"][wildcards.contrast]


def get_final_output():
	final_output = []	
	# z-score normalized bigwigs for individual replicates
	final_output.extend(expand(
					[
						"results/bigwigs/zscore_normalized/individual/{sample}_{frag_size}.bw"
					],
					sample = units["sample_name"], frag_size = ["small","large","total"]
				)
			)

	
	# z-score normalized bigwigs for merged replicates
	final_output.extend(expand(
					[
						"results/bigwigs/zscore_normalized/merged/{sample}_{frag_size}.bw"
					],
					sample = units["sample_group"], frag_size = ["small","large","total"]
				)
			)





	# add indvidual peak output
	final_output.extend(expand(
				[
					"results/peaks/individual/{sample}{ext}"
				],
				sample = pd.unique(units["sample_name"]),
				ext = ["_peaks.xls", "_peaks.narrowPeak","_summits.bed"]
			)
		)

	# add merged peak output
	final_output.extend(expand(
				[
					"results/peaks/merged_all/{experiment}{ext}"
				],
				experiment = pd.unique(samples["experiment"]),
				ext = ["_peaks.xls", "_peaks.narrowPeak","_summits.bed"]
			)
		)
	
	final_output.extend(expand(
			[
				"results/peaks/merged_by_sample/{sample_group}{ext}"
			],
			sample_group = pd.unique(units["sample_group"]),
			ext = ["_peaks.xls", "_peaks.narrowPeak","_summits.bed"]
		)
	)
	# count_tables
	final_output.extend(expand(
							[
								"results/count_tables/{experiment}.featureCounts"
							],
							experiment = pd.unique(samples["experiment"])
						)
					)
	# DEseq results
	if config["run_diff_accessibility"]:
		experiments = pd.unique(samples["experiment"])
		for e in experiments:
			final_output.extend(expand(
							[
								"results/DEseq2/{experiment}_{contrast}_results.tsv"
							],
							experiment = e, contrast = config["diff_accessibility"]["experiments"][e]["contrasts"]
						)
					)

	return final_output