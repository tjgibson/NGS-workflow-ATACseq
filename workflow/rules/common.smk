import pandas as pd


# def get_genome_fn():
# 	if config["use_spikeIn"]:
# 		return "resources/genome.fasta"
# 	else:
# 		return "resources/ref_genome.fasta"
# 		
# def get_genome_bn():
# 	if config["use_spikeIn"]:
# 		return "resources/genome"
# 	else:
# 		return "resources/ref_genome"

# read in table with sample metadata
# samples = (
#     pd.read_csv(config["samples"], sep="\t", dtype={"sample_name": str})
#     .set_index("sample_name", drop=False)
#     .sort_index()
# )



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


def get_bowtie2_input(wildcards):
	if not is_activated("mergeReads"):
		unit = units.loc[wildcards.sample]
		if all(pd.isna(unit["fq1"])):
			# SRA sample (always paired-end for now)
			accession = unit["sra"]
			if all(unit["read_format"] == "SE"):
				return expand("data/sra/se/{accession}.fastq.gz", accession=accession)
			else:
				return expand("data/sra/pe/{accession}_{read}.fastq.gz", accession=accession, read=[1,2])
		fastqs = units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
		if len(fastqs) == 2:
			return [fastqs.fq1, fastqs.fq2]
		return fastqs.fq1
	unit = units.loc[wildcards.sample]
	if all(unit["read_format"] == "SE"):
		return ["data/merged/{sample}_single.fastq.gz"]
	return ["data/merged/{sample}_1.fastq.gz", "data/merged/{sample}_2.fastq.gz"]

def get_bam_merge(wildcards):
	unit =  units[units["sample_group"] == wildcards.sample_group]
	group = unit["sample_name"]
	return expand(
		"results/aligned_reads/filtered/{group}.bam", group=group)


def get_macs2_input_narrow_se(wildcards):
	unit = units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		if all(unit["read_format"] == "SE"):
			if all(unit["peak_type"] == "narrow"):
				if all(pd.isna(unit["input"])):
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam"}
				else:
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam", "control": "results/aligned_reads/filtered/{input}.bam".format(input=unit.iloc[0].input)}

def get_macs2_input_broad_se(wildcards):
	unit = units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		if all(unit["read_format"] == "SE"):
			if all(unit["peak_type"] == "broad"):
				if all(pd.isna(unit["input"])):
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam"}
				else:
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam", "control": "results/aligned_reads/filtered/{input}.bam".format(input=unit.iloc[0].input)}

def get_macs2_input_narrow_pe(wildcards):
	unit = units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		if all(unit["read_format"] == "PE"):
			if all(unit["peak_type"] == "narrow"):
				if all(pd.isna(unit["input"])):
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam"}
				else:
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam", "control": "results/aligned_reads/filtered/{input}.bam".format(input=unit.iloc[0].input)}

def get_macs2_input_broad_pe(wildcards):
	unit = units.loc[wildcards.sample]
	if all(unit["call_peaks"]):
		if all(unit["read_format"] == "PE"):
			if all(unit["peak_type"] == "broad"):
				if all(pd.isna(unit["input"])):
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam"}
				else:
					return {"treatment": "results/aligned_reads/filtered/{sample}.bam", "control": "results/aligned_reads/filtered/{input}.bam".format(input=unit.iloc[0].input)}

def get_spikeIn_input():
	stat_files = expand(
				["results/aligned_reads/stats/{sample}_unireads.idxstats"],
				sample = units["sample_name"]
			)
	return stat_files


def get_final_output():
	final_output = []

	# coverage bigwigs for individual replicates
	final_output.extend(expand(
					[
						"results/bigwigs/coverage/individual/{sample}.bw"
					],
					sample = units["sample_name"]
				)
			)

	# coverage bigwigs for merged replicates
	final_output.extend(expand(
					[
						"results/bigwigs/coverage/merged/{sample}.bw"
					],
					sample = units["sample_name"]
				)
			)
	
	# z-score normalized bigwigs for merged replicates
	final_output.extend(expand(
					[
						"results/bigwigs/zscore_normalized/merged/{sample}.bw"
					],
					sample = units["sample_group"]
				)
			)

	 # peaks
	if any(units["call_peaks"]):
	 
		# add narrow peak output
		if any((units["peak_type"] == "narrow") & (units["read_format"] == "SE")):
			out_samples =  units[(units["peak_type"] == "narrow") & (units["read_format"] == "SE")]
			final_output.extend(expand(
					[
						"results/narrow_peaks/se/{sample}{ext}"
					],
					sample = out_samples["sample_name"],
					ext = ["_peaks.xls", "_peaks.narrowPeak","_summits.bed"]
				)
			)

		if any((units["peak_type"] == "narrow") & (units["read_format"] == "PE")):
			out_samples =  units[(units["peak_type"] == "narrow") & (units["read_format"] == "PE")]
			final_output.extend(expand(
					[
						"results/narrow_peaks/pe/{sample}{ext}"
					],
					sample = out_samples["sample_name"],
					ext = ["_peaks.xls", "_peaks.narrowPeak","_summits.bed"]
				)
			)

		if any((units["peak_type"] == "broad") & (units["read_format"] == "SE")):
			out_samples =  units[(units["peak_type"] == "broad") & (units["read_format"] == "SE")]
			final_output.extend(expand(
					[
						"results/broad_peaks/se/{sample}{ext}"
					],
					sample = out_samples["sample_name"],
					ext = ["_peaks.xls", "_peaks.broadPeak","_peaks.gappedPeak"]
				)
			)

		if any((units["peak_type"] == "broad") & (units["read_format"] == "PE")):
			out_samples =  units[(units["peak_type"] == "broad") & (units["read_format"] == "PE")]
			final_output.extend(expand(
					[
						"results/broad_peaks/pe/{sample}{ext}"
					],
					sample = out_samples["sample_name"],
					ext = ["_peaks.xls", "_peaks.broadPeak","_peaks.gappedPeak"]
				)
			)

	return final_output