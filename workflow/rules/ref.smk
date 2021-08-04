if config["ref_genome_source"] == "ensembl":
rule get_ensembl_ref_genome:
	output:
		temp("resources/ref_genome.fasta"),
	log:
		"logs/get_ref_genome.log",
	conda:
		"../envs/curl.yaml"
	params:
		database=config["ensembl_ref_genome"]["database"]
		species=config["ensembl_ref_genome"]["species"],
        datatype="dna",
        build=config["ensembl_ref_genome"]["build"],
        release=config["ensembl_ref_genome"]["release"],
	cache: True
	script:
		"scripts/retrieve_ref.py"


if config["ref_genome_source"] == "other":
rule get_other_ref_genome:
	output:
		temp("resources/ref_genome.fasta"),
	log:
		"logs/get_ref_genome.log",
	conda:
		"../envs/curl.yaml"

	params:
		link=config["other_ref_genome"]["link"],
	cache: True
	shell:
		"curl {params.link} | pigz -dc > {output} 2> {log}"


if config["use_spikeIn"]:
	if config["spikeIn_genome_source"] == "ensembl":
		rule get_ensembl_spikeIn_genome:
			output:
				temp("resources/spikeIn_genome.fasta"),
			log:
				"logs/get_ref_genome.log",
			conda:
				"../envs/curl.yaml"
			params:
				database=config["ensembl_spikeIn_genome"]["database"]
				species=config["ensembl_spikeIn_genome"]["species"],
				datatype="dna",
				build=config["ensembl_spikeIn_genome"]["build"],
				release=config["ensembl_spikeIn_genome"]["release"],
			cache: True
			script:
				"scripts/retrieve_ref.py"

	if config["spikeIn_genome_source"] == "other":
		rule get_other_spikeIn_genome:
			output:
				temp("resources/spikeIn_genome.fasta"),
			log:
				"logs/get_spikeIn_genome.log",
			conda:
				"../envs/curl.yaml"
			params:
				link=config["other_spikeIn_genome"]["link"],
			cache: True
			shell:
				"curl {params.link} | pigz -dc > {output} 2> {log}"

	rule combine_genomes:
		input:
			ref="resources/ref_genome.fasta",
			spikeIn="resources/spikeIn_genome.fasta",
		output:
			temp("resources/genome.fasta"),
		log:
			"logs/combine_genomes.log",
		cache: True
		shell:
			"""
			sed -e 's/>/>ref_/' {input.ref} > {input.ref} 2> {log}
			sed -e 's/>/>spikeIn_/' {input.spikeIn} > {input.spikeIn} 2>> {log}
			cat {input.ref} {input.spikeIn} > {output} 2>> {log}
			rm {input.ref} {input.spikeIn} 2>> {log}
			"""
else:
		rule rename_genome:
			input:
				"resources/ref_genome.fasta",
			output:
				"resources/genome.fasta",
			log:
				"logs/rename_genome.log",
			cache: True
			shell:
				"mv {input} {output} 2> {log}"

rule bowtie2_index:
	input:
		reference="resources/genome.fasta"
	output:
		multiext(
			"resources/genome",
			".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
		),
	log:
		"logs/bowtie2_build/build.log"
	params:
		extra=""  # optional parameters
	threads: 8
	wrapper:
		"0.77.0/bio/bowtie2/build"
		
