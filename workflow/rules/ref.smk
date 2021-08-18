rule get_ref_genome:
	output:
		temp("resources/ref_genome.fasta.gz"),
	log:
		"logs/get_ref_genome.log",
	conda:
		"../envs/curl.yaml"

	params:
		link=config["ref_genome"]["link"],
	cache: True
	shell:
		"curl {params.link} > {output} 2> {log}"


if config["use_spikeIn"]:
	rule get_spikeIn_genome:
		output:
			temp("resources/spikeIn_genome.fasta.gz"),
		log:
			"logs/get_spikeIn_genome.log",
		conda:
			"../envs/curl.yaml"
		params:
			link=config["spikeIn_genome"]["link"],
		cache: True
		shell:
			"curl {params.link} > {output} 2> {log}"

	rule combine_genomes:
		input:
			ref="resources/ref_genome.fasta.gz",
			spikeIn="resources/spikeIn_genome.fasta.gz",
		output:
			temp("resources/genome.fasta.gz"),
		log:
			"logs/combine_genomes.log",
		conda:
			"../envs/seqkit.yaml"
		cache: True
		shell:
			"""
			seqkit replace -p "(.+)" -r "spikeIn_\$1" -o resources/tmp_spikeIn.fasta.gz {input.spikeIn} 2> {log}
			zcat {input.ref} resources/tmp_spikeIn.fasta.gz | gzip > {output} 2>> {log}
			rm resources/tmp_spikeIn.fasta.gz
			"""
else:
		rule rename_genome:
			input:
				temp("resources/ref_genome.fasta.gz"),
			output:
				"resources/genome.fasta.gz",
			log:
				"logs/rename_genome.log",
			cache: True
			shell:
				"mv {input} {output} 2> {log}"

				
if config["filter_chroms"]:
	rule define_keep_chroms:
		input:
			genome="resources/genome.fasta.gz",
			keep_chroms=config["keep_chroms"]
		output:
			"resources/keep_chroms.bed",
		log:
			"logs/define_keep_chroms.log",
		conda:
			"../envs/seqkit.yaml"
		cache: True
		shell:
			"seqkit grep -f {input.keep_chroms} {input.genome}"
			" | seqkit fx2tab -nil"
			" | awk '{{print $1, 1, $2}}' > {output}"
				
rule bowtie2_index:
	input:
		reference="resources/genome.fasta.gz"
	output:
		multiext(
			"resources/genome",
			".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2",
		),
	log:
		"logs/bowtie2_build/build.log"
	params:
		extra="" # optional parameters
	threads: 8
	wrapper:
		"0.77.0/bio/bowtie2/build"
		
