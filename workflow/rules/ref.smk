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


rule rename_genome:
	input:
		"resources/ref_genome.fasta.gz",
	output:
		temp("resources/genome.fasta.gz"),
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
			" |  awk -v OFS='\t' '{{print $1, 1, $2}}' > {output}"
				
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
		"v1.1.0/bio/bowtie2/build"
		
