rule trim_reads:
	input:
		get_NGmerge_input,
	output:
		temp("data/trimmed/{sample}_1.fastq.gz"),
		temp("data/trimmed/{sample}_2.fastq.gz"),
	log:
		"logs/trim_NGmerge/{sample}.log",
	conda:
		"../envs/NGmerge.yaml"
	params:
		extra=""
	shell:
		"NGmerge -a -e 20 -u 41 -n 8 -v -1 {input[0]} -2 {input[1]} -o data/trimmed/{wildcards.sample} 2> {log}"