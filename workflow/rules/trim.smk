rule trim_reads:
	input:
        get_NGmerge_input,
	output:
		temp("data/trimmed/{sample}_{read}.fastq.gz"),
	log:
		"logs/trim_NGmerge/{sample}.log",
	conda:
		"../envs/NGmerge.yaml"
	params:
		extra=""
	wildcard_constraints:
        read="single|1|2",
	shell:
		"NGmerge -a -e 20 -u 41 -n 8 -v -1 {input[0]} -2 {input[1]} -o data/trimmed/{wildcards.sample} 2> {log}"