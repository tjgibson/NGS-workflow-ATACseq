rule macs2_call_atac_peaks_ind:
	input:
		treatment="results/aligned_reads/split_fragments/{sample}_small.bam"
	output:
		multiext("results/peaks/individual/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/macs2_call_atac_peaks_ind_{sample}.log"
	params: config["params"]["macs2_call_atac_peaks_ind"]
	wrapper:
		"v1.1.0/bio/macs2/callpeak"
		
rule macs2_call_atac_peaks_merged:
	input:
		treatment=get_macs2_merged_input
	output:
		multiext("results/peaks/merged_all/{experiment}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/macs2_call_atac_peaks_merged_{experiment}.log"
	params: config["params"]["macs2_call_atac_peaks_merged"]
	wrapper:
		"v1.1.0/bio/macs2/callpeak"

rule macs2_call_atac_peaks_merged_by_sample:
	input:
		treatment=get_macs2_merged_input_by_sample
	output:
		multiext("results/peaks/merged_by_sample/{sample_group}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/callpeak_merged_bySample_{sample_group}.log"
	params: 
		config["params"]["macs2_call_atac_peaks_merged"],
	wrapper:
		"v1.1.0/bio/macs2/callpeak"

rule extend_peak_summits:
	input:
		"results/peaks/merged/{experiment}_peaks.narrowPeak"
	output:
		"results/peaks/final/{experiment}.bed",
		"results/peaks/final/{experiment}.narrowPeak",
		"results/peaks/final/{experiment}.saf"
	conda:
		"../envs/zscore_normalize_bw.yaml"
	params:
		extend_width=201
	script:
		"../scripts/extend_peak_summits.R"