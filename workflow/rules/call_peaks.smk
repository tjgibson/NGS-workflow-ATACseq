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
		multiext("results/peaks/merged/{experiment}",
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