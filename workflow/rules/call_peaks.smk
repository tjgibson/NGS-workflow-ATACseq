rule macs2_call_peaks_narrow_pe:
	input:
		unpack(get_macs2_input_narrow_pe)
	output:
		multiext("results/narrow_peaks/pe/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/callpeak_narrow_{sample}.log"
	params: config["params"]["macs2_call_peaks_narrow_pe"]
	wrapper:
		"v1.1.0/bio/macs2/callpeak"
