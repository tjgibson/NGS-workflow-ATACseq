rule macs2_call_peaks_narrow_se:
	input:
		unpack(get_macs2_input_narrow_se)
	output:
		multiext("results/narrow_peaks/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/callpeak_narrow_{sample}.log"
	params: "-f BAM -g dm"
	wrapper:
		"0.77.0/bio/macs2/callpeak"

rule macs2_call_peaks_broad_se:
	input:
		unpack(get_macs2_input_broad_se)
	output:
		multiext("results/broad_peaks/{sample}",
                 "_peaks.xls",   ### required
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
	log:
		"logs/macs2/callpeak_broad_{sample}.log"
	params: "-f BAM -g dm"
	wrapper:
		"0.77.0/bio/macs2/callpeak"

rule macs2_call_peaks_narrow_pe:
	input:
		unpack(get_macs2_input_narrow_pe)
	output:
		multiext("results/narrow_peaks/{sample}",
                 "_peaks.xls",   ### required
                 ### optional output files
                 "_peaks.narrowPeak",
                 "_summits.bed"
                 )
	log: 
		"logs/macs2/callpeak_narrow_{sample}.log"
	params: "-f BAMPE -g dm"
	wrapper:
		"0.77.0/bio/macs2/callpeak"

rule macs2_call_peaks_broad_pe:
	input:
		unpack(get_macs2_input_broad_pe)
	output:
		multiext("results/broad_peaks/{sample}",
                 "_peaks.xls",   ### required
                 # these output extensions internally set the --broad option:
                 "_peaks.broadPeak",
                 "_peaks.gappedPeak"
                 )
	log:
		"logs/macs2/callpeak_broad_{sample}.log"
	params: "-f BAMPE -g dm"
	wrapper:
		"0.77.0/bio/macs2/callpeak"


