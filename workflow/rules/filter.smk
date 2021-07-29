rule filter_multireads:
    input:
        "results/mapped/{sample}.bam"
    output:
        temp("results/aligned_reads/filtered/{sample}.bam")
    log:
        "{sample}.log"
    params:
        extra="-bh -q 30" # optional params string
    wrapper:
        "0.77.0/bio/samtools/view"

# rule filter_chroms:
# 
rule samtools_sort:
    input:
       "results/aligned_reads/filtered/{sample}.bam"
    output:
        "results/aligned_reads/sorted/{sample}.bam"
    params:
        extra = "",
    threads:  # Samtools takes additional threads through its option -@
        8     # This value - 1 will be sent to -@.
    wrapper:
        "0.77.0/bio/samtools/sort"

rule samtools_index:
    input:
        "results/aligned_reads/sorted/{sample}.bam"
    output:
        "results/aligned_reads/sorted/{sample}.bam.bai"
    log:
        "logs/samtools_index/{sample}.log"
    params:
        "" # optional params string
    threads:  # Samtools takes additional threads through its option -@
        4     # This value - 1 will be sent to -@
    wrapper:
        "0.77.0/bio/samtools/index"

# rule align_stats:
# 
# rule filter_stats: 
# 
# rule summarize_read_processing:
