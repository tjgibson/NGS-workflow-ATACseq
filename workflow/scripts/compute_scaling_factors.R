# setup ------------------------------------------------------------------------
# load required packages
library(tidyverse)

# import data ------------------------------------------------------------------
# read in all relevant stat files
stat_files <- snakemake@input[["mapping_stats"]]
names(stat_files) <- stat_files

mapping_stats <- stat_files %>%
  map_dfr(read_tsv, .id = "source", col_names = c("chromosome", "chromosome_size", "mapped_reads", "unmapped_reads")) %>%
  mutate(sample_name = gsub("_unireads.idxstats", "", basename(source)))


# add column marking reference or spike-in chromosomes
mapping_stats <- mapping_stats %>%
  mutate(reference_genome = case_when(str_detect(chromosome, "spikeIn")  ~ "spikeIn",
                                      !str_detect(chromosome, "spikeIn")  ~ "reference")) 


# compute % reads mapping to reference and spike-in genomes
mapping_stats <- mapping_stats %>% 
  group_by(sample_name, reference_genome) %>% 
  summarise(total_reads = sum(mapped_reads)) %>%
  pivot_wider(names_from = reference_genome, values_from = total_reads) %>%
  mutate(total_mapped_reads = reference + spikeIn) %>%
  mutate(percent_reference = reference / total_mapped_reads * 100, percent_spikeIn = spikeIn / total_mapped_reads * 100)


# read in unit table
sample_table <- snakemake@config[["units"]] %>% 
  read_tsv() %>% 
  
  # remove redundant technical replicates from sample_table
  distinct(sample_name, .keep_all = TRUE)


# mark samples as inputs or IPs
sample_table <- sample_table %>% 
  mutate(chip_sample_type = case_when(
    (sample_name %in% .$input) ~ "input",
    !(sample_name %in% .$input) ~ "IP"
  )) %>% 
  mutate(has_input = !is.na(input))


# create new column indicating input/IP pairs
tmp <- sample_table %>% select(input, sample_name) %>% dplyr::rename(IP = sample_name, sample_name = input)
sample_table <- sample_table %>% 
  left_join(tmp, by = "sample_name") %>% 
  mutate(IP_group = case_when(
    is.na(IP) ~ sample_name,
    !is.na(IP) ~ IP
  ))




# prepare sample_table and mapping_stats tables for join
individual_mapping_stats <- mapping_stats %>% 
  select(sample_name, percent_reference, percent_spikeIn)

# join sample_table and mapping stats
individual_scaling_factors <- sample_table %>%
  left_join(individual_mapping_stats, by = "sample_name")


# compute normalization factors for individual samples -------------------------
sample_has_input <- individual_scaling_factors %>% 
  filter(chip_sample_type == "IP") %>% 
  select(has_input) %>% 
  pull()

if (all(sample_has_input)) {
  
  individual_scaling_factors <- individual_scaling_factors %>% 
    pivot_wider(names_from = chip_sample_type, values_from = percent_spikeIn, id_cols = IP_group) %>% 
    mutate(enrichment = IP / input) %>%
    mutate(scaling_factor = 1 / enrichment) %>%
    mutate(norm_scaling_factor = scaling_factor / max(scaling_factor))
  
} else {
  warning("not all IPs have corresponding input sample, calculating scaling factors based on %spikeIn for IP samples")
  individual_scaling_factors <- individual_scaling_factors %>% 
    pivot_wider(names_from = chip_sample_type, values_from = percent_spikeIn, id_cols = IP_group) %>% 
    mutate(scaling_factor = 1 / IP) %>%
    mutate(norm_scaling_factor = scaling_factor / max(scaling_factor))
}

individual_scaling_factors %>% 
  write_tsv(snakemake@output[[1]])

# prepare data for merged samples -----------------------------
merged_mapping_stats <- mapping_stats %>% 
  select(sample_name, reference, spikeIn)

# combine replicate reads
merged_mapping_stats <- sample_table %>% 
  left_join(merged_mapping_stats, by = "sample_name") %>% 
  group_by(sample_group) %>% 
  summarise(spikeIn = sum(spikeIn), reference = sum(reference))




# create new column indicating input/IP pairs
tmp <- sample_table %>% select(input, sample_group) %>% dplyr::rename(IP = sample_group, sample_name = input)
sample_table <- sample_table %>% 
  select(-IP, -IP_group) %>% 
  left_join(tmp, by = "sample_name") %>% 
  mutate(IP_group = case_when(
    is.na(IP) ~ sample_group,
    !is.na(IP) ~ IP
  ))



# join sample_table and mapping stats for merged reads
merged_scaling_factors <- merged_mapping_stats %>% 
  left_join(select(sample_table, sample_group, chip_sample_type,IP_group, has_input), by = "sample_group") %>% 
  distinct(.keep_all = TRUE) %>%
  mutate(total_mapped_reads = reference + spikeIn) %>%
  mutate(percent_reference = reference / total_mapped_reads * 100, percent_spikeIn = spikeIn / total_mapped_reads * 100)

# compute normalization factors for merged samples -------------------------
sample_has_input <- merged_scaling_factors %>% 
  filter(chip_sample_type == "IP") %>% 
  select(has_input) %>% 
  pull()
if (all(sample_has_input)) {
  
  merged_scaling_factors <- merged_scaling_factors %>% 
    pivot_wider(names_from = chip_sample_type, values_from = percent_spikeIn, id_cols = IP_group) %>% 
    mutate(enrichment = IP / input) %>%
    mutate(scaling_factor = 1 / enrichment) %>%
    mutate(norm_scaling_factor = scaling_factor / max(scaling_factor))
  
} else {
  warning("not all IPs have corresponding input sample, calculating scaling factors based on %spikeIn for IP samples")
  merged_scaling_factors <- merged_scaling_factors %>% 
    pivot_wider(names_from = chip_sample_type, values_from = percent_spikeIn, id_cols = IP_group) %>% 
    mutate(scaling_factor = 1 / IP) %>%
    mutate(norm_scaling_factor = scaling_factor / max(scaling_factor))
}

merged_scaling_factors %>% 
  write_tsv(snakemake@output[[2]])
