library(dplyr)
args = commandArgs(trailingOnly=TRUE)

cluster_meta <- read.table(args[1],header=T,sep="\t")
genome_length <- read.table(args[2],header=F)
names(genome_length)<-c("gb_accession","length")
protein_meta <- read.table(args[3],header=T,sep="\t")

table2_with_gene <- genome_length %>%
  left_join(protein_meta %>% select(-gene), by = "gb_accession") %>%
  left_join(cluster_meta %>% select(protein_id, gene), by = "protein_id") %>%
  select(gb_accession, length, gene)

gene_stats <- table2_with_gene %>%
group_by(gene) %>%
summarise(
	n = n(),
	median_length = median(length, na.rm = TRUE),
	mean_length = mean(length, na.rm = TRUE),
	sd_length = sd(length, na.rm = TRUE),
	min_length = min(length, na.rm = TRUE),
	max_length = max(length, na.rm = TRUE),
	.groups = "drop"
)

write.table(gene_stats, file = "", sep = "\t", row.names = FALSE, quote = FALSE)



