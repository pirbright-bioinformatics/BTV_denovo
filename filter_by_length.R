library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

candidates <- read.delim(args[1])
stats <- read.delim(args[2])

# Extract final number after last underscore
candidates <- candidates %>%
  mutate(
    contig_length = as.numeric(str_extract(contig, "[0-9]+$"))
)

filtered <- candidates %>%
  left_join(
    stats %>% select(gene, mean_length, sd_length),
    by = c("genes" = "gene")
  ) %>%
  mutate(
    mean_length = round(mean_length, 0),
    sd_length = round(sd_length, 0),
    z = round(abs(contig_length - mean_length) / sd_length, 2)
  )

write.table(
  filtered,
  args[3],
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
