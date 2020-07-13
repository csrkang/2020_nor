library(purrr)
library(dplyr)
library(here)

# Import all counts files in the counts/ folder
fileNames <- Sys.glob (here::here("counts", "*.txt"))

RNACounts <- lapply(fileNames, function(i) {
  read.table(i, sep = '\t', header = T, stringsAsFactors = F)
})

# Join data by locus tag
RNACounts <- RNACounts %>% 
  reduce(full_join, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))

# Rename columns to be SRX numbers
for (i in 7:length(colnames(RNACounts))){
  colnames(RNACounts)[i] <- 
    gsub('[^A-Z0-9]', '', colnames(RNACounts)[i])
}

# Calculate FPKM = X * 1e9 / (L * sum(X))
RNACounts <- RNACounts %>% 
  mutate_at(.funs = list(FPKM = ~.*(10^9)/(as.double(Length)*sum(.))), 
            .vars = 7:length(colnames(RNACounts)))

# Save FPKM values
write.table(RNACounts, file = here::here('counts', 'RNA_Counts.txt'), 
            quote = F, row.names = F, col.names = T, sep = '\t')
