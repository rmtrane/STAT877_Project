library(tidyverse)

#setwd("Documents/UW-Madison/STAT877/project/")

all_files <- tibble(file = list.files("GSE63818_RAW", full.names = T))

all_data <- all_files %>% 
  mutate(Expr = map(file, read_delim, delim = "\t")) %>% 
  #head() %>% 
  filter(row_number() < 320) %>% 
  print()

# all_data %>% 
#   mutate(nrows = map_dbl(Expr, nrow)) %>% 
#   pull(nrows) %>% 
#   unique()

joined_data <- plyr::join_all(all_data$Expr, type = "left")

data_matrix <- joined_data %>% 
  column_to_rownames('tracking_id') %>% 
  as.matrix()

library(SingleCellExperiment)
library(SCnorm)
library(scDD)

covars <- tibble(coln = colnames(joined_data)[-1]) %>% 
  separate(coln, into = c("Gender", "condition", "Age", "Embryo", "cell"))

SCE <- SingleCellExperiment(assays = list(counts = data_matrix),
                            colData = covars)

conds <- colData(SCE)$condition

pdf(file = "count_depth.pdf", height = 5, width = 7)
count_depth_ests <- plotCountDepth(Data = SCE, Conditions = conds,
                                   FilterCellProportion = 0.1, NCores = 3)
dev.off()

count_data <- SingleCellExperiment::counts(SCE)
# Total Count normalization, Counts Per Million, CPM.
normalized_data <- t((t(count_data) / colSums(count_data)) *
                            mean(colSums(count_data)))

pdf(file = "normalized_count_depth.pdf", height = 5, width = 7)
count_depth_norm_ests <- plotCountDepth(Data = SCE, 
                                        NormalizedData = normalized_data,
                                        Conditions = conds,
                                        FilterCellProportion = 0.1, NCores = 3)
dev.off()

## SCnorm

pdf("SCnorm_K_evaluation.pdf")
par(mfrow = c(2,2))

SCnormalized_data <- SCnorm(Data = SCE,
                            Conditions = conds,
                            PrintProgressPlots = T,
                            FilterCellNum = 10,
                            K = 1,
                            NCores = 3,
                            reportSF = T)
dev.off()

prior_param=list(alpha=0.01, mu0=0, s0=0.01, a0=0.01, b0=0.01)

scDD_res <- scDD(SCdat = SCnormalized_data, 
                 prior_param = prior_param,
                 testZeroes = FALSE)

write_rds(c("scDD_res", "SCnormalized_data", "SCE"), path = "data.Rds")

## 
all_DD <- metadata(scDD_res)$Genes %>% 
  arrange(nonzero.pvalue.adj) %>% 
  filter(row_number() <= 30) %>% 
  arrange(DDcategory)

SVs <- list()

for (gg in as.character(all_DD$gene)){
  SVs[[gg]] <- sideViolin(normcounts(scDD_res)[gg,], cond = conds,
                            title.gene = paste0(gg, " (", filter(all_DD, gene == gg)$DDcategory, ")"))
}

pdf("violin_plots.pdf", width = 5, height = 6, onefile = T)
SVs
dev.off()

