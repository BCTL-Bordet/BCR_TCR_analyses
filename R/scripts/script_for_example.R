if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr")
if (!require("MatrixGenerics", quietly = TRUE))
  install.packages("MatrixGenerics")
if (!require("reldist", quietly = TRUE))
  install.packages("reldist")

library(dplyr)
library(MatrixGenerics)

# load example files
list_examples <- readRDS("list_examples.RDS")
tot_number_reads <- readRDS("tot_number_reads.RDS")
# "list_examples" is a named list of MiXCR outputs for each sample, reorganized for the downstream analyses. 
# The column "Proportion" was calculated as the sum of all clones (column "Clones") divided by the total number of clones for that sample

# "tot_number_reads" is a dataframe with sample names and the total number of reads mapping to genes for each sample

########### 
# First run functions stored in "R/functions" folder
###########

#####################################
# calculate BCR/TCR measures
measures <- lapply(list_examples, FUN= function(x) tcrParamsAll(x, todo = "all")) 

# reorganize results in a dataframe, re-assigning column names (all BCR/TCR measures) and column "sample" based on sample names
df_measures <- data.frame(matrix(unlist(measures), nrow=length(measures), byrow=T))
colnames_bcr_tcr <- gsub(" |[-]", "_", names(measures[[1]]))
colnames(df_measures) <- colnames_bcr_tcr  
df_measures$sample <- names(list_examples)


df_measures <- df_measures %>%
  dplyr::left_join(tot_number_reads, by = "sample")
df_measures <- df_measures %>%
  select(sample, tot_number_reads_mapping_to_genes, everything())

# compute number of IG/TR reads normalized by total number of reads mapping to genes (multiplied by a factor of 1000)

df_measures$IG_Nreads_NORM <- (df_measures$IG_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000
df_measures$IGH_Nreads_NORM <- (df_measures$IGH_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000 
df_measures$IGK_Nreads_NORM <- (df_measures$IGK_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000 
df_measures$IGL_Nreads_NORM <- (df_measures$IGL_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000 

df_measures$TR_Nreads_NORM <- (df_measures$TR_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000
df_measures$TRA_Nreads_NORM <- (df_measures$TRA_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000 
df_measures$TRB_Nreads_NORM <- (df_measures$TRB_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000 
df_measures$TRD_Nreads_NORM <- (df_measures$TRD_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000 
df_measures$TRG_Nreads_NORM <- (df_measures$TRG_Nreads / df_measures$tot_number_reads_mapping_to_genes) * 1000 