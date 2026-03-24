if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("vsn")
library("vsn")
library(optparse)


option_list <- list(
  make_option(c("--input_data"), type="character", help="Path to input CSV"),
  make_option(c("--file_name"), type="character", help="Base file name"),
  make_option(c("--output"), type="character", help="Output CSV file name")
 
)

opt_parser <- OptionParser(option_list=option_list)
args <- parse_args(opt_parser)


VSN = function (data_matrix){
 

  #data_matrix <- replace(data_matrix, data_matrix == 0,1)

  library("vsn")
  fit <- vsn2(data_matrix)
  ynorm <- predict(fit, data_matrix)
  return (ynorm)
}

df_in <- read.csv(args$input_data,header = TRUE, sep = ",",row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
metadata_rows <- c("Group", "TECH_replica","qc")
df_metadata <- df_in[rownames(df_in) %in% metadata_rows,,drop = FALSE]
df_data <- df_in[!(rownames(df_in) %in% metadata_rows),,drop = FALSE]
df_data[] <- lapply(df_data, function(x) as.numeric(x))
df_data <- as.matrix(df_data)
df_normalized <- VSN(df_data)
df_out = rbind(df_normalized,df_metadata)


write.csv(df_out, file = args$output, row.names = TRUE)
