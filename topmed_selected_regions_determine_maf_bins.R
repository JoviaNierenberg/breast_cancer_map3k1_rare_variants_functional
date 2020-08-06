library(data.table)
library(dplyr)
library(epiDisplay)

# read and concatenate files with topmed variants for each region
setwd("/zivlab/data3/shuntsman/jovia/")
files = list.files()
list_of_read_files = list()
for (i in 1:length(files)) {
  curr_file = fread(files[i])
  list_of_read_files[[i]] = curr_file
}
files_concatenated = bind_rows(list_of_read_files)

# create MAF bins
maf_bins = files_concatenated %>% mutate(maf_bins = ifelse(MAF<0.0004, "singleton",
                                                           ifelse(MAF<0.005, ">singleton but <0.005", 
                                                                  ifelse(MAF<0.01, ">=0.005-0.01",
                                                                         ifelse(MAF<0.05, ">=0.01-0.05", 
                                                                                ifelse(MAF<0.1, ">=0.05-0.1", ">=0.1"))))))

# determine bin percentages
maf_percentages = tab1(maf_bins$maf_bins, sort.group = "decreasing", graph = FALSE, cum.percent = FALSE)

# write file
fwrite(maf_bins, file = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/topmed_variants_in_selected_regions.csv")
