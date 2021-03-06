# libraries
library(data.table)
library(dplyr)
library(stringr)

# create vector of all positions in enhancer regions in MCF-7 cells using Enhancer Atlas (https://academic.oup.com/nar/article/48/D1/D58/5628925)
mcf7 = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/MAP3K1_1M_Upstream_Downstream_MCF_7_EnhancerAtlas.txt", header = FALSE) %>%
  filter(startsWith(V1, ">")) %>%
  mutate(V1 = str_remove_all(V1, "^.*:|_.*$"),
         min_pos = str_remove(V1, "^.*-"),
         max_pos = str_remove(V1, "-.*$"))
         #V1 = str_replace(V1, "-", ":")) 
create_vect = function(pos_str){
  min_pos = as.numeric(str_remove(pos_str, "^.*-"))
  max_pos = as.numeric(str_remove(pos_str, "-.*$"))
  c(min_pos:max_pos)
}
enhancer_pos = unlist(lapply(mcf7$V1, create_vect))

# read files
sig_results_map3k1 = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/map3k1_gwsig_with_hap_blocks_1M_upstream_downstream.csv")
all_results_map3k1 = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/chr5_map3k1_1M_upstream_downstream.txt")
caviar_posteriors = list.files(path = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/Upstream_Downstream_1M/", pattern = "caviar_out_") %>%
  str_subset(., "post")

# include all variants from region, even without genome wide significance
results_map3k1 = full_join(all_results_map3k1, sig_results_map3k1) %>% 
  mutate(minus_log10_p = -log10(bcac_onco_icogs_gwas_P1df))

# determine whether variants are in enhancer regions in MCF-7 cells
results_map3k1 = results_map3k1 %>% mutate(mcf7_enhancer = position_b37 %in% enhancer_pos)

# function to read and row bind all caviar posteriors for an rsq threshold
combine_caviar_posts = function(rsq_threshold, mult_causal_string=""){
  post_vect = str_subset(caviar_posteriors, paste0("out_", mult_causal_string, "rsq_", rsq_threshold, "_"))
  post_list = list()
  # read caviar results for all files named in post_vect
  for (i in seq_along(post_vect)) {
    post_list[[i]] = fread(post_vect[i])
  }
  # create one file with all posteriors
  post_col = paste0("caviar_posts_", mult_causal_string, "rsq_", rsq_threshold)
  all_posts = bind_rows(post_list) %>% select(snp_id=SNP_ID, !!post_col:=`Causal_Post._Prob.`)
  print(dim(all_posts))
  return(all_posts)
}

# join to main file with genome wide significant map3k1 variants to caviar results
setwd("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/Upstream_Downstream_1M/")
map3k1_caviar = results_map3k1 %>% 
  left_join(combine_caviar_posts("02")) %>% 
  left_join(combine_caviar_posts("04")) %>%
  left_join(combine_caviar_posts("02", "mult_causal_")) %>%
  left_join(combine_caviar_posts("04", "mult_causal_")) 

head(map3k1_caviar)

# Determine priority
map3k1_caviar_prioritized = map3k1_caviar %>% 
  mutate(avg_caviar = rowMeans(select(., caviar_posts_rsq_04, caviar_posts_mult_causal_rsq_04, 
                                      caviar_posts_rsq_02, caviar_posts_mult_causal_rsq_02), na.rm = TRUE),
         avg_caviar = ifelse((is.na(avg_caviar) & (!is.na(rsq_04) | !is.na(rsq_02))), 0.051, 
                             ifelse(is.na(avg_caviar), 0, avg_caviar)),
         priority = ifelse(mcf7_enhancer & avg_caviar>0.05, 1, 
                           ifelse(mcf7_enhancer & avg_caviar>0, 2, 
                                  ifelse(avg_caviar>0.05, 3, 
                                         ifelse(avg_caviar>0, 4, 
                                                ifelse(minus_log10_p>-log10(5e-8), 5, 
                                                       ifelse(mcf7_enhancer, 6, NA))))))) %>%
  arrange(priority)
head(map3k1_caviar_prioritized)
table(map3k1_caviar_prioritized$priority)

# Write file
fwrite(map3k1_caviar_prioritized, file = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/map3k1_gw_sig_caviar_enhancers_prioritized.csv")
