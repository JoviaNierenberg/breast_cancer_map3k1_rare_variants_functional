# libraries
library(data.table)
library(dplyr)
library(stringr)

# read files
results_map3k1 = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/map3k1_gwsig_with_hap_blocks.csv")
caviar_posteriors = list.files(path = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/", pattern = "caviar_out_rsq_") %>%
  str_subset(., "post")
all_enhancers = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/fuma/BCAC_MAP3K1_All_GW_Sig_All_Enhancers_2020_07_28/ciSNPs.txt") %>%
  mutate(enh_present = 1) %>%
  reshape(., idvar = c("uniqID", "rsID", "chr", "pos", "reg_region", "type"), 
          timevar = "tissue/cell",
          direction = "wide") %>%
  rename_with(~ str_replace_all(., "enh_present.", "")) %>%
  mutate(breast_enh = select(., E027:E028) %>%  rowSums(na.rm = TRUE), 
         num_enh = select(., 7:117) %>% rowSums(na.rm = TRUE)) %>%
  select(chr, position_b37=pos, reg_region, type, breast_enh, num_enh, everything())

# function to read and row bind all caviar posteriors for an rsq threshold
combine_caviar_posts = function(rsq_threshold){
  post_vect = str_subset(caviar_posteriors, paste0("_", rsq_threshold, "_"))
  post_list = list()
  # read caviar results for all files named in post_vect
  for (i in seq_along(post_vect)) {
    post_list[[i]] = fread(post_vect[i])
  }
  # create one file with all posteriors
  post_col = paste0("caviar_posts_rsq_", rsq_threshold)
  all_posts = bind_rows(post_list) %>% select(snp_id=SNP_ID, !!post_col:=`Causal_Post._Prob.`)
  print(dim(all_posts))
  return(all_posts)
}

# join to main file with genome wide significant map3k1 variants to caviar results and enhancers
setwd("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/")
map3k1_caviar_fuma = results_map3k1 %>% 
  left_join(combine_caviar_posts("02")) %>% 
  left_join(combine_caviar_posts("04")) %>%
  left_join(all_enhancers) 

# write file
fwrite(map3k1_caviar_fuma, file = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/map3k1_gw_sig_caviar_fuma.csv")


