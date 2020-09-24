# libraries
library(data.table)
library(dplyr)

# function to create caviar input files for one block at a given rsq threshold
create_caviar_files_for_block = function(rsq_threshold, block_num, results_map3k1, ld_map3k1, analysis_directory=""){
  block_num_col = paste0("block_num_", rsq_threshold)
  # count variants in block
  current_block = results_map3k1 %>% filter(get(block_num_col)==block_num) 
  print(paste0("block: ", block_num, ", number of variants: ", nrow(current_block)))
  # for blocks with >1 variant, create caviar files
  if (nrow(current_block)>1){
    # prepare z-score file
    z_file = current_block %>% 
      mutate(z_score = bcac_onco_icogs_gwas_beta/bcac_onco_icogs_gwas_se) %>%
      select(snp_id, z_score)
    # create r (not rsq) matrix - sign is important
    sign_betas_ordered = z_file %>% select(snp_id) %>% left_join(results_map3k1) %>% select(sign_beta) %>% pull()
    sign_matrix = sign_betas_ordered %*% t(sign_betas_ordered)
    r_matrix_unsigned = z_file %>% select(RS_number = snp_id) %>%
      left_join(ld_map3k1) %>% 
      select(z_file$snp_id) %>%
      mutate_all(sqrt) 
    r_matrix = r_matrix_unsigned*sign_matrix # multiply these together so that the correlations are signed
    # write files
    fwrite(z_file, file = paste0("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/", analysis_directory, "z_file_rsq_", rsq_threshold, "_block_", block_num, ".txt"), sep="\t",  col.names = FALSE)
    fwrite(r_matrix, file = paste0("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/", analysis_directory, "r_matrix_rsq_", rsq_threshold, "_block_", block_num, ".txt"), sep="\t",  col.names = FALSE)
  }
}

# function to create caviar input files for an rsq threshold
create_all_caviar_files = function(rsq_threshold, results_map3k1, ld_map3k1, analysis_directory=""){
  # identify number of blocks
  block_num_col = paste0("block_num_", rsq_threshold)
  highest_block_num = results_map3k1 %>% select({{block_num_col}}) %>% pull() %>% max(., na.rm = TRUE) 
  # loop through blocks and create caviar files for those with >1 variant
  for (block_number in 1:highest_block_num) {
    create_caviar_files_for_block(rsq_threshold, block_number, results_map3k1, ld_map3k1, analysis_directory)
  }
}

# function to create files for rsq thresholds of 0.2 and 0.4 for a set of input and output files
create_caviar_both_rsq_thresholds = function(results_path, ld_path, analysis_directory=""){
  # read files
  results_map3k1 = fread(results_path) %>%
    filter(!is.na(bcac_onco_icogs_gwas_beta)) %>%
    mutate(sign_beta = ifelse(bcac_onco_icogs_gwas_beta>0, 1,
                              ifelse(bcac_onco_icogs_gwas_beta<0, -1, NA)))
  ld_map3k1 = fread(ld_path)
  # create files
  create_all_caviar_files("04", results_map3k1, ld_map3k1, analysis_directory)
  create_all_caviar_files("02", results_map3k1, ld_map3k1, analysis_directory)
}

# 305kb rebion
results_path = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/map3k1_gwsig_with_hap_blocks.csv"
ld_path = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/ld_matrix_1kg_eur.csv"
create_caviar_both_rsq_thresholds(results_path, ld_path)

# 1 Mb upstream and downstream of promotor
results_path = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/map3k1_gwsig_with_hap_blocks_1M_upstream_downstream.csv"
ld_path = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/ld_matrix_1kg_eur_1M_upstream_downstream.csv"
create_caviar_both_rsq_thresholds(results_path, ld_path, "Upstream_Downstream_1M/") 
