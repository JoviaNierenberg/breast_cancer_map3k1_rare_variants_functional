# libraries
library(data.table)
library(dplyr)

# read files
results_map3k1 = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/map3k1_gwsig_with_hap_blocks.csv") %>%
  filter(!is.na(bcac_onco_icogs_gwas_beta)) %>%
  mutate(sign_beta = ifelse(bcac_onco_icogs_gwas_beta>0, 1,
                          ifelse(bcac_onco_icogs_gwas_beta<0, -1, NA)))
ld_map3k1 = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/ld_matrix_1kg_eur.csv")

# prepare z-score file
z_file = results_map3k1 %>% filter(block_num_04==1) %>% 
  mutate(z_score = bcac_onco_icogs_gwas_beta/bcac_onco_icogs_gwas_se) %>%
  select(snp_id, z_score)
dim(z_file)

# create r (not rsq) matrix - sign is important
sign_betas_ordered = z_file %>% select(snp_id) %>% left_join(results_map3k1) %>% select(sign_beta) %>% pull()
sign_matrix = sign_betas_ordered %*% t(sign_betas_ordered)
r_matrix_unsigned = z_file %>% select(RS_number = snp_id) %>%
  left_join(ld_map3k1) %>% 
  select(z_file$snp_id) %>%
  mutate_all(sqrt) 
r_matrix = r_matrix_unsigned*sign_matrix # multiply these together so that the correlations are signed

# write files
fwrite(z_file, file = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/z_file_rsq_04_block_1.txt", sep="\t",  col.names = FALSE)
fwrite(r_matrix, file = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/caviar/r_matrix_rsq_04_block_1.txt", sep="\t",  col.names = FALSE)

# will need to run the above for each haplotype with >1 variant

# CAVIAR command (without last few parameters), from project caviar directory
# CAVIAR -o caviar_out_rsq_04_block_1 -l r_matrix_rsq_04_block_1.txt -z z_file_rsq_04_block_1.txt

