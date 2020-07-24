# libraries
library(LDlinkR)
library(data.table)
library(dplyr)

# read MAP3K1 variants with genome-wide significance in BCAC, order by p-value, create snp_id variable
map3k1_gwsig = fread("/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/chr5_map3k1_gw_sig.txt") %>% 
  arrange(bcac_onco_icogs_gwas_P1df) %>%
  mutate(snp_id = ifelse(startsWith(phase3_1kg_id, "rs"), sub("\\:.*", "", phase3_1kg_id), paste0("chr", chr, ":", position_b37)),
         # variants without rsIDs in the original dataset looked up from dbSNP and manually added below:
         snp_id = replace(snp_id, snp_id=="chr5:56109436", "rs74762363"),
         snp_id = replace(snp_id, snp_id=="chr5:56070497", "rs948816056"),
         snp_id = replace(snp_id, snp_id=="chr5:56109723", "rs1383408172"),
         snp_id = replace(snp_id, snp_id=="chr5:56150961", "rs561074912"),
         snp_id = replace(snp_id, snp_id=="chr5:56121786", "rs558534116"),
         snp_id = replace(snp_id, snp_id=="chr5:56006428", "rs538602887"),
         snp_id = replace(snp_id, snp_id=="chr5:56228765", "rs542672129"),
         snp_id = replace(snp_id, snp_id=="chr5:56193299", "rs559760363")) 

# create LD matrix using European ancestry, write to file
ld_matrix = LDmatrix(map3k1_gwsig$snp_id, pop = "EUR", r2d = "r2", token = "e15519386932")
fwrite(ld_matrix, file = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/ld_matrix_1kg_eur.csv")

# Remove variants not in LD matrix
not_in_ld_matrix = setdiff(map3k1_gwsig$snp_id, ld_matrix$RS_number)
in_matrix_map3k1_gwsig = map3k1_gwsig %>% filter(!(snp_id %in% not_in_ld_matrix))

# function to create a single haplotype block
create_hap_block = function(start_num, rsq_threshold){
  curr_rsq=1
  next_num = start_num
  row_list = list()
  i = 1
  while(curr_rsq>rsq_threshold & next_num<=nrow(in_matrix_map3k1_gwsig)){
    start_variant = in_matrix_map3k1_gwsig$snp_id[start_num]
    next_variant = in_matrix_map3k1_gwsig$snp_id[next_num]
    curr_rsq = ld_matrix %>% filter(RS_number=={{start_variant}}) %>% select({{next_variant}}) %>% pull()
    if(curr_rsq>rsq_threshold){
      curr_row = data.frame(top_variant = start_variant, variant = next_variant, rsq = curr_rsq)
      row_list[[i]] = curr_row
    } 
    next_num = next_num + 1
    i = i+1
  }
  return(bind_rows(row_list))
}

# function to combine blocks for all variants
combine_blocks = function(rsq_threshold){
  block_list = list()
  top_variant_num=1
  j=1
  while(top_variant_num<=nrow(in_matrix_map3k1_gwsig)){
    curr_hap = create_hap_block(top_variant_num, rsq_threshold)
    curr_hap$block_num = j
    block_list[[j]] = curr_hap
    top_variant_num = top_variant_num + nrow(curr_hap)
    j=j+1
  }
  all_blocks = bind_rows(block_list)
  names(all_blocks)[c(1,3:4)] = paste0(names(all_blocks)[c(1,3:4)], "_", sub("\\.", "", rsq_threshold))
  return(all_blocks)
}

# determine haplotype blocks at rsq thresholds of 0.2 and 0.4
variants_blocked_02 = combine_blocks(0.2)
variants_blocked_04 = combine_blocks(0.4)

# Add -log10 p-value and haplotype blocks to variant list
sig_with_haps = map3k1_gwsig %>% mutate(minus_log10_p = -log10(bcac_onco_icogs_gwas_P1df)) %>% 
  left_join(variants_blocked_02, by=c("snp_id"="variant")) %>%
  left_join(variants_blocked_04, by=c("snp_id"="variant"))

# Write file
fwrite(sig_with_haps, file = "/zivlab/data3/jnierenberg/Data/rare_variant_functional_map3k1/independent_variants/map3k1_gwsig_with_hap_blocks.csv")

