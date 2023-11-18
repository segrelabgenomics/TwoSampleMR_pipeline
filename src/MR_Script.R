#!/usr/bin Rscript

# Authors: Skanda Rajasundaram, Puja Mehta 
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA, USA
# Date: November 16, 2023

#Load required packages
library("data.table")
library("dplyr")
library("tidyr")
library("foreign")
library("tibble")
library("metafor")
library("meta")
library("survival")
library("ggplot2")
library("plyr")
library("gridExtra")
library("gtable")
library("grid")
library("tidyverse")
library("stringr")
library("coloc")
library("devtools")
library("glmnet")
library("MendelianRandomization")
library("TwoSampleMR")


args = commandArgs(trailingOnly=TRUE)

check_args <- function(x) {
  n_args <- length(list(x))
  # Check if the number of arguments is as expected
  if (n_args != 9) {
    print("Please provide the required arguments in the following order: Trait name, Gene ID, Gene Symbol, Tissue, QTL type, p-Value cutoff, File name of the trait, path to the QTL file and QTL file extension")
    stop("Incorrect number of arguments. Please provide exactly 9 arguments.")
    
  }
}
check_args(args)

cat(args, sep = "\n")
#GWAS should be in build 38
gwas_name = args[1]
gene_ID = args[2]
gene_symb = args[3]
tissue = args[4]
qtl_type = args[5]
pvalue_cutoff = args[6] 
gwas_file = args[7]
path_to_QTL = args[8]
file_extension = args[9]

gwas <- fread(gwas_file)
setnames(gwas, old=c("effect", "gwas_p_value"), new=c("Effect", "pval"))
gwas$OA <- chartr("a-zA-Z", "A-Za-z", gwas$Other_allele)
gwas$EA <- gwas$effect_allele
head(gwas)
gwas <- gwas[!duplicated(gwas$SNP)]

write.csv(gwas, paste0("../data/",gene_symb,"_",gene_ID,"_",tissue,"_",qtl_type,"_",strsplit(gwas_name," ","_"),"_gwas.csv"), row.names = F,  quote = F)
gc()


QTL_file <-fread(paste0(path_to_QTL,tissue,file_extension))
gc()

#Filtering QTL file on the gene ID
QTL_file_subsest <- QTL_file %>% filter(gene_id == gene_ID)
gc()


split_variant <- data.frame(str_split(QTL_file_subsest$variant_id, "_", simplify = T)[,c(1:4)])
colnames(split_variant) <- c("chr","pos","OA","EA")

QTL_file_subsest <- data.frame(cbind(QTL_file_subsest,split_variant))
QTL_file_subsest$chrpos <- paste(QTL_file_subsest$chr, QTL_file_subsest$pos, sep=":")

gwas$chrpos <- paste(gwas$chr, gwas$pos, sep=":")
gwas_subset=subset(gwas, select = c(chrpos,SNP))

#Merge the QTL data with the GWAS summary statistics on Chr:pos
GWAS_QTL_merge <- data.frame(merge(QTL_file_subsest, gwas_subset, by = "chrpos", all.x = T))
gc()

set.seed(1010)

#Rename columns and save the file after removing duplicate SNPs
setnames(GWAS_QTL_merge, old=c("slope", "slope_se", "pval_nominal"), new=c("Effect", "StdErr", "Pval"))
GWAS_QTL_merge <- GWAS_QTL_merge[!duplicated(GWAS_QTL_merge$SNP),]

#Load the saved file and filter according to a P-value threshold
GWAS_QTL_merge_subset <- subset(GWAS_QTL_merge, Pval< pvalue_cutoff)
write.csv(GWAS_QTL_merge_subset, paste0("../tmp_data/",gene_symb,"_",gene_ID,"_",tissue,"_",qtl_type,"_",gwas_name,"_sig_pair_qtls_filtered.csv"), row.names = F,  quote = F)

options(warn = -1) 
gene_unclumped <- read_exposure_data(
  filename = paste0("../tmp_data/",gene_symb,"_",gene_ID,"_",tissue,"_",qtl_type,"_",gwas_name,"_sig_pair_qtls_filtered.csv"),
  sep = ",",
  snp_col = "SNP",
  chr_col = "chr",
  pos_col = "pos",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "EA", 
  other_allele_col = "OA",
  pval_col = "Pval")
system(paste0("rm ",paste0("../tmp_data/",gene_symb,"_",gene_ID,"_",tissue,"_",qtl_type,"_",gwas_name,"_sig_pair_qtls_filtered.csv")))

# Clump at r2 = 0.1 
gene_iv <- clump_data(
  gene_unclumped,
  clump_r2 = 0.1,
  pop = "EUR")

gwas_outcome <- read_outcome_data(
  filename = paste0("../data/",gene_symb,"_",gene_ID,"_",tissue,"_",qtl_type,"_",strsplit(gwas_name," ","_"),"_gwas.csv"),
  sep = ",",
  snp_col = "SNP",
  chr_col = "chr",
  pos_col = "pos",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "EA",
  other_allele_col = 'OA',
  pval_col = "pval")

system(paste0("rm ",paste0("../data/",gene_symb,"_",gene_ID,"_",tissue,"_",qtl_type,"_",strsplit(gwas_name," ","_"),"_gwas.csv")))


# Harmonise data
gene_gwas_harmonised <- harmonise_data(gene_iv, gwas_outcome)
header_num_variants <- as.numeric(length(unique(gene_gwas_harmonised$SNP)))
 
# MR analysis
gene_gwas_res <- mr(gene_gwas_harmonised, method_list=c("mr_ivw", "mr_simple_median", "mr_weighted_median", "mr_egger_regression"))

if(length(gene_gwas_res$method) == 0){
  gene_gwas_res <- mr(gene_gwas_harmonised)
}

#Calculate 95% CIs
gene_gwas_res$b
if(tolower(colnames(gwas)) %in% c("beta","effect")){
  print("Converting Beta to OR for IWD/Wald's test")
  gene_gwas_res$b = exp(gene_gwas_res$b)
}
gene_gwas_res$LCI <- gene_gwas_res$b-1.96*gene_gwas_res$se
gene_gwas_res$UCI <- gene_gwas_res$b+1.96*gene_gwas_res$se


#MR PRESSO Global test of heterogeneity Pval	
print("Obtaining MR PRESSO Global test of heterogeneity Pval")
tryCatch({heterogeneity <- mr_heterogeneity(gene_gwas_harmonised)},
         error = function(e){
	   print("MR PRESSO Global test of heterogeneity failed to run")
	   heterogeneity <- "NA"})

#Perform MR-Egger Intercept test and report P-value
print("Performing MR-Egger Intercept test and reporting P-value")
tryCatch({pleiotropy_res <- mr_pleiotropy_test(gene_gwas_harmonised)},
         error = function(e){
	 	print("MR-Egger Test failed to run")
	 	pleiotropy_res <- "NA"})
if(nrow(pleiotropy_res) == 0){
  intercept <- "NA"
}else{
  print("MR-Egger Intercept Test completed")
  intercept = pleiotropy_res$pval
}

#Perform MR-PRESSO as an additional sensitivity analysis
print("Performing MR-PRESSO as an additional sensitivity analysis")
tryCatch({presso_res <- run_mr_presso(gene_gwas_harmonised)},
         error = function(e){
	   print("MR-PRESSO Test failed to run")
	   })
if(exists("presso_res")){
  print("MR-PRESSO Test completed")
  }else{
    presso_res <- NA
   }

make_result_df <- data.frame(matrix(NA, nrow = 1, ncol = 6))
colnames(make_result_df) <- c("GWAS (outcome)","Gene","Tissue","molQTL type","Intron excision cluster ID (start:end sites)","No. of variants")
row <- c(gwas_name,gene_symb,tissue,qtl_type, 
          ifelse(qtl_type == "sQTL",gene_ID,""),
          header_num_variants)
make_result_df <- data.frame(rbind(make_result_df,row))
make_result_df<- na.omit(make_result_df)


gene_gwas_res$method <- ifelse((gene_gwas_res$method) == "Inverse variance weighted" | 
                                 (gene_gwas_res$method) == "Wald ratio", 
       "MR IVW/Wald ratio", (gene_gwas_res$method))
gene_gwas_res$method <- ifelse((gene_gwas_res$method) == "Simple median", 
              "MR Simple Median",(gene_gwas_res$method))
gene_gwas_res$method <- ifelse((gene_gwas_res$method) == "Weighted median", 
              "MR Weighted Median",(gene_gwas_res$method))
gene_gwas_res <- data.frame(gene_gwas_res[,c("method","b","se","pval","LCI","UCI")])

a <- unlist(gene_gwas_res)[c(1:nrow(gene_gwas_res))]
if(length(a) > 1){
  names <- paste0(rep(a,nrow(gene_gwas_res)), " ",names(unlist(gene_gwas_res)))
  names <- str_remove_all(names, "method[:digit:]")
  names <- str_replace_all(names, "b[:digit:]","Beta or OR")
  names <- str_replace_all(names, "se[:digit:]","SE")
  names <- str_replace_all(names, "pval[:digit:]","Pval")
  names <- str_replace_all(names, "LCI[:digit:]","or OR L 95% CI")
  names <- str_replace_all(names, "UCI[:digit:]","or OR U 95% CI")
  gene_gwas_res_unlist <- unlist(gene_gwas_res)
  names(gene_gwas_res_unlist) <- names
}else{
  names <- paste0(rep(a,nrow(gene_gwas_res)), " ",names(unlist(gene_gwas_res)))
  names <- str_remove_all(names, "method")
  names <- str_replace_all(names, "b","Beta or OR")
  names <- str_replace_all(names, "se","SE")
  names <- str_replace_all(names, "pval","Pval")
  names <- str_replace_all(names, "LCI","or OR L 95% CI")
  names <- str_replace_all(names, "UCI","or OR U 95% CI")
  gene_gwas_res_unlist <- unlist(gene_gwas_res)
  names(gene_gwas_res_unlist) <- names
}

gene_gwas_res_unlist <- gene_gwas_res_unlist[sort(names)]
make_result_df <- data.frame(cbind(make_result_df,t(gene_gwas_res_unlist)))
colnames(make_result_df) <- str_replace_all(colnames(make_result_df) , "\\.\\.","\\.")
colnames(make_result_df) <- str_replace_all(colnames(make_result_df) , "\\."," ")
colnames(make_result_df) <- trimws(colnames(make_result_df) , "right")
colnames(make_result_df)

make_result_df <- make_result_df %>% dplyr::select(any_of(c("GWAS outcome","Gene","Tissue","molQTL type",
                                   "Intron excision cluster ID start end sites",
                                   "No of variants",
                                   "MR IVW Wald ratio Beta or OR",
                                   "MR IVW Wald ratio or OR L 95 CI",
                                   "MR IVW Wald ratio or OR U 95 CI",
                                   "MR IVW Wald ratio Pval",
                                   "MR Simple Median Beta or OR",	
                                   "MR Simple Median SE", "MR Simple Median Pval",
                                   "MR Weighted Median Beta or OR", 
                                   "MR Weighted Median SE","MR Weighted Median Pval", 
                                   "MR Egger Beta or OR", 
                                   "MR Egger SE", "MR Egger Pval")))

ncol(make_result_df)
#The results file should have 19 columns
if(ncol(make_result_df) != 19){
  to_add <- 19 - ncol(make_result_df)
  df_to_add <- data.frame(matrix(NA,ncol = to_add, nrow = 1))
  make_result_df <- data.frame(cbind(make_result_df, df_to_add))
}
make_result_df$intercept <- intercept
if(is.na(presso_res) == FALSE){
  print(head(presso_res))
  presso_df <- data.frame()
  presso_df <- data.frame(rbind(presso_df,presso_res[[1]][["Main MR results"]]))
  presso_df <- presso_df[which(presso_df$`MR.Analysis` == "Outlier-corrected"),]
  outlier_pval <- str_remove_all(presso_res[[1]][["MR-PRESSO results"]][["Outlier Test"]]$Pvalue,"\\>")
  outlier_pval <- str_remove_all(presso_res[[1]][["MR-PRESSO results"]][["Outlier Test"]]$Pvalue,"\\<")
  outlier_pval <- as.numeric(outlier_pval)
  outlier_pval <- length(which(outlier_pval < 0.05))
  make_result_df$Global_test_heterogeneity_pval <- str_remove_all(presso_res[[1]][["MR-PRESSO results"]][["Global Test"]]$Pvalue,"\\>")
  make_result_df$Global_test_heterogeneity_pval <- str_remove_all(presso_res[[1]][["MR-PRESSO results"]][["Global Test"]]$Pvalue,"\\<")
  make_result_df$Global_test_heterogeneity_pval <- round(as.numeric(make_result_df$Global_test_heterogeneity_pval),3)

  if(outlier_pval > 0){
    make_result_df$outlier_estimate <- presso_df$`Causal.Estimate`
    make_result_df$outlier_estimate_sd <- presso_df$Sd
    make_result_df$outlier_estimate_pval <- presso_df$`P.value`
    make_result_df$n_outlier <- outlier_pval
  }else{
    make_result_df$outlier_estimate <- NA
    make_result_df$outlier_estimate_sd <- NA
    make_result_df$outlier_estimate_pval <- NA
    make_result_df$n_outlier <- NA
  }
}else{
  make_result_df$Global_test_heterogeneity_pval <- NA
  make_result_df$outlier_estimate <- NA
  make_result_df$outlier_estimate_sd <- NA
  make_result_df$outlier_estimate_pval <- NA
  make_result_df$n_outlier <- NA
}


if(header_num_variants > 4){
  if(is.na(heterogeneity)){
    print("MR PRESSO Global test of heterogeneity failed to run despite having minimum number of variants")}}
if(header_num_variants > 4 & is.na(make_result_df$n_outlier) == FALSE){
    if(is.na(unique(c(make_result_df[1,c(22:24)]))[[1]][1])){
      print("Outlier tests failed to run")
      }else{
      print("Outlier tests look good")}}

if(ncol(make_result_df) == 25){
  print("Results file successfully generated")
}else{
  print("Results file is missing columns")
}

print("Results File:")
print(paste0("../results/RESULTS_",tissue,"_",gene_symb,"_",gene_ID,"_",qtl_type,"_",str_replace(gwas_name," ","_"),"_",Sys.Date(),".txt"))

write.table(make_result_df,paste0("../results/RESULTS_",
	tissue,"_",gene_symb,"_",gene_ID,"_",qtl_type,"_",
	str_replace(gwas_name," ","_"),"_",Sys.Date(),".txt"))

gc()
