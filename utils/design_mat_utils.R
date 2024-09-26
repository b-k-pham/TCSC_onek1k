library(tidyverse)

generate_design_mat <- function(ct,source = 'onek1k'){
  if (source == '1k1k'){
    source = 'onek1k'
  }
  outdir = 'TCSC/predicted_expression/onek1k/'
  #tissues <- fread("/n/groups/price/tiffany/subpheno/AllGTExTissues_restore/v8_all_tissues.txt", header = F)$V1
  #tissues <- fread("/n/groups/price/tiffany/subpheno/AllGTExTissues_restore/v8_meta_tissues.txt", header = F)$V1
  #y <- fread("prep_TCSC/TissueGroups_onek1k.txt", header = T)
  #y <- fread('TCSC/analysis/TissueGroups_onek1k_and_gtex_trunc.txt')
  #groups <- unique(y$MetaTissues)
  #tissues <- groups[-which(groups == "Remove")] #23
  #tissues <- groups
  tissue1 <- ct #make this into a function and represent values smaller than 48 bits? 
  print(tissue1)
  
  if (source == 'onek1k'){
    files1 <- list.files(paste0("prep_TCSC/score_files/",ct),pattern = "\\.profile$",full.names = T)
  } else{
    files1 <- list.files(paste0("prep_TCSC/score_files_gtex_meta/",ct),pattern = "\\.profile$",full.names = T)
  }
  
  #files1 <- acc_TWAS_df[,4]
  #files1 <- list.files(paste0("score_files/",tissue1),pattern = "\\.profile$") #first file is error directory
  #files1 <- list.files(paste0("score_files/",tissue1),pattern = "\\.profile$")
  #files1 <- files1[-1]
  #transcripts1 <- sapply(1:length(files1), function(x) strsplit(files1[x], split = paste0("plink2.",tissue1,"."))[[1]][2]) 
  #transcripts1 <- sapply(1:length(files1), function(x) strsplit(transcripts1[x], split = ".profile")[[1]][1])
  #print(files1)
  transcripts1 <- sapply(1:length(files1), function(x) strsplit(files1[x], split = ".score.profile")[[1]][1])
  if (source == 'gtex'){
    transcripts1 <- sapply(1:length(files1), function(x) strsplit(transcripts1[x], split = paste0(tissue1,'.'))[[1]][2])
  } else{
    transcripts1 <- sapply(1:length(files1), function(x) strsplit(transcripts1[x], split = paste0(tissue1,'.'))[[1]][3])
  }
  
  #print(transcripts1)
  #write.table(transcripts1, file = paste0(outdir,sstat_abbrev,"_TranscriptsIn",tissue1,"Model.txt"), quote = F, row.names = F, col.names = F, sep = "\t")
  
  #mat1 <- lapply(1:length(files1), function(x) fread(paste0("/n/scratch3/users/t/tja2/GTEx_v8/1KG_res/",tissue1,"/",files1[x]), header = T)$SCORE)
  #mat1 <- lapply(1:length(files1), function(x) fread(paste0("/n/scratch3/users/t/tja2/GTEx_v8/1KG_res/",tissue1,"/",files1[x]), header = T)$SCORE)
  mat1 <- lapply(1:length(files1), function(x) fread(paste0(files1[x]), header = T)$SCORE)
  #mat1 <- lapply(1:length(files1), function(x) fread(paste0("score_files/",tissue1,"/",files1[x]), header = T)$SCORE)
  df1 <- do.call(cbind,mat1) # each row is a person, each column is a gene
  rm("mat1")
  print("Made design matrix!") 

  print(length(files1))
  print(length(transcripts1))
  save(df1, file = paste0("TCSC/predicted_expression/onek1k/",tissue1,"_onek1k.RData"))
  save(transcripts1, file = paste0("TCSC/predicted_expression/onek1k/",tissue1,"_onek1k_transcripts.RData")) # This is the col names of the matrix (the genes in the design matrix)
  print(object.size(df1))
}