library(data.table)
library(Hmisc)
library(tidyverse)
library(dplyr)

### NOT NEEDED ANYMORE

# TWAS_table_lookup <- function(row,trait){ #loads twas results and outputs a z based on gene and tissueassign group in row. The input is a dataframe of transcript and tissue assign.
#   gene = as.character(row[1])
#   i = as.numeric(row[2])
#   #print(gene)
#   #print(i)
#   source = as.character(row[3])
#   tissue = as.character(row[4])
#   
#   if (source == '1k1k'){
#     zs = fread(paste0("output/",trait,"/",tissue,"_twas.txt"),header = T)
#     out = unique(zs[match(gene,zs$ID),]$TWAS.Z)
#   } else{
#     zs = zs_df %>% filter(tissue == tissue)
#     out = unique(zs[match(gene,zs$transcript),]$TWAS.Z)
#   }
#   return(out)
# }


get_cov_alpha1alpha1_multitissue <- function(sources,cells,alpha_z,CoRegMat,nGWAS,weights1,weights2,weights3){ #add cells argument.
# weight4 is a different experiment.
# cells should be cells per ind
# sources check if 1k1k or gtex. If gtex, do not residualize things out.
y <- (alpha_z^2)/nGWAS
idx_sources = which(sources == '1k1k') # find which are from 1k1k
mod <- lm(y[idx_sources] ~ cells[idx_sources])
y[idx_sources] <- resid(mod)
w1 <- 1/weights1
w2 <- 1/weights2
w3 <- 1/weights3
mod <- summary(glm(y ~ CoRegMat, weights = w1*w2*w3)) #to constrain intercept, use summary(glm(I(y+1) ~ CoRegMat, weights = w1*w2*w3))
#mod <- summary(glm(I(y+1) ~ CoRegMat, weights = w1*w2*w3))
cov_b1b1 <- coef(mod)[-1,1]  
cov_b1b1 <- cov_b1b1[1:ncol(CoRegMat)]
return(cov_b1b1)
}


### MAIN FUNCTION TO RUN TCSC ON A SUMSTAT

run_TCSC_by_idx <- function(trait,h2g,N,coreg_path,print_statements = FALSE){
  #y <- fread("TCSC/analysis/TissueGroups_onek1k_and_gtex_trunc.txt", header = T) #one row per GTEx tissue analyzed in Amariuta et al 2022 bioRxiv
  #tissues <- unique(y$Tissues) #tissues with eQTL sample size < 320 are grouped together via meta-analysis if the correlation of marginal eQTL effects is > 0.93. 
  y <- fread('TCSC/analysis/onek1k_and_gtex_meta_tgs.txt', header = T)
  tissues <- y$MetaTissue
  #tg_sources <- y$Source
  
  #n_eqtl <- sapply(1:length(tissues), function(x) sum(y$N_EUR[which(y$MetaTissue == tissues[x])])) #find total available sample size from GTEx data after accounting for meta-analysis of tissues
  #small_tissues <- which(n_eqtl < 320)
  #normal_tissues <- c(1:length(tissues))[-small_tissues] #in primary analysis, these tissues are subsampled such that eQTL sample size = 320
  cells_per_ind = read.csv('onek1k_cells_per_ind.csv', header = F)
  colnames(cells_per_ind) = c('tissue','cells_per_ind')
  cells_per_ind = cells_per_ind %>% mutate(tissue = str_replace(tissue,' ','_'))
  
  
  load(coreg_path) #"TCSC/analysis/InputCoreg_onek1k_gtex_TCSC.RData"
  gtex <- fread("TCSC/analysis/gene_annotation.txt.gz", header = F, sep = "\t")
  
  gtex_tissues = y %>% filter(Source == 'gtex')
  gtex_tissues = gtex_tissues$MetaTissue
  
  zs_df = list()
  
  
  zs_df <- lapply(gtex_tissues, function(tissue) {
    twas_z <- fread(paste0('TCSC/twas_statistics/320EUR_GTEx/', trait, '/Marginal_alphas_', trait, '_', tissue, '.txt.gz'), header = FALSE)$V1
    transcriptsin <- fread(paste0('TCSC/weights/heritablegenes/N320/TranscriptsIn', tissue, 'Model.txt'), header = FALSE)$V1
    temp <- data.frame(transcripts = transcriptsin, TWAS.Z = twas_z)
    temp$tissue <- tissue
    return(temp)
  })
  
  zs_df = do.call(rbind,zs_df)
  
  transcript_tissue_assigns = data.frame(transcripts,tissueassign,sources) %>% mutate(tissue = tissues[tissueassign]) %>% left_join(cells_per_ind,by = 'tissue') %>% replace_na(list(x = NA, cells_per_ind = median(cells_per_ind$cells_per_ind)))
  
  dt1 = as.data.table(transcript_tissue_assigns %>% filter(sources == 'gtex'))
  dt2 = as.data.table(zs_df)
  
  
  gtex_dt = dt1 %>% left_join(dt2, by = c('transcripts','tissue'))
  
  dt1 = as.data.table(transcript_tissue_assigns %>% filter(sources == '1k1k'))
  
  twas_file_paths = list.files(paste0('output/',trait,'/'),pattern = '_twas.txt',full.names = T)
  twas_dt = list()
  for (i in 1:length(twas_file_paths)){
    twas_dt[[i]] = fread(twas_file_paths[i])
  }
  
  
  twas_dt = do.call('rbind',twas_dt) %>% select(PANEL,ID,TWAS.Z)
  colnames(twas_dt) = c('tissue','transcripts','TWAS.Z')
  
  onek1k_dt = dt1 %>% left_join(twas_dt, by = c('transcripts','tissue'))
  
  #(nrow(gtex_dt) + nrow(onek1k_dt)) == nrow(transcript_tissue_assigns) # TRUE
  
  z_dt = rbind(gtex_dt,onek1k_dt) # this takes 3 seconds
  
  colnames(z_dt)[6] = 'z'
  
  
  alpha_z = z_dt$z
  sources = z_dt$sources
  cells_per_ind = z_dt$cells_per_ind
  
  w <- which(alpha_z^2 > max(80,0.001*N) | is.na(alpha_z)) #trait specific qc
  if(length(w)>0){
    alpha_z <- alpha_z[-w]
    transcripts <- transcripts[-w]
    starts <- starts[-w]
    chrs <- chrs[-w]
    tissueassign <- tissueassign[-w]
    sources <- sources[-w]
    X <- X[-w,]
    cells_per_ind <- cells_per_ind[-w]
  }
  N_tissuespecific <- as.numeric(table(tissueassign))
  a_transcripts <- table(transcripts)
  weights3 <- sapply(1:length(transcripts), function(x) as.numeric(a_transcripts)[match(transcripts[x], names(a_transcripts))]) #tissue redundancy weight to update genes that are regulated in more tissue-specific contexts
  weights2 <- sapply(1:nrow(X), function(x) sum(X[x,-tissueassign[x]])) + 1 #total co-regulation weight to prevent double counting of signal from co-regulated genes 
  expected_true_cish2_genes <- N_tissuespecific
  
  #### set up jackknife #### 
  a <- cbind(1:length(starts),as.numeric(chrs),as.numeric(starts),unlist(sapply(1:length(N_tissuespecific), function(x) rep(x,N_tissuespecific[x]))))
  a <- a[order(a[,2],a[,3],decreasing = F),]
  chunks <- 200 
  size_groups <- floor(length(starts)/chunks)
  #if (length(size_groups) == 0){
  #    counter = counter + 1
  #    bad_sstats[[counter]] = idx
  #  }
  size_group5 <- length(starts) - (chunks-1)*size_groups
  
  group_assignment <- c()
  
  for (j in 1:chunks){
    if(j == chunks){
      group_assignment <- c(group_assignment, rep(j,size_group5))
    }else{
      group_assignment <- c(group_assignment, rep(j,size_groups))
    }
  }
  
  a <- cbind(a, group_assignment)
  a <- a[order(a[,1], decreasing = F),]
  a <- cbind(a,transcripts)
  
  mean_chisq <- mean(alpha_z^2, na.rm = T)
  totcoreg <- rowSums(X)
  crude_h2_est <- (mean_chisq - 1)/(N*mean(totcoreg)) #should this be average co-reg score or average sum across tissues co-reg score? 
  weights1 <- (1 + N*crude_h2_est*totcoreg)^2 #bias corrected but not scaled X; heteroscedasticity weight based on s-ldsc
  
  variance_mat <- matrix(0,1+length(tissues),7) #tissue, h2ge_t, jackknife SE, P, FDR across tissues, proph2, proph2_se
  variance_mat[,1] <- c(tissues,"AllTissues")
  # here, I put tissueassign as cells.
  variance_mat[1:length(tissues),2] <- get_cov_alpha1alpha1_multitissue(sources,cells_per_ind,alpha_z,X,N,weights1,weights2,weights3) * expected_true_cish2_genes 
  variance_mat[(length(tissues)+1),2] <- sum(as.numeric(variance_mat[1:length(tissues),2]))
  
  jk <- matrix(0,nrow = chunks,ncol = length(tissues))
  jk_weights <- matrix(0,nrow = chunks,ncol = 1)
  jk_sum <- matrix(0,nrow = chunks,ncol = 1)
  
  for (chunk in 1:chunks){
	if (print_statements){
		print(paste0("processing jackknife chunk ",chunk," out of 200"))
	}
    remove_genes <- which(a[,5] == chunk)
    tab <- table(a[remove_genes,4]) #freq of tissues
    subtract_genes <- rep(0,length(tissues))
    m <- match(1:length(tissues), names(tab))
    w <- which(!is.na(m))
    subtract_genes[w] <- as.numeric(tab)[m[w]]
    N_tissuespecific_jk <- sapply(1:length(N_tissuespecific), function(x) N_tissuespecific[x] - subtract_genes[x])
    alpha_z_jk <- alpha_z[-remove_genes]
    tissueassign_jk <- tissueassign[-remove_genes]
    sources_jk <- sources[-remove_genes]
    cells_per_ind_jk <- cells_per_ind[-remove_genes]
    mean_chisq <- mean(alpha_z_jk^2, na.rm = T)
    X_jk <- X[-remove_genes,]
    totcoreg <- rowSums(X_jk)
    crude_h2_est <- (mean_chisq - 1)/(N*mean(totcoreg))
    weights1_jk <- (1 + N*crude_h2_est*totcoreg)^2
    weights2_jk <- weights2[-remove_genes] 
    weights3_jk <- weights3[-remove_genes] 
    jk[chunk,] <- get_cov_alpha1alpha1_multitissue(sources_jk,cells_per_ind_jk,alpha_z_jk,X_jk,N,weights1_jk,weights2_jk,weights3_jk) * N_tissuespecific_jk 
    jk_sum[chunk,1] <- sum(as.numeric(jk[chunk,]))
    jk_weights[chunk,1] <- sum(1/(weights1_jk*weights2_jk*weights3_jk))
  } 
  
  variance_mat[1:ncol(jk),3] <- sapply(1:ncol(jk), function(x) sqrt(wtd.var(jk[,x], jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0)))) 
  variance_mat[nrow(variance_mat),3] <- sqrt(wtd.var(jk_sum[,1],jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0)))
  
  variance_mat[,4] <- pnorm(q = 0, mean = as.numeric(variance_mat[,2]), sd = as.numeric(variance_mat[,3]))
  variance_mat[1:ncol(jk),5] <- p.adjust(as.numeric(variance_mat[1:ncol(jk),4]), method = "fdr")
  variance_mat[nrow(variance_mat),5] <- NA
  variance_mat[1:nrow(variance_mat),6] <- as.numeric(variance_mat[,2]) / h2g
  variance_mat[1:nrow(variance_mat),7] <- as.numeric(variance_mat[,3]) / h2g
                                      
  colnames(variance_mat) <- c("Tissue","h2_ge_t","h2_ge_t_se","NomP","FDRP","Proph2","Proph2_se")
 return(variance_mat) 
}