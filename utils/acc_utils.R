library(tidyverse)



gtex_get_acc_dfs <- function(ct){
  h2_acc_link = paste0('TCSC/weights/accuracies/N320_smalltoo_accuracies_cish2_fromall_',ct,'.txt.gz')
  keep_link = paste0('TCSC/weights/heritablegenes/N320/TranscriptsIn',ct,'Model_keep.txt')
  keeps = fread(keep_link, header = F)$V1
  h2_acc_keeps = fread(h2_acc_link, header = F) %>% filter(V1 %in% keeps)
  return(h2_acc_keeps)
}


gtex_get_acc <-function(ct,transcript){
  temp = gtex_get_acc_dfs(ct) %>% filter(V1 == transcript)
  out = temp$V2
  if (length(out) == 0){
    out = NA
  }
  return(out)
}



get_acc <- function(weight){
  load(weight)
  performance = cv.performance
  performance = performance[,which(is.na(colSums(performance)) == F)] # check which column is NA, only grab columns that are all non NA.
  h2 = hsq[1]
  i = which(performance[2,] == min(performance[2,]))
  if (length(i) > 1){ # if there are multiple with the same p-value, just take the first one
    i = 1
  }
  out = performance[1,i]/h2 # R2.cv/h2 #crossvalidation r2 of model over true h2.
  # If NA, then shouldn't use the model in the first place.
  if (is.na(out) == T){ # if either the R2.cv is NA or h2 is NA, flag as bad with NA
    out = NA
  }
  if (performance[1,i] < 0){ # if R2.cv is < 0, then flag as bad with NA
    out = NA
  }
  return(out)
}


#args = commandArgs(trailingOnly=TRUE)


make_acc_df <- function(ct,source = 'onek1k'){
  if (source == '1k1k'){
    source = 'onek1k'
  }
  if (source == 'onek1k'){
    dir_name = paste0('prep_TCSC/RA_only/',ct,'_5pcs_all/')
  } else {
    #dir_name = paste0() # must be gtex
    dir_name = paste0('E:/PhD/manual_pipeline/tcsc_original_weights/v8_320EUR_meta_only/',ct,'/') # must be gtex
  }
  
  
  
  weights = list.files(dir_name,full.names = T)
  #
  weight_accs = list()
  gene_names = list()
  valid_files = list()
  counter = 0
  for (i in 1:length(weights)){ # for each CT, load weight file. Compute accuracy of said weight.
    if (source == 'onek1k'){
      gene_names[[i]] = str_split(str_split(weights[i],'_1Mwindow.wgt.RDat')[[1]][1],'[.]')[[1]][2]
      weight_acc = get_acc(weights[i])
      weight_accs[[i]] = weight_acc
    } else{
      gene_name = gsub("^\\.|\\.$", "", str_split(str_split(weights[i],ct)[[1]][3],'wgt')[[1]][1])
      gene_names[[i]] = gene_name
      #weight_acc = gtex_get_acc(ct,gene_name)
    }
    
  }
  
  gene_names = unlist(gene_names)
  weight_accs = unlist(weight_accs)
  if (source == 'onek1k'){
    weight_accs = unlist(weight_accs)
  } else{
    gtex_acc_df = data.frame(V1 = gene_names) %>% left_join(gtex_get_acc_dfs(ct)) %>% filter(V1 %in% gene_names)
    weight_accs = gtex_acc_df$V2
  }
  
  df = na.omit(data.frame(gene_names = gene_names,weight_accs = weight_accs,file_paths = weights))
  df = df %>% mutate(weight_accs = ifelse(weight_accs > 1,1,weight_accs))
  
  
  # get valid score.profile files and filter based off that
  if (source == 'onek1k'){
    scores_dir = paste0('prep_TCSC/score_files/',ct)
  } else{
    scores_dir = paste0('prep_TCSC/score_files_gtex_meta/',ct)
  }
  
  scores = list.files(scores_dir,pattern = '.score.profile')
  scores = unlist(str_split(scores,paste0(ct,'.')))
  scores = scores[scores != '']
  scores = str_split(scores,'.score')
  genes = list()
  for (i in 1:length(scores)){
    genes[[i]] = scores[[i]][1]
  }
  
  genes = unique(unlist(genes))
  
  df1 = df %>% filter(gene_names %in% genes) %>% mutate(score_profile_paths = paste0(scores_dir,'/',ct,'.',gene_names,'.score.profile'))
  colnames(df1)[1] = 'ID'
  write.table(df1,paste0('TCSC/acc_dfs/',ct,'_onek1k_accs.txt'),quote = F, row.names = F, col.names = F,sep = '\t')
  return(df1)
}
