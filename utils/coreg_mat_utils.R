library(tidyverse)
# Function to load the required data
load_data <- function(tissues, sources, i) {
  tissue = tissues[i]
  return(list(df1 = get(load(paste0("TCSC/predicted_expression/onek1k/", tissue, "_onek1k.RData"))), transcripts1 = get(load(paste0('TCSC/predicted_expression/onek1k/', tissue, '_onek1k_transcripts.RData'))), tissue = tissue,source = sources[i]))
}

# Function to merge gene_ref and gtex ref
common_generef_gtex <- function(gene_ref,gtex){
  colnames(gene_ref)[1] = 'chr'
  colnames(gene_ref)[2] = 'start'
  colnames(gene_ref)[3] = 'end'
  
  m <- match(gene_ref$symbol, gtex$V4)
  
  #length(which(is.na(m))) #1,415 genes are missing
  genetype <- gtex$V8[m]
  gene_name <- gtex$V4[m]
  transcript_name <- gtex$V7[m]
  gene_ref <- cbind(gene_ref,genetype,gene_name,transcript_name)
  
  return(gene_ref)
}

# Function to filter protein-coding genes
filter_protein_coding_genes <- function(transcript_key, df1_focal, gene_ref, source) {
  if (source == '1k1k') {
    m <- match(transcript_key, gene_ref$gene_name)
  } else {
    m <- match(transcript_key, gene_ref$transcript_name)
  }
  
  focal_tissue_gene_ref <- gene_ref[m,]
  keepgenes_focal <- which(focal_tissue_gene_ref[,5] == "protein_coding")
  transcript_key <- transcript_key[keepgenes_focal]
  df1_focal <- df1_focal[, keepgenes_focal]
  focal_tissue_gene_ref <- focal_tissue_gene_ref[keepgenes_focal,]
  
  return(list(transcript_key = transcript_key, df1_focal = df1_focal, focal_tissue_gene_ref = focal_tissue_gene_ref))
}

# Function to calculate co-regulation scores for one tissue
calculate_coregulation <- function(i, j,focal_filtered,sources,gene_ref, tissues, ciswindow) {
  source_j <- sources[j]
  if (i == j) {
    acc <- fread(paste0("TCSC/acc_dfs/", tissues[i], "_onek1k_accs.txt"), header = F)
    acc <- acc[match(focal_filtered$transcript_key, acc$V1)]$V2
    #genegenecoreg_pc <- sapply(1:ncol(focal_filtered$df1_focal), function(k) {
	genegenecoreg_pc <- c() #vector of length g in focal tissue 
	for (k in 1:ncol(focal_filtered$df1_focal)){
      gene_k_chr <- focal_filtered$focal_tissue_gene_ref$chr[k]
      gene_k_pos <- focal_filtered$focal_tissue_gene_ref$start[k]
      cis_genes <- focal_filtered$focal_tissue_gene_ref[which(focal_filtered$focal_tissue_gene_ref$chr == gene_k_chr & focal_filtered$focal_tissue_gene_ref$start > gene_k_pos - ciswindow & focal_filtered$focal_tissue_gene_ref$start < gene_k_pos + ciswindow),] # equivalent of making focal_tissue_gene_ref_k_cislocus
	  cis_genes_pc <- match(cis_genes$symbol,focal_filtered$focal_tissue_gene_ref$symbol)
      if(length(cis_genes_pc) == 0){cor_est_pc <- 0}else{cor_est_pc <- sapply(cis_genes_pc, function(x) cor(as.numeric(focal_filtered$df1_focal[,k]),as.numeric(focal_filtered$df1_focal[,x])))}
      genegenecoreg_pc  <- c(genegenecoreg_pc,sum(cor_est_pc^2) - 1 + acc[k])
    }#)
  } else {
    other_data <- load_data(tissues, sources, j)
    other_filtered <- filter_protein_coding_genes(other_data$transcripts1, other_data$df1, gene_ref, source_j)
    #genegenecoreg_pc <- sapply(1:ncol(focal_filtered$df1_focal), function(k) {
	genegenecoreg_pc <- c() #vector of length g in focal tissue 
	for (k in 1:ncol(focal_filtered$df1_focal)){
      gene_k_chr <- focal_filtered$focal_tissue_gene_ref$chr[k]
      gene_k_pos <- focal_filtered$focal_tissue_gene_ref$start[k]
      cis_genes <- other_filtered$focal_tissue_gene_ref[which(other_filtered$focal_tissue_gene_ref$chr == gene_k_chr & other_filtered$focal_tissue_gene_ref$start > gene_k_pos - ciswindow & other_filtered$focal_tissue_gene_ref$start < gene_k_pos + ciswindow),] # equivalent of making focal_tissue_gene_ref_k_cislocus
	  cis_genes_pc <- match(cis_genes$symbol,other_filtered$focal_tissue_gene_ref$symbol)
      if(length(cis_genes_pc) == 0){cor_est_pc <- 0}else{cor_est_pc <- sapply(cis_genes_pc, function(x) cor(as.numeric(focal_filtered$df1_focal[,k]),as.numeric(other_filtered$df1_focal[,x])))}
      genegenecoreg_pc <- c(genegenecoreg_pc,sum(cor_est_pc^2))
    }#)
  }
  return(genegenecoreg_pc)
}

# Main function to calculate co-regulation score matrix X
calculate_coregulation_matrix <- function(ntissues, tissues, sources, gene_ref, ciswindow) {
  X <- matrix(0, nrow = 2, ncol = ntissues)
  
  for (i in 1:ntissues) {
    print(paste0("Working on tissue ", i, " out of ", ntissues))
	focal_data <- load_data(tissues,sources,i)
	source_i = sources[i]
	focal_filtered <- filter_protein_coding_genes(focal_data$transcripts1, focal_data$df1, gene_ref, source_i)
	#return(list(transcript_key = transcript_key, df1_focal = df1_focal, focal_tissue_gene_ref = focal_tissue_gene_ref))
    coregmat <- matrix(0, nrow = length(focal_filtered$transcript_key), ncol = ntissues)
    for (j in 1:ntissues) {
      #coregmat[, j] <- calculate_coregulation(i, j, filtered$df1_focal, filtered$transcript_key, filtered$focal_tissue_gene_ref, sources, tissues, ciswindow)
	  coregmat[, j] <- calculate_coregulation(i, j,focal_filtered, sources,gene_ref, tissues, ciswindow)
    }
    
    X <- rbind(X, coregmat)
  }
  
  return(X[3:nrow(X), ])
}


# Function to format coreg matrix (get X, transcripts, starts, chrs, tissueassign, and sources per tissue.)

process_tissue <- function(i,tissues,sources,gene_ref){
  tissue = tissues[i]
  source = sources[i]
  data = get(load(paste0("TCSC/predicted_expression/onek1k/",tissue,"_onek1k.RData"))) # for now, this doesn't do anything but does it matter if number of individuals are small?
  transcript_key = get(load(paste0('TCSC/predicted_expression/onek1k/',tissue,'_onek1k_transcripts.RData')))
  keep = transcript_key
  
  w <- which(transcript_key %in% keep)
  transcript_key <- transcript_key[w]
  if (source == '1k1k'){
    m <- match(transcript_key,gene_ref$gene_name) # do genes instead of ENSG transcripts.
  } else {
    m <- match(transcript_key,gene_ref$transcript_name) # do ENSG transcripts instead
  }
  
  genetype <- gene_ref$genetype[m]
  transcript_key <- transcript_key[which(genetype == "protein_coding")]
  
  #transcripts <- c(transcripts,transcript_key)
  #tissueassign <- c(tissueassign, rep(i,length(transcript_key)))
  #sources <- c(sources,rep(source,length(transcript_key)))
  
  if (source == '1k1k'){
    temp = gene_ref %>% filter(gene_name %in% transcript_key) %>% filter(genetype == "protein_coding") %>% select(chr,start,end,gene_name,transcript_name) %>% arrange(match(gene_name,transcript_key))
  } else {
    temp = gene_ref %>% filter(transcript_name %in% transcript_key) %>% filter(genetype == "protein_coding") %>% select(chr,start,end,gene_name,transcript_name) %>% arrange(match(transcript_name,transcript_key))
  }
  
  
  #starts <- c(starts,temp$start)
  #chrs <- c(chrs,temp$chr)
  
  res = list(transcripts = transcript_key,tissueassign_labels = rep(tissue,length(transcript_key)),tissueassign = rep(i,length(transcript_key)),sources = rep(source,length(transcript_key)),starts = temp$start,chrs = temp$chr)
  
  return(res)
}