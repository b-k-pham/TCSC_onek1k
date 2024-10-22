sim_vegene_vesnp <- commandArgs(trailingOnly=TRUE)
s <- strsplit(sim_vegene_vesnp, split = "_")[[1]]
sim <- as.numeric(s[1])
ve_gene <- s[2]
ve_snp <- s[3]
ngwas <- as.character(s[4])

library(data.table)
library(Hmisc) 
tissues <- 0:2
#tissues <- 1:3
#samplesizes <- c(100, 200, 300, 500, 1000, 1500)
samplesizes <- c(100)
ntiss <- length(tissues)
ntraits <- 2 
nsamp <- length(samplesizes)
ngenes <- 1000
#ngwas <- 10000
rg <- 0.75

## set up jackknife across ntiss tissues ## 
chunks <- 20
starts <- c(1:ngenes)
a <- cbind(1:length(starts),as.numeric(starts)) #col1 = gene index from 1 to genes *tissues, col2 = gene start position, col3 = group assignment, col4 = tissue
size_groups <- floor(length(starts)/chunks)
size_group5 <- length(starts) - (chunks-1)*size_groups
group_assignment <- c()
for (j in 1:chunks){
    if(j == chunks){
    group_assignment <- c(group_assignment, rep(j,size_group5))
    }else{group_assignment <- c(group_assignment, rep(j,size_groups))}
}
a <- cbind(a, group_assignment)
a <- a[order(a[,1], decreasing = F),]
aa <- do.call("rbind", replicate(ntiss, a, simplify = FALSE))
aa[,1] <- 1:(ntiss * ngenes)
aa <- cbind(aa, sort(rep(1:ntiss,ngenes))) #already the case that all of the same transcript per tissue is in the same block. 
#double check same number of genes in each chunk. table(aa[,3]) confirmed 

## functions 
source("TCSC_fun_modJuly.R")

## load data 
mastermat_h2 <- matrix(0,nrow = ngenes*length(samplesizes), ncol = ntiss+1) #2085 originally, all chunks ready for tissues 0-11 
mastermat_h2p <- mastermat_h2
mastermat_R2CV <- mastermat_h2
for (ti in 1:ntiss){
    print(ti)
    h2 <- as.matrix(fread(paste0("weights/Nov_",rg,"_ggcoreg/h2_r2_Nov_Group1_Tissue",ti,"_sim",sim,".txt.gz"), header = F)) #[nqtl, gene, h2g, hsq_p, r2all, r2train]
    #h2 <- h2[,-c(1,3)] #remove chunk and gene columns, left: nqtl, h2g, hsq_p, r2all, r2train #update this #[nqtl, gene, h2g, hsq_p, r2all, r2train]
    if(ti == 1){
        h2est <- h2[,c(1,3)]
        h2p <- h2[,c(1,4)]
        R2CV <- h2[,c(1,5)]
        mastermat_h2[,c(1:2)] <- h2est
        mastermat_h2p[,c(1:2)] <- h2p
        mastermat_R2CV[,c(1:2)] <- R2CV
    }else{
        h2est <- h2[,3]
        h2p <- h2[,4]
        R2CV <- h2[,5]
        mastermat_h2[,ti+1] <- h2est
        mastermat_h2p[,ti+1] <- h2p
        mastermat_R2CV[,ti+1] <- R2CV
    }
}

for (samp in 1:nsamp){ 
    w <- which(mastermat_h2[,1] == samplesizes[samp]) #h2 estimate
    gctah2_samp <- mastermat_h2[w,-1]
    w <- which(mastermat_h2p[,1] == samplesizes[samp]) #h2 p 
    gctah2p_samp <- mastermat_h2p[w,-1]
    w <- which(mastermat_R2CV[,1] == samplesizes[samp]) #r2cv 
    r2cv_samp <- mastermat_R2CV[w,-1]
    genenames <- unlist(sapply(1:ncol(gctah2p_samp), function(x) which(gctah2p_samp[,x] < 0.01 & gctah2_samp[,x] > 0))) #added second term 4/12
    assign(paste0("genenames_",samp),genenames)
    X_PP_gene <- as.matrix(fread(paste0("TissueCoReg/Nov_",rg,"_ggcoreg/",samplesizes[samp],"/GeneCoReg_PredPred_Nov_Group1_nqtl",samplesizes[samp],"_sim",sim,".txt.gz"), header= F)) 
    assign(paste0("X_PP_gene_",samp),X_PP_gene)
    allw_h2genes <- c()
    N_tissue_specific <- c()
    N_tissue_specific_unique <- c()
    count <- 0
    for (ti in tissues){
        w_h2genes <- which(gctah2p_samp[,ti+1] < 0.01 & gctah2_samp[,ti+1] > 0) #added second term 4/12
        allw_h2genes <- c(allw_h2genes, w_h2genes+count)
        N_tissue_specific_unique <- c(N_tissue_specific_unique, length(w_h2genes))
        N_tissue_specific <- c(N_tissue_specific, rep(length(w_h2genes),length(w_h2genes)))
        count <- count + ngenes
    }
    assign(paste0("N_tissue_specific_",samp),N_tissue_specific)
    assign(paste0("N_tissue_specific_unique_",samp),N_tissue_specific_unique)
    assign(paste0("allw_h2genes_",samp),allw_h2genes)
    tissuelabel <- aa[allw_h2genes,4]
    X_PP_gene_bc <- X_PP_gene
    for (ti in tissues){
        w <- which(tissuelabel == ti + 1)
        w_h2genes <- which(gctah2p_samp[,ti+1] < 0.01 & gctah2_samp[,ti+1] > 0) #added second term 4/12
        r2 <- r2cv_samp[w_h2genes,ti+1]
        h2 <- gctah2_samp[w_h2genes,ti+1]
        r2[r2 < 0] <- 0
        r2[r2 > 1] <- 1
        h2[h2 > 1] <- 1
        h2[h2 < 0] <- 0
        ww <- which(r2 > h2)
        if(length(ww) > 0){r2[ww] <- h2[ww]}
        acc <- r2/h2
        acc[acc == "Inf" | acc == "-Inf"] <- 0
        X_PP_gene_bc[w,ti+1] <- X_PP_gene[w,ti+1] - 1 + acc
    }
    assign(paste0("X_PP_gene_bc_",samp),X_PP_gene_bc)
    assign(paste0("X_PP_gene",samp),X_PP_gene)
}

twas_sumstat_dir <- c()
if(as.numeric(ve_gene)!=0){twas_sumstat_dir <- c(twas_sumstat_dir, paste0("vgene",ve_gene))}
if(as.numeric(ve_snp)!=0){twas_sumstat_dir <- c(twas_sumstat_dir, paste0("vgene",ve_gene,"_vsnp",ve_snp))} 

for (dir in 1:length(twas_sumstat_dir)){ 
## single trait analysis
    varalphamat <- matrix(0,ntraits*nsamp,3+ntiss*2+1) 
    varcount <- 1
    for (tr in 1:ntraits){
        for (samp in 1:nsamp){
            X_PP_h2genes <- get(paste0("X_PP_gene_bc_",samp))
            X_PP_gene_uncorrected <- get(paste0("X_PP_gene",samp))
            genenames <- get(paste0("genenames_",samp)) #names, e.g. 1:999 , therefore table(genenames) = this gene is expressed in how many tissues? 
            N_tissue_specific <- get(paste0("N_tissue_specific_",samp))
            N_tissue_specific_unique <- get(paste0("N_tissue_specific_unique_",samp))
            allw_h2genes <- get(paste0("allw_h2genes_",samp))
            alpha_z <- c()
            for (ti in tissues){
                alpha_z <- c(alpha_z,fread(paste0("TWASsumstats/Nov_",twas_sumstat_dir[dir],"_",ngwas,"/",samplesizes[samp],"/TWAS_Group1_T",ti+1,"_Trait",tr,"_sim",sim,".txt.gz"), header = F)$V3)
            }
            varalphamat <- run_TCSC_2steps(alpha_z, varalphamat, varcount, ngenes, N_tissue_specific, as.numeric(ngwas), tr, ve_gene, samp, X_PP_h2genes, genenames,X_PP_gene_uncorrected,ntiss, allw_h2genes, N_tissue_specific_unique) #allw_h2genes after ntiss?? 
            varcount <- varcount + 1
            print(varcount)
        }
    }
	dir.create('Group1res_Nov/',showWarnings = F, recursive = T)
    write.table(varalphamat, file = paste0("Group1res_Nov/TCSC_Group1_variance_T10_sim",sim,"_010323_",twas_sumstat_dir[dir],"_",ngwas,".txt"), row.names = F, col.names = F, quote = F, sep = "\t")


    gencorralphamat <- matrix(0,nsamp,ntiss*2+1)
    gencorrcount <- 1
    for (samp in 1:length(samplesizes)){
        allw_h2genes <- get(paste0("allw_h2genes_",samp))
        X_PP_h2genes <- get(paste0("X_PP_gene_bc_",samp))
        X_PP_gene_uncorrected <- get(paste0("X_PP_gene",samp))
        genenames <- get(paste0("genenames_",samp)) #names, e.g. 1:999 
        N_tissue_specific <- get(paste0("N_tissue_specific_",samp))
        N_tissue_specific_unique <- get(paste0("N_tissue_specific_unique_",samp))
        count <- 0
        alpha_z1 <- c()
        alpha_z2 <- c()
        for (ti in tissues){
            for (tr in 1:2){
                alpha_z <- fread(paste0("TWASsumstats/Nov_",twas_sumstat_dir[dir],"_",ngwas,"/",samplesizes[samp],"/TWAS_Group1_T",ti+1,"_Trait",tr,"_sim",sim,".txt.gz"), header = F)$V3
                alpha_z_tr <- get(paste0("alpha_z",tr))
                alpha_z_tr <- c(alpha_z_tr,alpha_z)
                assign(paste0("alpha_z",tr),alpha_z_tr)
            }     
            count <- count + ngenes
         }
         gencorralphamat <- run_TCSC_Rg(alpha_z1,alpha_z2,gencorralphamat,gencorrcount,ngenes, as.numeric(ngwas), samp, allw_h2genes, X_PP_h2genes, genenames, aa, X_PP_gene_uncorrected, ntiss, N_tissue_specific_unique)
         gencorrcount <- gencorrcount + 1
         print(gencorrcount)
    }
    write.table(gencorralphamat, file = paste0("Group1res_Nov/TCSC_Group1_covariance_T10_sim",sim,"_010323_",twas_sumstat_dir[dir],"_",ngwas,".txt"), row.names = F, col.names = F, quote = F, sep = "\t")
} #over trait files (ve_gene vs ve_snp)
