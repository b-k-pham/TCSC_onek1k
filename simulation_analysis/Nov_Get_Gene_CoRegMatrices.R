sim <- as.numeric(commandArgs(trailingOnly=TRUE))
library(data.table)
#samplesizes <- c(100,200,300,500,1000,1500)
samplesizes <- c(100)
geneannot <- fread("1000genes_Nov.txt", header = T)
ngenes <- 1000
#tissues <- c(0:9)
tissues <- c(1:3)
othertissues <- tissues

for (samp in 1:length(samplesizes)){

    for (focal in tissues){
        focal_T <- as.matrix(fread(paste0("predExp/Nov_0.75_ggcoreg/",samplesizes[samp],"/predExp_Nov_Group1_sim",sim,"_Tissue",focal,".txt.gz"), header = F)) 
        h2 <- as.matrix(fread(paste0("weights/Nov_0.75_ggcoreg/h2_r2_Nov_Group1_Tissue",focal,"_sim",sim,".txt.gz"), header = F))
        w <- which(h2[,4] < 0.01 & h2[,3] > 0 & h2[,1] == samplesizes[samp]) #[nqtl, gene, h2g, hsq_p, r2all, r2train] 
        h2_focal <- h2[w,] 

        ww <- which(h2[,1] == samplesizes[samp])
        index_for_genes_focal <- which(h2[ww,4] < 0.01 & h2[ww,3] > 0) #need this to be out of 1000 

        PartialCoReg <- matrix(0,nrow(focal_T),length(tissues))
        for (oth in othertissues){
            print(oth)
            if(oth == focal){
            genegenecoreg <- c()
            for (k in 1:nrow(focal_T)){
                start <- as.numeric(geneannot$START[h2_focal[k,2]])
                stop <- as.numeric(geneannot$STOP[h2_focal[k,2]]) #col 2 is the gene index. 
                cis_genes <- which(abs(stop-geneannot$START[index_for_genes_focal]) < 1000000 | abs(start-geneannot$STOP[index_for_genes_focal]) < 1000000)
                if(length(cis_genes) == 0){s_corr <- 0}else{s_corr <- sapply(cis_genes,function(x) cor(as.numeric(focal_T[k,]),as.numeric(focal_T[x,])))}
                coreg_before_sum <- s_corr^2 
                genegenecoreg <- c(genegenecoreg, sum(coreg_before_sum, na.rm = T))
            }
            }else{
            oth_T <- as.matrix(fread(paste0("predExp/Nov_0.75_ggcoreg/",samplesizes[samp],"/predExp_Nov_Group1_sim",sim,"_Tissue",oth,".txt.gz"), header = F)) 
            h2 <- as.matrix(fread(paste0("weights/Nov_0.75_ggcoreg/h2_r2_Nov_Group1_Tissue",oth,"_sim",sim,".txt.gz"), header = F))
            w <- which(h2[,4] < 0.01 & h2[,3] > 0 & h2[,1] == samplesizes[samp])
            h2_oth <- h2[w,]

            ww <- which(h2[,1] == samplesizes[samp])
            index_for_genes <- which(h2[ww,4] < 0.01 & h2[ww,3] > 0)

            genegenecoreg <- c() #vector of length g in focal tissue 
            for (k in 1:nrow(focal_T)){
                start <- as.numeric(geneannot$START[h2_focal[k,2]])
                stop <- as.numeric(geneannot$STOP[h2_focal[k,2]]) #col 2 is the gene index. 
                cis_genes <- which(abs(stop-geneannot$START[index_for_genes]) < 1000000 | abs(start-geneannot$STOP[index_for_genes]) < 1000000) 
                if(length(cis_genes) == 0){s_corr <- 0}else{s_corr <- sapply(cis_genes,function(x) cor(as.numeric(focal_T[k,]),as.numeric(oth_T[x,])))}
                coreg_before_sum <- s_corr^2
                genegenecoreg <- c(genegenecoreg, sum(coreg_before_sum, na.rm = T))
            }
            } #else: oth != focal

           #store genegenecoreg in matrix 
           PartialCoReg[,oth] <- genegenecoreg
       } #other tissues loop 
    assign(paste0("CoRegMat_",focal),PartialCoReg)
    } #focal tissues loop 

    FullCoReg <- data.frame()
    for (focal in tissues){
    a <- get(paste0("CoRegMat_",focal))
    FullCoReg <- rbind(FullCoReg,a)
    }
	dir.create(paste0("TissueCoReg/Nov_0.75_ggcoreg/",samplesizes[samp],"/"),showWarnings = F, recursive = T)
    write.table(FullCoReg, file = paste0("TissueCoReg/Nov_0.75_ggcoreg/",samplesizes[samp],"/GeneCoReg_PredPred_Nov_Group1_nqtl",samplesizes[samp],"_sim",sim,".txt"), row.names = F, col.names = F, sep = "\t", quote = F)
    system(paste0("gzip -f TissueCoReg/Nov_0.75_ggcoreg/",samplesizes[samp],"/GeneCoReg_PredPred_Nov_Group1_nqtl",samplesizes[samp],"_sim",sim,".txt"))
}
