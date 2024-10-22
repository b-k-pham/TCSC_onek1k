sim_vegene_vesnp <- commandArgs(trailingOnly=TRUE)
s <- strsplit(sim_vegene_vesnp, split = "_")[[1]]
sim <- as.numeric(s[1])
ve_gene <- s[2]
ve_snp <- s[3]
ngwas <- s[4]
#ct <- s[4] #this will be mode; 1) +1 within module, 2) +1 outside module, 3) + one of each

library(data.table)
library(glmnet)
#tissues <- c(0:9)
tissues <- c(1:3)
ntiss <- length(tissues)
ngenes <- 1000
#samplesizes <- c(100,200,300,500,1000,1500)
samplesizes <- c(100)
traitfiles <- c()
dirs <- c() #outputs
if(as.numeric(ve_gene)!=0){traitfiles <- c(traitfiles, paste0("ComplexTraits_Nov_Group1_vgene",ve_gene,"_sim",sim,".txt.gz"))
dirs <- c(dirs,paste0("Nov_vgene",ve_gene,"_",ngwas,"/"))
}
if(as.numeric(ve_snp)!=0){traitfiles <- c(traitfiles, paste0("ComplexTraits_Nov_Group1_vgene",ve_gene,"_vsnp",ve_snp,"_sim",sim,".txt.gz"))
dirs <- c(dirs,paste0("Nov_vgene",ve_gene,"_vsnp",ve_snp,"_",ngwas,"/"))
}

system(paste0("mkdir TWASsumstats/Nov_vgene",ve_gene,"_",ngwas,"/"))
for (samp in 1:length(samplesizes)){
    if(as.numeric(ve_gene)!=0){system(paste0("mkdir -p TWASsumstats/Nov_vgene",ve_gene,"_",ngwas,"/",samplesizes[samp]))
    system(paste0("rm TWASsumstats/Nov_vgene",ve_gene,"_",ngwas,"/",samplesizes[samp],"/*_sim",sim,".txt"))
    system(paste0("rm TWASsumstats/Nov_vgene",ve_gene,"_",ngwas,"/",samplesizes[samp],"/*_sim",sim,".txt.gz"))
    }
    if(as.numeric(ve_snp)!=0){system(paste0("mkdir -p TWASsumstats/Nov_vgene",ve_gene,"_vsnp",ve_snp,"_",ngwas,"/",samplesizes[samp]))
    system(paste0("rm TWASsumstats/Nov_vgene",ve_gene,"_vsnp",ve_snp,"_",ngwas,"/",samplesizes[samp],"/*_sim",sim,".txt"))
    system(paste0("rm TWASsumstats/Nov_vgene",ve_gene,"_vsnp",ve_snp,"_",ngwas,"/",samplesizes[samp],"/*_sim",sim,".txt.gz"))
    }
}


#for multiple causal tissues 
#modes <- c("2CausalTissues_withinmodule","2CausalTissues_outsidemodule","3CausalTissues")
#ComplexTrait_Y1 <- ComplexTraits$V1
#ComplexTrait_Y2 <- ComplexTraits$V2
#n <- length(ComplexTrait_Y1)
#print(n)
gwas_geno <- as.matrix(fread(paste0("Simulated_GWAS_Cohort_",ngwas,".txt.gz"), header = F)) #need to index into this one below to get the right snps.  
genenum <- sort(rep(1:ngenes,5))

for (traitname in 1:length(traitfiles)){
    print(traitfiles[traitname])
    ComplexTraits <- fread(paste0("ComplexTraits/",ngwas,"_Nov/",traitfiles[traitname]),header = F) 
    ComplexTrait_Y1 <- ComplexTraits$V1
    ComplexTrait_Y2 <- ComplexTraits$V2
    n <- length(ComplexTrait_Y1)

    for (ti in tissues){ 
        print(paste0("working on tissue ",ti))
        eqtl <- fread(paste0("eQTLs_1000sims_Nov/Nov_sims_eQTLeffectsizes_0.75rg_SNPpolygen5_tissue",ti,"_sim",sim,".txt.gz"), header = F)
        coef <- as.matrix(read.table(paste0("weights/Nov_0.75_ggcoreg_coef/coef_Nov_Group1_Tissue",ti,"_sim",sim,".txt.gz"), header = F, fill = T, col.names = paste0("V",seq_len(45))))
        weights <- as.matrix(fread(paste0("weights/Nov_0.75_ggcoreg/h2_r2_Nov_Group1_Tissue",ti,"_sim",sim,".txt.gz"), header = F)) #need this to be 6000 rows. 
        #they are different dimensions; only weights has chunk, coef is only for those that are cis heritable

        for (samp in 1:length(samplesizes)){
            print(paste0("working on sample size", samplesizes[samp]))
            w1 <- which(weights[,1] == samplesizes[samp] & weights[,4] < 0.01 & weights[,3] > 0) #[nqtl, gene, h2g, hsq_p, r2all, r2train]
            weights_samp <- weights[w1,] #need this to get the chunk right. 
            w2 <- which(coef[,1] == samplesizes[samp]) #these are the cis heritable genes 
            coef_samp <- coef[w2,] #this matrix has all different numbers of fields
            if(length(w1) == length(w2)){print("True")}else{print("Problem")}
  
            for (gene in 1:ngenes){ #should be ok to put here because separate files for samp and tissue 

                w <- which(weights_samp[,2] == gene)
                if(length(w) > 0){ 
                    minsnp <- as.numeric(eqtl[which(genenum == gene)[1],3])
                    maxsnp <- as.numeric(eqtl[which(genenum == gene)[1],4])
                    index <- c(minsnp:maxsnp)
                    print("there are cis heritable genes")
                    coef_samp_chunk <- coef_samp[w,]
                    coef_samp_chunk <- coef_samp_chunk[which(!is.na(coef_samp_chunk))] #added this

                    #run linear models 
                    for (tr in 1:2){
                        trait <- get(paste0("ComplexTrait_Y",tr))
                        marginal_twas <- matrix(0,nrow = length(w), ncol = 3)
                        for (j in 1:length(w)){ #genes
                            if(length(w) == 1){
                                predExp <- gwas_geno[,index] %*% matrix(coef_samp_chunk[-1], ncol = 1) #if vector 
                            }else{ #when does this happen? 
                                print("something weird is happening")
                                predExp <- gwas_geno[,index] %*% matrix(coef_samp_chunk[j,-1], ncol = 1) #if matrix 
                            }
                            tryCatch({
                                lmod <- summary(lm(trait ~ predExp))
                                marginal_twas[j,1:3] <- coef(lmod)[2,1:3]
                            },error=function(cond){message(paste("Invalid regression ", j))})
                        } #over genes
						dir.create(paste0("TWASsumstats/",dirs[traitname],"/",samplesizes[samp],"/"),showWarnings = F,recursive = T)
                        write.table(marginal_twas, file = paste0("TWASsumstats/",dirs[traitname],"/",samplesizes[samp],"/TWAS_Group1_T",ti,"_Trait",tr,"_sim",sim,".txt"), row.names = F, col.names = F, sep = "\t", quote = F, append = T)
                     } #traits
                } #if there are any cis h2 genes for this tissue, sample pair
            } #over genes (new)
        } #sample size
    } #tissues
} #types of traits 

for (traitname in 1:length(traitfiles)){
    for (samp in 1:length(samplesizes)){
        system(paste0("gzip TWASsumstats/",dirs[traitname],"/",samplesizes[samp],"/*_sim",sim,".txt"))
    }
}

