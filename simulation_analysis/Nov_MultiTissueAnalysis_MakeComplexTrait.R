sim_vegene_vesnp <- commandArgs(trailingOnly=TRUE)
s <- strsplit(sim_vegene_vesnp, split = "_")[[1]]
sim <- as.numeric(s[1])
ve <- s[2]
g <- s[3] #can be a list separated by commas
vsnp <- strsplit(g, split = ",")[[1]]
ngwas <- s[4]

library(data.table)

alpha_y1 <- as.matrix(fread("Nov_Y1_alpha_nCausalGenes100.txt.gz", header = F))
alpha_y2 <- as.matrix(fread("Nov_Y2_alpha_nCausalGenes100.txt.gz", header = F)) #found that this said Y1 on Dec 30

#GENE-TRAIT EFFECT: 
if(as.numeric(ve) != 0){
genetrait1 <- as.matrix(fread(paste0("genetraits/Nov_0.75_ggcoreg/",ngwas,"_vgene",ve,"/genetraitY1_sim",sim,".txt.gz"), header = F))
genetrait2 <- as.matrix(fread(paste0("genetraits/Nov_0.75_ggcoreg/",ngwas,"_vgene",ve,"/genetraitY2_sim",sim,".txt.gz"), header = F))

#add alpha here as multiplying factor.
w1 <- which(alpha_y1[,sim] != 0)
s <- t(sapply(1:nrow(genetrait1), function(x) genetrait1[x,]*alpha_y1[w1[x],sim]))
w1 <- which(alpha_y2[,sim] != 0)
s2 <- t(sapply(1:nrow(genetrait2), function(x) genetrait2[x,]*alpha_y2[w1[x],sim]))

ComplexTrait1 <- rowSums(t(s), na.rm = T)
ComplexTrait2 <- rowSums(t(s2), na.rm = T)
a <- cbind(ComplexTrait1,ComplexTrait2)
print(dim(a))
dir.create(paste0("ComplexTraits/",ngwas,"_Nov/"),showWarnings = F, recursive = T)
filename <- paste0("ComplexTraits/",ngwas,"_Nov/ComplexTraits_Nov_Group1_vgene",ve,"_sim",sim,".txt")
write.table(a, file = filename, row.names =F, col.names =F, quote=F, sep = "\t")
system(paste0("gzip -f ", filename))
}

#SNP-TRAIT EFFECT: 
if(all(as.numeric(vsnp) != 0)){
#if(length(vsnp) > 0){
for (i in 1:length(vsnp)){
genetrait1 <- as.matrix(fread(paste0("genetraits/Nov_0.75_ggcoreg/",ngwas,"_vsnp",vsnp[i],"/snptraitY1_sim",sim,".txt.gz"), header = F))
genetrait2 <- as.matrix(fread(paste0("genetraits/Nov_0.75_ggcoreg/",ngwas,"_vsnp",vsnp[i],"/snptraitY2_sim",sim,".txt.gz"), header = F))
ComplexTrait3 <- rowSums(t(genetrait1), na.rm = T) 
ComplexTrait4 <- rowSums(t(genetrait2), na.rm = T)
if(as.numeric(ve) != 0){a <- cbind(ComplexTrait1+ComplexTrait3,ComplexTrait2+ComplexTrait4)}else{a <- cbind(ComplexTrait3,ComplexTrait4)}
print(dim(a))
filename <- paste0("ComplexTraits/",ngwas,"_Nov/ComplexTraits_Nov_Group1_vgene",ve,"_vsnp",vsnp[i],"_sim",sim,".txt")
write.table(a, file = filename, row.names =F, col.names =F, quote=F, sep = "\t")
system(paste0("gzip -f ", filename))
}
}
