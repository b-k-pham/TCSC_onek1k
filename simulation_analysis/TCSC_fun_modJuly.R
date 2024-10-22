

#these functions account for recent updates to method 
#covariates (# cis h2 genes per tissue, # genes expressed per tissue) 
#solving for h2ge_t' / Gt' for each tissue t, then multiply by Gt' 
#co-reg bug fix to not utilize co-regulation with unobserved cis h2 genes. 

library(nnls)

get_cov_alpha1alpha1_multitissue_withcrude <- function(alpha_z,gtot,n,x,genenames,x_uncorrected){
y <- (alpha_z^2)*(gtot/n) #modified mar 8 from -1; mod Jan 31 commented out 
weights <- get_weights_0827(x,genenames,alpha_z,n,x_uncorrected)
mod <- summary(glm(y ~ x, weights = weights)) #modified mar 8 from + covars
cov_b1b1 <- coef(mod)[-1,1]
cov_b1b1 <- cov_b1b1[1:ntiss]
cov_b1b1 <- cov_b1b1
return(cov_b1b1)
}

get_cov_alpha1alpha1_multitissue_withcrude_SigCish2Genes <- function(alpha_z,gtot,n,x,genenames,x_uncorrected,gt){
y <- (alpha_z^2)*(1/n) 
weights <- get_weights_0827(x,genenames,alpha_z,n,x_uncorrected)
mod <- summary(glm(y ~ x, weights = weights)) #modified mar 8 from + covars
cov_b1b1 <- coef(mod)[-1,1]
cov_b1b1 <- cov_b1b1[1:ntiss]
cov_b1b1 <- cov_b1b1*gt
return(cov_b1b1)
}

get_cov_alpha1alpha2_multitissue_withcrude_SigCish2Genes <- function(alpha_z1,alpha_z2,gtot,n,x,genenames,x_uncorrected,Gt){

mean_coreg <- mean(rowSums(x), na.rm = T) #modified mar 8 from just x
mean_chisq <- mean(alpha_z1^2, na.rm = T)
crude_h2_1 <- (mean_chisq - 1)/(n*mean_coreg) #per-SNP heritability = quantity h2/M, therefore no M in equation; mean tau.
mean_chisq <- mean(alpha_z2^2, na.rm = T)
crude_h2_2 <- (mean_chisq - 1)/(n*mean_coreg)

y <- (alpha_z1*alpha_z2)*(1/n)
#1.rho_pheno*NS=0, rho_g_est = the below
rho_g_est <- sum(alpha_z1*alpha_z2)/n #the first regression uses weights with rho_g_est and rho_pheno = 0. 
weights <- get_weights_0118(x,genenames,n,crude_h2_1,crude_h2_2,rho_g_est,0,gtot,x_uncorrected)
if(is.null(nrow(x))){rs <- x}else{rs <- rowSums(x)}
mod <- summary(glm(y ~ rs, weights = weights)) #modified mar 8 from + covars
#get estimates for rho_pheno*Ns and rho_g
intercept <- coef(mod)[1,1] #rho_pheno*Ns/sqrt(N1N2)
rhop_NS <- intercept*n #cannot be greater in magnitude than n*1 and it is if we do multi variate regression. using rowsums, then we get -777, which means pheno corr = -0.78, but we know it is 0.03. therefore something is wrong. we get 3K > 1K. 
rho_g_est2 <- coef(mod)[2,1]*sum(Gt) #need a single estimate, not one per tissue, so maybe first regression should not be multi variate. or just take the estimate for the causal tissue? 

weights <- get_weights_0118(x,genenames,n,crude_h2_1,crude_h2_2,rho_g_est2,rhop_NS,gtot,x_uncorrected)
mod <- summary(glm(y ~ x, weights = weights)) ##modified mar 8 from + covars
cov_b1b1 <- coef(mod)[-1,1]
cov_b1b1 <- cov_b1b1[1:ntiss]
cov_b1b1 <- cov_b1b1*Gt
return(list("covariance" = cov_b1b1, "weights" = weights))
}

get_cov_alpha1alpha2_multitissue_withcrude <- function(alpha_z1,alpha_z2,gtot,n,x,genenames,x_uncorrected){

mean_coreg <- mean(rowSums(x), na.rm = T) #modified mar 8 from just x
mean_chisq <- mean(alpha_z1^2, na.rm = T)
crude_h2_1 <- (mean_chisq - 1)/(n*mean_coreg) #per-SNP heritability = quantity h2/M, therefore no M in equation; mean tau.
mean_chisq <- mean(alpha_z2^2, na.rm = T)
crude_h2_2 <- (mean_chisq - 1)/(n*mean_coreg)

y <- (alpha_z1*alpha_z2)*(gtot/n) 
#1.rho_pheno*NS=0, rho_g_est = the below
rho_g_est <- sum(alpha_z1*alpha_z2)/n #the first regression uses weights with rho_g_est and rho_pheno = 0. 
weights <- get_weights_0118(x,genenames,n,crude_h2_1,crude_h2_2,rho_g_est,0,gtot,x_uncorrected)
if(is.null(nrow(x))){rs <- x}else{rs <- rowSums(x)}
mod <- summary(glm(y ~ rs, weights = weights)) #modified mar 8 from + covars
#get estimates for rho_pheno*Ns and rho_g
intercept <- coef(mod)[1,1] #rho_pheno*Ns/sqrt(N1N2)
rhop_NS <- intercept*n #cannot be greater in magnitude than n*1 and it is if we do multi variate regression. using rowsums, then we get -777, which means pheno corr = -0.78, but we know it is 0.03. therefore something is wrong. we get 3K > 1K. 
rho_g_est2 <- coef(mod)[2,1] #need a single estimate, not one per tissue, so maybe first regression should not be multi variate. or just take the estimate for the causal tissue? 

weights <- get_weights_0118(x,genenames,n,crude_h2_1,crude_h2_2,rho_g_est2,rhop_NS,gtot,x_uncorrected)
mod <- summary(glm(y ~ x, weights = weights)) ##modified mar 8 from + covars
cov_b1b1 <- coef(mod)[-1,1]
cov_b1b1 <- cov_b1b1[1:ntiss] 
return(list("covariance" = cov_b1b1, "weights" = weights))
}

get_weights_0827 <- function(x,genenames,alpha_z,ngwas,x_uncorrected){
t_genename <- table(genenames)
tiss_per_gene <- as.numeric(t_genename)[match(genenames, names(t_genename))]
if(is.null(nrow(x))){rs <- x}else{rs <- rowSums(x)}
if(is.null(nrow(x_uncorrected))){rs_uncorrected <- x_uncorrected}else{rs_uncorrected <- rowSums(x_uncorrected)}
mean_chisq <- mean(alpha_z^2, na.rm = T)
mean_coreg <- mean(rs, na.rm =T) #this is correct; was wrong in get_cov_alpha1alpha2_multitissue_withcrude
crude_h2_est <- (mean_chisq - 1)/(ngwas*mean_coreg)

heteroscedasticity <- (1 + ngwas*crude_h2_est*rs)^2
w1 <- 1/(heteroscedasticity) #variance/heteroskedasticity
w2 <- 1/rs_uncorrected #(1 + rs) #total co-reg across tissues; March 8 - this must be un-biascorrected sum, do we have something like "tissueassign" variable?  
w3 <- 1/tiss_per_gene #num tissues representing this gene. 
weights <- w1*w2*w3
weights[weights < 0] <- 0
return(weights)
}

get_weights_0118 <- function(x,genenames,ngwas,crude_h2_1,crude_h2_2,crude_rg,rhop_NS,gtot,x_uncorrected){
t_genename <- table(genenames)
tiss_per_gene <- as.numeric(t_genename)[match(genenames, names(t_genename))]
if(is.null(nrow(x))){rs <- x}else{rs <- rowSums(x)}
if(is.null(nrow(x_uncorrected))){rs_uncorrected <- x_uncorrected}else{rs_uncorrected <- rowSums(x_uncorrected)}
het_prod1 <- 1 + ngwas*crude_h2_1*rs #maybe shouldn't have /gtot here because our crude_h2_1 is a perSNP estimate; in cross trait supp, I think h2 is a genomewide estimate and they divide by M. 
het_prod2 <- 1 + ngwas*crude_h2_2*rs
het_prod3 <- ngwas*crude_rg*rs/gtot #leave gtot beause crude_rg is not a per snp quantity
het_prod4 <- rhop_NS/ngwas
heteroscedasticity <- het_prod1*het_prod2 + (het_prod3 + het_prod4)^2
w1 <- 1/(heteroscedasticity) #variance/heteroskedasticity
w2 <- 1/rs_uncorrected #total co-reg across tissues # march 8 update 
w3 <- 1/tiss_per_gene #num tissues representing this gene. 
weights <- w1*w2*w3
weights[weights < 0] <- 0
return(weights)
}

run_TCSC_Rg <- function(alpha_z1, alpha_z2, gencorralphamat,gencorrcount, gtot, n, samp, allw_h2genes, x, genenames, aa, x_uncorrected, ntiss, cish2genes){
gencorralphamat[gencorrcount,1] <- paste0("VariedParam",samp) #usually sample size
if(is.null(nrow(x))){x <- matrix(x,ncol = 1)
x_uncorrected <- matrix(x_uncorrected,ncol=1)}
gencorralphamat[gencorrcount,2:(ntiss + 2- 1)] <- get_cov_alpha1alpha2_multitissue_withcrude(alpha_z1,alpha_z2,gtot,n,x,genenames,x_uncorrected)$covariance
#gencorralphamat[gencorrcount,2:(ntiss + 2- 1)] <- get_cov_alpha1alpha2_multitissue_withcrude_SigCish2Genes(alpha_z1,alpha_z2,gtot,n,x,genenames,x_uncorrected,cish2genes)$covariance

#jackknife se 
jk <- matrix(0,nrow = chunks,ncol = ntiss)
jk_weights <- matrix(0,nrow = chunks,ncol = 1)
jkrefmat <- aa[allw_h2genes,]
for (chunk in 1:chunks){
remove_genes <- which(jkrefmat[,3] == chunk) #these will correspond to different tissues though. So different tissues will have different changes in ngenes

if(length(remove_genes) > 0){ #often no genes removed. 
alpha_z1_update <- alpha_z1[-remove_genes]
alpha_z2_update <- alpha_z2[-remove_genes]
#if(is.null(nrow(x))){x_update <- x[-remove_genes]
#x_update <- matrix(x_update,ncol=1)
#}else{
x_update <- x[-remove_genes,]#}
#if(is.null(nrow(x_uncorrected))){x_uncorrected_update <- x_uncorrected[-remove_genes]
#x_uncorrected_update <- matrix(x_uncorrected_update,ncol=1)
#}else{
x_uncorrected_update <- x_uncorrected[-remove_genes,]#}
genenames_upd <- genenames[-remove_genes]
}else{
alpha_z1_update <- alpha_z1
alpha_z2_update <- alpha_z2
x_update <- x
x_uncorrected_update <- x_uncorrected
genenames_upd <- genenames
}
#covars <- rep(gtot,length(alpha_z1_update))
if(is.null(nrow(x_update))){x_update <- matrix(x_update,ncol=1)}
if(is.null(nrow(x_uncorrected_update))){x_uncorrected_update <- matrix(x_uncorrected_update,ncol=1)}
cov_alpha1alpha2 <- get_cov_alpha1alpha2_multitissue_withcrude(alpha_z1_update,alpha_z2_update,gtot,n,x_update,genenames_upd,x_uncorrected_update)
#cov_alpha1alpha2 <- get_cov_alpha1alpha2_multitissue_withcrude_SigCish2Genes(alpha_z1_update,alpha_z2_update,gtot,n,x_update,genenames_upd,x_uncorrected_update,cish2genes)
jk[chunk,] <- cov_alpha1alpha2$covariance
if(length(remove_genes) > 0){jk_weights[chunk,1] <- sum(cov_alpha1alpha2$weights,na.rm = T)}else{jk_weights[chunk,1] <- 0}
} #chunks 

jackknife_se <- sapply(1:ntiss, function(x) sqrt(wtd.var(jk[,x], jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0))))
gencorralphamat[gencorrcount,(ntiss + 2):ncol(gencorralphamat)] <- jackknife_se
return(gencorralphamat)
}

run_TCSC_2steps <- function(alpha_z, varalphamat, varcount, gtot, ngenes, n, tr, ve, samp, x, genenames, x_uncorrected, ntiss, allw_h2genes, cish2genes){ #actually provide genenames and not allw_h2genes
varalphamat[varcount,1] <- paste0("Trait",tr)
varalphamat[varcount,2] <- paste0("VariedParam",ve)
varalphamat[varcount,3] <- paste0("nQTL",samplesizes[samp])
varalphamat[varcount,4:(ntiss + 4 - 1)] <- get_cov_alpha1alpha1_multitissue_withcrude(alpha_z,gtot,n,x,genenames,x_uncorrected) #don't rep, beacuse already repeated.  #was withcrude
#varalphamat[varcount,4:(ntiss + 4 - 1)] <- get_cov_alpha1alpha1_multitissue_withcrude_SigCish2Genes(alpha_z,gtot,n,x,genenames,x_uncorrected,cish2genes)

jk <- matrix(0,nrow = chunks,ncol = ntiss)
jk_weights <- matrix(0,nrow = chunks,ncol = 1)
jkrefmat <- aa[allw_h2genes,]
for (chunk in 1:chunks){
remove_genes <- which(jkrefmat[,3] == chunk) #these will correspond to different tissues though. So different tissues will have different changes in ngenes

if(length(remove_genes) > 0){ #often no genes removed. 
alpha_z_update <- alpha_z[-remove_genes]
if(ntiss == 1){x_update <- x[-remove_genes]}else{x_update <- x[-remove_genes,]}
genenames_upd <- genenames[-remove_genes]
if(ntiss == 1){x_uncorrected_update <- x_uncorrected[-remove_genes]}else{x_uncorrected_update <- x_uncorrected[-remove_genes,]}
weights <- get_weights_0827(x_update, genenames_upd, alpha_z_update, n, x_uncorrected_update) #of genes in regression, not in block. 
jk_weights[chunk,1] <- sum(weights, na.rm =T)
}else{
alpha_z_update <- alpha_z
x_update <- x
x_uncorrected_update <- x_uncorrected
genenames_upd <- genenames
jk_weights[chunk,1] <- 0 #if no genes in chunk, that chunk has weight = 0 in the sd
}

#covars_jk <- rep(gtot,length(alpha_z_update)) 
jk[chunk,] <- get_cov_alpha1alpha1_multitissue_withcrude(alpha_z_update,gtot, n, x_update, genenames_upd,x_uncorrected_update)} 
#jk[chunk,] <- get_cov_alpha1alpha1_multitissue_withcrude_SigCish2Genes(alpha_z_update,gtot, n, x_update, genenames_upd,x_uncorrected_update,cish2genes)}

varalphamat[varcount,(ntiss + 4):(ncol(varalphamat)-1)] <- sapply(1:ntiss, function(x) sqrt(wtd.var(jk[,x], jk_weights[,1]))*sqrt(length(which(jk_weights[,1] != 0))))
return(varalphamat)
}
