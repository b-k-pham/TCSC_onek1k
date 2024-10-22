# Overarching Questions:
How many simulations are enough?
How many tissues are enough?
How large should the sample size be? Is 100 enough? It takes much longer when n >> 100.
Which version of the code is newer - April or November? For now, I've been using the november sims since those are released on GitHub.

# Processing Steps I Did:

```
0. Make sure 1KG_HM3_chr1 bed/bim/fam files are in directory. This is from Tiffany's cluster.
1. Rscript Nov_Sim_eQTL_EffectSizes.R
Changed paths from TCSC/Simulation_Analyis to ./
Create directory when does not exist
Only do 100 simulations.
2. python Nov_sim_get_imp_geno.py (simulate imputed genotype only once otherwise it will keep adding to the same file and wont over write) 
change from n = 1500 to n = 100.
Changed paths from TCSC/Simulation_Analyis to ./
3. git clone fusion_twas, copy paste it to simulation_analysis/ folder, rename to fusion_twas-master, 
In python Nov_sim_genemodels.py, change samplesizes from [100,200,300,500,1000,1500]
 to just [100] for brevity.
 Also, add function to make path when doesn't exist.
 
1 simulation for 1 tissue takes 2 minutes.
100 simulations for 1 tissue would take 200 minutes.
Run concurrently t1,t2,t3 - 100 simulations
Opened screens:
run_t1
run_t2
run_t3
for each screen, run 100 simulations per ti.
4. Run Nov_Get_Gene_CoRegMatrices.R from Tiffany's cluster.
Changed samplesizes to just 100 and 3 tissues.
Changed indicies for PartialCoReg[,oth + 1] <- genegenecoreg to PartialCoReg[,oth] <- genegenecoreg because index needs to go from 1:t instead of 0:t-1.
Run all 100sims for 3 tissues with get_all_gene_coreg.sh in custom_scripts/
5. Run Nov_sim_get_gwas_geno.py from Tiffany's cluster.
Since the sample size is 100, ensure that nGWAS  = 100. Just need to run python Nov_sim_get_gwas_geno.py 100 (like with the eQTL one, probably do only once!)
6. Run Nov_sim_genetraits.py from Tiffany's cluster. Here, sim ranges from 1-100, ve_gene = 0.1, and ngwas = 100. 
question to ask here is why b_qtls is fixed for first tissue?
I'm guessing this is the variance explained by gene expression.
Also added function to create filedirs if does not exist.
Use sim_gene_traits.sh script in custom_scripts/
7. Run Nov_sim_snptraits.py from Tiffany's cluster. sim ranges from 1-100, ve_snp = {0,0.1,0.25} and ngwas = 100
Also added function to create filedirs if does not exist.
Changed Y3_direct from Nov_Y3_direct_nCausalGenes100.txt.gz to Nov_Y1_alpha_nCausalGenes100.txt.gz. What is the difference here?
I'm guessing this is the variance explained by SNPs.
Use sim_snp_traits.sh script in custom_scripts/
8. Run Rscript Nov_MultiTissueAnalysis_MakeComplexTrait.R ${sim}_${ve_gene}_${ve_snp}_${ngwas} for each ve_snp
Added function to create filedirs if does not exist.

Use do_complex_traits.sh in custom_scripts/ to iterate all sims, ve_snps, with ve_gene = 0.1 and ngwas = 100


9. Use do_TCSC_Tiffany.sh  in custom_scripts/. sims 1:100, ve_gene = 0.1, ve_snp = {0,0.1,0.25}, ngwas = 100


do_TCSC_Tiffany.sh is running below as sims 1:100, ve_gene = 0.1, ve_snp = {0,0.1,0.25}, ngwas = 100
changed tissues from 0:2 instead of 0:9 and 100 sims.
Rscript Nov_sim_TWASsumstats.R ${sim}_${ve_gene}_${ve_snp}_${ngwas} #depends on complex trait and gene models
Rscript TCSC_Nov_Sims_Group1_varynQTL_ggcoreg.R ${sim}_${ve_gene}_${ve_snp}_${ngwas}


Variance results are in Group1res_Nov/. However, how do I parse it?

```
