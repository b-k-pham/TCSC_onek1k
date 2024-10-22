conda activate onek1k

cd /home/benjamin/Documents/TCSC_onek1k/simulation_analysis/

ve_gene=0.1
ngwas=100

for sim in {1..100};do
  for ve_snp in 0 0.1 0.25;do
    Rscript Nov_sim_TWASsumstats.R ${sim}_${ve_gene}_${ve_snp}_${ngwas} #depends on complex trait and gene models
    Rscript TCSC_Nov_Sims_Group1_varynQTL_ggcoreg.R ${sim}_${ve_gene}_${ve_snp}_${ngwas}
  done
done