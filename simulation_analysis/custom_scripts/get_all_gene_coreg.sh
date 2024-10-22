conda activate onek1k

cd /home/benjamin/Documents/TCSC_onek1k/simulation_analysis/

for sim in {1..100};do
  Rscript Nov_Get_Gene_CoRegMatrices.R $sim
done