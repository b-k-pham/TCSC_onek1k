conda activate onek1k

cd /home/benjamin/Documents/TCSC_onek1k/simulation_analysis/

ti=1
for sim in {1..100};do
  python Nov_sim_genemodels.py $sim $ti
done