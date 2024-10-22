conda activate onek1k

cd /home/benjamin/Documents/TCSC_onek1k/simulation_analysis/

ve_gene=0.1
ngwas=100

for sim in {1..100};do
  python Nov_sim_genetraits.py $sim $ve_gene $ngwas
done