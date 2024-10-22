conda activate onek1k

cd /home/benjamin/Documents/TCSC_onek1k/simulation_analysis/

ngwas=100

for sim in {1..100};do
  for ve_snp in 0 0.1 0.25;do
    python Nov_sim_snptraits.py $sim $ve_snp $ngwas
  done
done