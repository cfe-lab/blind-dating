# runs test for simulated data with latent tips 
# MUST be run with data/latent_sim/ as the working directory

for l in 1 0.01 0.0001; do
for u in 1 0.01 0.0001; do
	sh ../common/1_simulate_sequences.sh $u $l 91409891 91409891
	sh ../common/2_build_trees.sh HIV 1
	sh ../common/3_analysis.sh trees 2 2
	mv latency\.csv latency\.${l}\.${u}\.csv
	mv plot\.pdf plot\.${l}\.${u}\.pdf
	mv hist\.pdf hist\.${l}\.${u}\.pdf
	mv stats\.csv stats\.${l}\.${u}\.csv
	mv data\.csv data\.${l}\.${u}\.csv
done;
done;