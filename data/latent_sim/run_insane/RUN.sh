# runs test for simulated data with latent tips 
# MUST be run with data/latent_sim/ as the working directory
sh ../../common/1_simulate_sequences.sh 0 1 91409891 91409891
sh ../../common/2_build_trees.sh HIV_run1 1
sh ../../common/3_analysis.sh trees 2 2 91409891 '' "Simulated w/ latency 1"