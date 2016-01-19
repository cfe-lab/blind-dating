# runs test for simulated data with latent tips 
# MUST be run with data/latent_sim/ as the working directory
sh ../common/1_simulate_sequences.sh 1 19041989 19041989
sh ../common/2_build_trees.sh HIV 0
sh ../common/3_analysis.sh trees 2 1
