# runs test for simulated data
# MUST be run with data/simulated/ as the working directory
sh ../common/1_simulate_sequences.sh 0 0 91409891 91409891
sh ../common/2_build_trees.sh HIV 1
sh ../common/3_analysis.sh trees 1 2 91409891 '' "Simulated w/o latency"
