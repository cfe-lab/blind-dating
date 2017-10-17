# runs test for simulated data
# SBATCH --cpus=4
# SBATCH --mem=2G

dataset=${SIM_9_TRUE}
rep=${SLURM_ARRAY_INDEX_ID}
src=~/blind-dating/src/
raxml=raxmlHPC-PTHREADS-AVX

echo "Building"
R --slave --silent -f ${src}/build.tree.R --args aligned/${pat_id}.fasta trees/${pat_id}.nwk ${raxml} 1 ${rep}
echo "Rooting"
R --slave --silent -f ${src}/root.tree.R --args trees/${pat_id}.nwk info/patient_1943.csv trees.rooted/${pat_id}.rtt.nwk 1
echo "Regressing"
R --slave --silent -f ${src}/regression.R --args --tree=trees.rooted/${pat_id}.rtt.nwk --info=info/patient_1943.csv --patid=${pat_id}
