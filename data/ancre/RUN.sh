# runs test for real mixed data
# MUST be run with data/ancre/ as the working directory
# currently does not work; ancre data needs to be reprocessed
sh ../common/2_build_trees.sh patient 1
read -p "Copy good trees to trees.good.ogr/ and trees.good.rtt/ and then press [Enter] to continue"
sh ../common/3_analysis.sh trees.good.ogr 0 0 91409891 .ogr
sh ../common/3_analysis.sh trees.good.rttall 0 2 91409891 .rttall "Plasma"
