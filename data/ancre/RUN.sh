# runs test for real mixed data
# MUST be run with data/ancre/ as the working directory
# currently does not work; ancre data needs to be reprocessed
sh ../common/2_build_trees.sh patient 0
read -p "Copy good trees to trees.good/ and then press [Enter] to continue"
sh ../common/3_analysis.sh trees.good 1 0
mv stats.csv stats.ogr.csv
mv plot.pdf plot.ogr.csv
mv hist.pdf hist.ogr.csv
sh ../common/3_analysis.sh trees.good 1 1
mv stats.csv stats.rtt.csv
mv plot.pdf plot.rtt.csv
mv hist.pdf hist.rtt.csv

