# runs test for real mixed data
# MUST be run with data/lanl/ as the working directory
sh ../common/2_build_trees.sh patient 1
read -p "Copy good trees to trees.good/ and then press [Enter] to continue"
sh ../common/3_analysis.sh trees.good 1 0
mv stats.csv stats.ogr.csv
mv plot.pdf plot.ogr.pdf
mv hist.pdf hist.ogr.pdf
sh ../common/3_analysis.sh trees.good 1 2
mv stats.csv stats.rttall.csv
mv plot.pdf plot.rttall.pdf
mv hist.pdf hist.rttall.pdf
