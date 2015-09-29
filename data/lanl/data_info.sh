#!/bin/sh

# A small script to automate some of the data retreival for 
# table 1
for f in aligned/*; do
	g=${f#aligned/patient_}
	echo ${g%.fasta}
	echo " & "

	grep -Eo  PLASMA_[0-9]+$ ${f} | wc -l
	echo " & "
	grep -Eo  PBMC_[0-9]+$ ${f} | wc -l
	echo " & "
	grep -Eo  [0-9]+$ ${f} | wc -l
	echo " & "

	grep -Eo  PLASMA_[0-9]+$ ${f} | sort -u | wc -l
	echo " & "
	grep -Eo  PBMC_[0-9]+$ ${f} | sort -u | wc -l
	echo " & "
	grep -Eo [0-9]+$ ${f} | sort -u | wc -l
	echo " & Reference & Notes \\\\\\\\ "
done;