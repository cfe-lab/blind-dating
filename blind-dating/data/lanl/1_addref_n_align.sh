#!/bin/sh

for f in fasta_unaligned/*.fasta; do 
	echo ${f#fasta_unaligned/}
	
	# ## HXB2 aligment (ENV only)
	# cat reference/B.FR.83.HXB2_LAI_IIIB_BRU_K03455.fasta ${f} > ${f#fasta_unaligned/}.unaligned ;
	# python trim.py  ${f#fasta_unaligned/}.unaligned  > ${f#fasta_unaligned/}.trim 6045 8795 ;

	## B concensus 
	python trim.py  ${f} > ${f#fasta_unaligned/}.unaligned 6045 8795 ;
	cat reference/B.ENV.ANCESTOR.fasta ${f#fasta_unaligned/}.unaligned > ${f#fasta_unaligned/}.trim ;

	muscle -in ${f#fasta_unaligned/}.trim -out aligned/${f#fasta_unaligned/} -maxiters 32
	# MAFFT
	#/usr/local/bin/mafft  --auto --inputorder ${f#fasta_unaligned/}.trim > aligned/${f#fasta_unaligned/} ;

	rm ${f#fasta_unaligned/}.unaligned ${f#fasta_unaligned/}.trim ;
done;