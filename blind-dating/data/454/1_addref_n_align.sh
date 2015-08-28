#!/bin/sh

for f in fasta/*.fasta; do 
	echo ${f#fasta/}
	
	# ## HXB2 aligment (ENV only)
	# cat reference/B.FR.83.HXB2_LAI_IIIB_BRU_K03455.fasta ${f} > ${f#fasta_unaligned/}.unaligned ;
	# python trim.py  ${f#fasta_unaligned/}.unaligned  > ${f#fasta_unaligned/}.trim 6045 8795 ;

	## B Ancestor 
	python trim.py  reference/B.ENV.ANCESTOR.fasta 700 1450 > ref.fasta  ;
	python trim.py ${f} 0 8795  > ${f#fasta/}.unaligned ;
	cat ref.fasta ${f#fasta/}.unaligned > ${f#fasta/}.trim ;

	muscle -in ${f#fasta/}.trim -out aligned/${f#fasta/} -maxiters 32
	
	# # MAFFT
	# #/usr/local/bin/mafft  --auto --inputorder ${f#fasta_unaligned/}.trim > aligned/${f#fasta_unaligned/} ;

	rm ${f#fasta_wo_ref/}.unaligned ${f#fasta_wo_ref/}.trim ;
done;
