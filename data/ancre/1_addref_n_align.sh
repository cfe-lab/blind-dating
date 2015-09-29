#!/bin/sh
r -f gen_sequences.R

for f in fasta_wo_ref/*.fasta; do 
	echo ${f#fasta_wo_ref/}
	
	# ## HXB2 aligment (ENV only)
	# cat reference/B.FR.83.HXB2_LAI_IIIB_BRU_K03455.fasta ${f} > ${f#fasta_unaligned/}.unaligned ;
	# python trim.py  ${f#fasta_unaligned/}.unaligned  > ${f#fasta_unaligned/}.trim 6045 8795 ;

	## B concensus 
	python trim.py  reference/B.ENV.ANCESTOR.fasta 700 1450 > ref.fasta  ;
	python trim.py ${f} 0 8795  > ${f#fasta_wo_ref/}.unaligned ;
	cat ref.fasta ${f#fasta_wo_ref/}.unaligned > ${f#fasta_wo_ref/}.trim ;

	muscle -in ${f#fasta_wo_ref/}.trim -out aligned/${f#fasta_wo_ref/} -maxiters 32
	
	# # MAFFT
	# #/usr/local/bin/mafft  --auto --inputorder ${f#fasta_unaligned/}.trim > aligned/${f#fasta_unaligned/} ;

	rm ${f#fasta_wo_ref/}.unaligned ${f#fasta_wo_ref/}.trim ;
done;


python trim.py  reference/B.ENV.ANCESTOR.fasta 370 1770 > ref.fasta  ;
python trim.py fasta_wo_ref/patient_427.fasta 0 8795  > patient_427.fasta.unaligned ;
cat reference/B.ENV.ANCESTOR.fasta patient_427.fasta.unaligned > patient_427.fasta.trim ;

muscle -in patient_427.fasta.trim -out aligned/patient_427.fasta -maxiters 32
	rm patient_427.fasta.trim patient_427.fasta.unaligned ;