#!/usr/bin/env python3

# Read through all the XML files fetched from Genbank, filtering out
# only those with a year, and write them in FASTA format. Help from
# here:
# http://boscoh.com/programming/reading-xml-serially.html

import xml.etree.ElementTree as etree
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import re
import os.path
import sys

infile = sys.argv[1]
quallistin = sys.argv[2]
outfile = sys.argv[3]

out_fasta = open(outfile, "w")
year_ptn = re.compile("(19|20)[0-9]{2}")

def findtext(parent, tag):
    child = parent.find(tag)
    if child is not None:
        return child.text.strip()
    return None

quallist = quallistin.split(':')

records_processed = 0
for event, record in etree.iterparse(open(infile), events=("start", "end")):
	if event != "end":
		continue

	if record.tag != "gbseq":
		continue

	accession = findtext(record, "gbseq_locus")
	header = findtext(record, "gbseq_definition")
	comment = findtext(record, "gbseq_comment")
	all_text = header + comment if comment else header

	qualvalues = ["" for x in quallist]
	for table in record.iter("gbseq_feature-table"):
		for feature in table.iter("gbfeature"):
			qualifiers = feature.find("gbfeature_quals")
			if qualifiers is None: continue
			for qualifier in qualifiers.iter("gbqualifier"):
				qual_name = findtext(qualifier, "gbqualifier_name")
				qual_value = findtext(qualifier, "gbqualifier_value")
				if qual_name is None or qual_value is None: continue
				if qual_name in quallist:
					index = quallist.index(qual_name)
					qualvalues[index] = qual_value
					
				all_text += qual_value
	
	seq = record.find("gbseq_sequence").text.strip().upper()
	header = accession
	for x in qualvalues:
		header += "_" + x
	seqrecord = SeqRecord(Seq(seq), id=header, description="")
	SeqIO.write(seqrecord, out_fasta, "fasta")
	out_fasta.flush()
        
	record.clear()
	records_processed += 1
	if records_processed % 1000 == 0:
		print("Processed {} records".format(records_processed))
		
out_fasta.close()
