#!/usr/bin/env python3

# This is a script to fetch all the HCV sequences in Genbank. It's
# adapted from Application 3 in Sayers, E. (2010). Entrez programming
# utilities help.  Bethesda (MD): National Center for Biotechnology
# Information (US).

from urllib.request import urlopen
from bs4 import BeautifulSoup
import os.path
import sys
import time

query = sys.argv[1]
email = sys.argv[2]
outfile = sys.argv[3]
retmax = 100

base_url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
params = "db=nucleotide&email=" + email + "&tool=python3"
esearch_params = "&usehistory=y"
efetch_params = "&rettype=gb&retmode=xml&retmax={}".format(retmax)
esearch_url = base_url + "esearch.fcgi?" + params + esearch_params
efetch_url = base_url + "efetch.fcgi?" + params + efetch_params

soup = BeautifulSoup(urlopen(esearch_url + "&term={}".format(query)))
count = int(soup.count.text)
web = soup.webenv.text
key = soup.querykey.text

print(str(count) + ' records found')

efetch_url += "&query_key={}&WebEnv={}".format(key, web)

outfile = open(outfile, "w")
retstart = 0
while retstart < count:
    print("Fetching records {} through {}".format(retstart, retstart+retmax))
    soup = BeautifulSoup(urlopen(efetch_url + "&retstart={}".format(retstart)))
    for record in soup.find_all("gbseq"):
        outfile.write(record.prettify())
    retstart += retmax
    outfile.flush()
    time.sleep(10)

