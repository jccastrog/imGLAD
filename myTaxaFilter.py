#!/usr/bin/env python
'''
@name: myTaxaFilter.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 16-Mar-2017
@version: 1.0.0
@license: GNU General Public License v3.0.
please type "./myTaxaFilter.py -h" for usage help
'''

import os, sys, subprocess
import argparse
import random
import numpy as np
try:
	from Bio import SeqIO
	from Bio.Seq import Seq
except:
	sys.stderr.write('ERROR! Module BioPython (Bio) is missing please install it before running fitModel.py\n')
	sys.exit()

#=====================1.2 Initialize variables======================
#1.2.1 Parser variables=============================================
parser = argparse.ArgumentParser(description="myTaxaFilter.py: Use MyTaxa to filter un-specific regions of the target genome [jccastrog@gatech.edu]")
group = parser.add_argument_group('Required arguments')
group.add_argument('-t', action='store', dest='target', required=True, help='The target genome to in FASTA format.')
group.add_argument('-mt', action='store', dest='my_Taxa', required=True, help = 'Path to the MyTaxa binaries.')
group.add_argument('-mdb', action='store', dest='db', required=False, help='Path to the MyTaxa protein database.')
group = parser.add_argument_group('Alignment arguments')
group.add_argument('-p', action='store', dest='threads', required=False, default='1', help='Number of threads. (default : %(default)s)')

args = parser.parse_args()
#1.2.2 Global variables==============================================
fragmentLength = 1000
fragsName = 'fragments.fa'
alignName = 'myTAxaAlign.tbl'
path = os.path.dirname(os.path.realpath(__file__))
#=======================1.3 Define functions=========================
def genome_fragments(genomeFile,size,outFile):
	with open(genomeFile) as fastaFile:
		fragsFile = open(outFile,'w')
		for fastaParse in SeqIO.parse(fastaFile,"fasta"):
			ID = fastaParse.id
			seq = fastaParse.seq
			seqLen = len(seq)
			numFrags = int(seqLen/size)
			for i in range(1,numFrags+2):
				fragSeq = seq[size*(i-1):fragmentLength*i]
				fragsFile.write(ID+str(i)+'\n')
				fragsFile.write(fraqSeq+'\n')
		fragsFile.close()

'''2.0 Align the genome segments to the protein databaseDownload the gneomes from NCBI'''
#=======================2.1 Create fragments from the genome file========================
genome_fragments(args.target,fragmentLength,fragsName)
#==============2.2 Align the fragments to the MyTaxa database with diamond===============
os.system("diamond blastx -q "+args.target+" -d "+args.db+" -p "+args.threads+" --out "+alignName+" -f 6")
#================2.3 Format the allignment table with BlastTab.catsbj.pl=================





