#!/usr/bin/env python
'''
@name: myTaxaFilter.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 31-Mar-2017
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
parser = argparse.ArgumentParser(description="myTaxaFilter.py: Use MyTaxa to filter un-specific regions of the target genome. Before running myTaxaFilter you should install MyTaxa (https://github.com/luo-chengwei/MyTaxa) as well as download the AllGenome database in diamond format (http://enve-omics.ce.gatech.edu/data/public_mytaxa/AllGenomes.faa.dmnd) [jccastrog@gatech.edu]")
group = parser.add_argument_group('Required arguments')
group.add_argument('-t', action='store', dest='target', required=True, help='The target genome to in FASTA format.')
group.add_argument('-mt', action='store', dest='my_taxa', required=True, help = 'Path to the MyTaxa binaries.')
group.add_argument('-mdb', action='store', dest='db', required=False, help='Path to the MyTaxa protein database.')
group = parser.add_argument_group('Alignment arguments')
group.add_argument('-p', action='store', dest='threads', required=False, default='1', help='Number of threads. (default : %(default)s)')

args = parser.parse_args()
#1.2.2 Global variables==============================================
fragmentLength = 1000
fragsName = 'fragments.fa'
alignName = 'myTaxaAlign.tbl'
myTaxaName = 'myTaxaIn.txt'
outName = 'myTaxaOut.txt'
regNameFile = 'regions.txt'
vals = []
path = os.path.dirname(os.path.realpath(__file__))
#=======================1.3 Define functions=========================
def genome_fragments(genomeFile,size,outFile):
	fragsFile = open(outFile,'w')
	with open(genomeFile) as fastaFile:
		for fastaParse in SeqIO.parse(fastaFile,"fasta"):
			ID = fastaParse.id
			seq = fastaParse.seq
			seqLen = len(seq)
			numFrags = int(seqLen/size)
			for i in range(1,numFrags+2):
				fragSeq = seq[size*(i-1):fragmentLength*i]
				fragsFile.write(">{}\n".format(i))
				fragsFile.write(str(fragSeq)+'\n')
	fragsFile.close()
#def 
def remove_regions(fragsFile,regionsFile,outFile):
	finalSeq = ""
	regsArr = []
	with open(regionsFile) as regs:
		lines = regs.readlines()
		for line in lines:
			line = line.rstrip('\n')
			regsArr.append('line')
	with open(fragsFile) as fastaFile:
		for fastaParse in SeqIO.parse(fastaFile):
				ID = fastaParse.id
				seq = fastaParse.seq
				if ID in regsArr:
					continue
				else:
					finalSeq = finalSeq+seq
	updatedGenome = open(outFile,'w')
	updatedGenome.write(">updatedGenome\n")
	updatedGenome.write(finalSeq+'\n')
	updatedGenome.close()
	os.remove(fragsFile)
	os.remove(regionsFile)

'''2.0 Align the genome segments to the protein databaseDownload the gneomes from NCBI'''
#=======================2.1 Create fragments from the genome file========================
genome_fragments(args.target,fragmentLength,fragsName)
#==============2.2 Align the fragments to the MyTaxa database with diamond===============
subprocess.call(["diamond", "blastx", "-q", args.target, "-d ", args.db, "-p", args.threads, "--out", alignName, "-f", "6"])
#================2.3 Format the allignment table with BlastTab.catsbj.pl=================
myTaxaFile = open(myTaxaName, 'w')
prepScript = args.my_taxa+"/BlastTab.metaxaPrep.pl"
subprocess.call([prepScript, "-q", "-f", "no genes.txt", alignName], stdout=myTaxaFile, stderr=FNULL)
os.remove(alignName)
'''3.0 Run myTaxa and select the regions to be removed'''
#=====================3.1 Run MyTaxa=====================
myTaxaBin = args.my_taxa+"/MyTaxa "
subprocess.call(myTaxaBin, myTaxaName, outName, "5") 
'''4.0 Remove the regions selected from the original genome'''
#============4.1 Check the values for each region=============
vals = []
with open(outName) as myTaxaFile:
	lines = myTaxaFile.readlines()
	for line in lines:
		line = line.rstrip('\n')
		if line.startswith('>'):
			continue
		else:
			fields = line.split('\t')
			if fields[2]=='NA':
				vals.append('NaN')
			else:
				vals.append(fields[2])
os.remove(myTaxaFile)
valsDist = np.array(vals)
removeTh = np.percentile(valsDist,5)
#=======4.2 Select the bottom 5 pecentile and remove it=======
regFile = open(regNameFile,'w')
lnNumber = 0
for i in range(0,len(vals)):
	lnNumber += 1
	if vals[i] <= removeTh:
		continue
	else :
		regFile.write("{}\n".format(i))
regFile.close()
remove_regions(fragsName, regNameFile, 'updatedGenome.fa')
#=============================================================