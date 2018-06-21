#!/usr/bin/env python

'''
@name: probEstimate.py
@author: Juan C. Castro <jccastrog at gatech dot edu> & William T. Harvey <wharvey31@gatech.edu>
@update: 09-Mar-2017
@version: 1.0.4
@license: GNU General Public License v3.0.
please type "./protEstimate.py -h" for usage help
'''

'''1.0 Import modules, define functions, and initialize variables'''
#========================1.1 Import modules=========================
try:
	import os, sys, argparse #interact with the system, use regular expressions
	from os.path import basename
	import numpy as np #make calculations base
	from Bio import SeqIO #parse and manipulate fasta files
	from scipy import stats
except:
	sys.stderr.write('ERROR! Cannot import required modules remember probEstimate.py requires os, sys, argparse, numpy, scipy, and Bio')

#=====================1.2 Initialize variables======================
parser = argparse.ArgumentParser(description="protEstimate.py: Estimation of bacterial genomes in biological samples [jccastrog@gatech.edu]")
group = parser.add_argument_group('Required arguments') #Required
group.add_argument('-t', action='store', dest='target', required=True, default='', help='Subject sequences (ref) in FastA format.')
group.add_argument('-m', action='store', dest='map', required=True, default='', help='One or more Tabular BLAST files of reads vs genes (or contigs).',nargs='+')
group.add_argument('-p', action='store', dest='param', required=True,help="Parameters file obtained from fitModel.py")
group = parser.add_argument_group('Optional arguments') #Required
group.add_argument('-l', action='store', dest='mode', required=False, default='single', choices=["single","general"],help='Number of parameters to be used to estimate the prescence probability: single mode will use only sequencing breadth whereas general mode will use sequencing depth and sequencing breadth. (default : %(default)s)')
group.add_argument('-i', action='store', dest='perc_identity', required=False, default='95', help='Percentage of identity to recruit a read to the genome. (default: %(default)s)')
group.add_argument('-a', action= 'store', dest='aln_length', required=False, default='135', help='Alignment length to recruit a read to the genome. (default: %(default)s)')
args = parser.parse_args()
genomeSize = 0
thetaB = []
thetaBD = []

#=======================1.3 Define functions========================
def predProb(X,theta):
	X = np.array(X)
	theta = np.array(theta)
	z = np.dot(X,theta)
	g = 1/(1 + np.exp(-z))
	return(g)

'''2.0 Load the blast results to calculate sequencing depth and breadth'''
#====================2.1 Sequencing depth and breadth=====================
#2.1.1 Parse reference====================================================
with open(args.target) as fastaFile:
	for fastaParse in SeqIO.parse(fastaFile,"fasta"):
		ID = fastaParse.id
                seq = fastaParse.seq
		genomeSize = genomeSize+len(seq)
	
#2.1.2 Parse parameters===================================================
with open(args.param) as paramFile:
	lines = paramFile.readlines()
	for line in lines:
		line = line.rstrip('\n')
		fields = line.split(',')
		if len(fields)==2:
			thetaB = []
			for field in fields:
				thetaB.append(float(field))
		elif len(fields)==3:
			for field in fields:
				thetaBD.append(float(field))
#2.1.3 Parse mapping=====================================================
print 'File'+'\t'+'Depth'+'\t'+'Breadth'+'\t'+'p-Value'
for file in args.map:
	genPos = dict()
	depthPos = np.zeros(genomeSize)
	wholeDepth = 0 
	with open(file) as map:
		lines = map.readlines()
		for line in lines:
			line = line.rstrip('\n')
			fields = line.split('\t')
			alnLength = int(fields[3])
			perIden = float(fields[2])
			if alnLength>=int(args.aln_length) and perIden>=float(args.perc_identity):
				subSta = int(fields[8])
				subEnd = int(fields[9])
				wholeDepth = wholeDepth + alnLength
				if subSta < subEnd:
					keys = range(subSta,subEnd+1)
					for key in keys:
						if key in genPos:
							depthPos[key]+=1
							continue
						else:
							genPos[key]=1
							depthPos[key]+=1
				else :
					keys = range(subEnd,subSta+1)
					for key in keys:
						if key in genPos:
							depthPos[key]+=1
							continue
						else:
							genPos[key]=1
							depthPos[key]+=1
	#2.3.2 Calculate sequencing depth and breadth====================
	seqDepth = wholeDepth/float(genomeSize)
	seqDepth = stats.trim_mean(depthPos,0.025)
	seqBreadth = sum(genPos.values())/float(genomeSize)
	'''3.0 Calculate the probability of presence'''
	if args.mode=='single':
		X = [1,seqBreadth]
		predict=predProb(X,thetaB)
	else:
		X = [1,seqBreadth,seqDepth]
		predict=predProb(X,thetaBD)
	print basename(file)+'\t'+str(seqDepth)+'\t'+str(seqBreadth)+'\t'+str(1-predict)
#======================================================================================
