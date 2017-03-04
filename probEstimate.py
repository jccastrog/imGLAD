#!/usr/bin/env python

'''
@name: probEstimate.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 21-Apr-2016
@version: 1.0.4
@license: GNU General Public License v3.0.
please type "./protEstimate.py -h" for usage help
'''

'''1.0 Import modules, define functions, and initialize variables'''
#========================1.1 Import modules=========================
import os, sys, argparse, subprocess, re  #interact with the system, use regular expressions
from os.path import basename
import numpy as np #make calculations base
from Bio import SeqIO #parse and manipulate fasta files

#=====================1.2 Initialize variables======================
parser = argparse.ArgumentParser(description="protEstimate.py: Estimation of bacterial genomes in biological samples [jccastrog@gatech.edu]")
group = parser.add_argument_group('Required arguments') #Required
group.add_argument('-m', action='store', dest='map', required=True, default='', help='One or more Tabular BLAST files of reads vs genes (or contigs).',nargs='+')
parser.add_argument('-s', action='store', dest='seq', required=True, default='', help='Subject sequences (ref) in FastA format.')
roup = parser.add_argument_group('Optional arguments') #Required
group.add_argument('-p', action='store', dest='param', required=False, default=[-34.738273,550.229,1080.350],help="Parameters file obtained from fitModel.py")
group.add_argument('-l', action='store', dest='mode', required=False, default='single', choices=["single","general"],help='Number of parameters to be used to estimate the prescence probability: single mode will use only sequencing breadth whereas general mode will use sequencing depth and seuqencing breadth. (default : %(default)s)')
args=parser.parse_args() #parse all the arguments into the array "args"
genomeSize=0 #the genome size
thetaB = []
thetaBD = []

#=======================1.3 Define functions========================
def predProb(vars,theta):
	X = [1]
	for i in range(0,len(vars)):
		X.append(float(vars[i]))
	X = np.array(X)
	theta = np.array(theta)
	z = np.dot(X,theta)
	g = 1/(1 + np.exp(-z))
	return(g)

'''2.0 Load the blast results to calculate sequencing depth and breadth'''
#====================2.1 Sequencing depth and breadth=====================
#2.1.1 Parse reference====================================================
with open(args.seq) as fastaFile:
	for fastaParse in SeqIO.parse(fastaFile,"fasta"):
		ID=fastaParse.id
                seq=fastaParse.seq
		genomeSize=genomeSize+len(seq)
	
#2.2.1 Parse parameters===================================================
with open(args.param) as paramFile:
	lines = paramFile.readlines()
	for line in lines:
		line = line.rstrip('\n')
		fields = line.split('\t')
		if len(fields)==2:
			thetaB = []
			for field in fields:
				thetaB.append(float(field))
		elif len(fields)==3:
			for field in fields:
				thetaBD.append(float(field))
#2.3.1 Parse mapping=====================================================
print 'File'+'\t'+'Depth'+'\t'+'Breadth'+'\t'+'p-Value'
for file in args.map:
	genPos = dict() #position directory
	wholeDepth = 0 #the depth of the overall genome
	with open(file) as map: #open Mapping file
		lines = map.readlines() #read the file line by line
		for line in lines: #for each line
			line = line.rstrip('\n') #chomp
			fields = line.split('\t') #split the fields of the blast file
			alnLength = int(fields[3]) #the alignment length
			perIden = float(fields[2]) #the percentage of identity
			if alnLength>=20 and perIden>=95:
				subSta = int(fields[8]) #the subject start of the alignment
				subEnd = int(fields[9]) #the subjec end of the alignment
				wholeDepth = wholeDepth + alnLength
				if subSta < subEnd:
					keys = range(subSta,subEnd+1) #create an arrray of the corrdinates being mapped
				else :
					keys = range(subEnd,subSta+1) #create an arrray of the corrdinates being mapped
				for key in keys: #for each of the positions ask if it is in the position directory or not
					if key in genPos:
						continue
					else:
						genPos[key]=1
	#2.3.2 Calculate sequencing depth and breadth====================
	seqDepth=wholeDepth/float(genomeSize)
	seqBreadth=sum(genPos.values())/float(genomeSize)
	'''3.0 Calculate the probability of presence'''
	if args.mode=='single':
		X = [seqBreadth]
		predict=predProb(X,thetaB)
	else:
		X = [seqBreadth,seqDepth]
		predict=predProb(X,thetaBD)
	print basename(file)+'\t'+str(seqDepth)+'\t'+str(seqBreadth)+'\t'+str(1-predict)
#======================================================================================
