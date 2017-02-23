#!/usr/bin/env python

'''
@name: probEstimate.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 21-Apr-2016
@version: 1.0
@license: artistic license 2.0
'''

#======================================================1.0 Import modules and initialize variables=========================================================
#1.1 Import modules
import os, sys, argparse, subprocess, re  #interact with the system, use regular expressions
from os.path import basename
import numpy #make calculations base
from Bio import SeqIO #parse and manipulate fasta files

#1.2 Initialize variables
parser = argparse.ArgumentParser(description="protEstimate.py: Estimation of bacterial genomes in biological samples [jccastrog@gatech.edu] Aplha optimized for E. coli") #EDIT FOR THE ACTUAL NAME
parser.add_argument('-m', action='store', dest='Map', required=True, default='', help='One or more Tabular BLAST files of reads vs genes (or contigs).',nargs='+')
parser.add_argument('-s', action='store', dest='Seq', required=False, default='', help='Subject sequences (ref) in FastA format.')
args=parser.parse_args() #parse all the arguments into the array "args"
genomeSize=0 #the genome size

#1.3 Define functions
def predProb(depth,breadth):
#	theta=numpy.array([-42.492028,855.165048,403.992113])
	theta=numpy.array([-48.576425,923.055095,433.755612])
	z=numpy.dot(numpy.array([1,depth,breadth]),theta)
	g= 1/(1+numpy.exp(-z))
	return(g)

#============================================2.0 Load the blast results to calculate sequencing depth and breadth==========================================
#2.1 Sequencing depth and breadth
#2.1.1 Parse reference
with open(args.Seq) as fastaFile:
	for fastaParse in SeqIO.parse(fastaFile,"fasta"):
		ID=fastaParse.id
                seq=fastaParse.seq
		genomeSize=genomeSize+len(seq)
		
#2.2.1 Parse mapping
print 'File'+'\t'+'Depth'+'\t'+'Breadth'+'\t'+'p-Value'
for file in args.Map:
	genPos=dict() #position directory
	wholeDepth=0 #the depth of the overall genome
	with open(file) as map: #open Mapping file
		lines=map.readlines() #read the file line by line
		for line in lines: #for each line
			line=line.rstrip('\n') #chomp
			fields=line.split('\t') #split the fields of the blast file
			alnLength=int(fields[3]) #the alignment length
			perIden=float(fields[2]) #the percentage of identity
			if alnLength>=20 and perIden>=95:
				subSta=int(fields[8]) #the subject start of the alignment
				subEnd=int(fields[9]) #the subjec end of the alignment
				wholeDepth=wholeDepth+alnLength
				if subSta < subEnd:
					keys = range(subSta,subEnd+1) #create an arrray of the corrdinates being mapped
				else :
					keys = range(subEnd,subSta+1) #create an arrray of the corrdinates being mapped
				for key in keys: #for each of the positions ask if it is in the position directory or not
					if key in genPos:
						continue
					else:
						genPos[key]=1
#		print sum(genPos.values())
#2.2.2 Calculate sequencing depth and breadth
	seqDepth=wholeDepth/float(genomeSize)
	seqBreadth=sum(genPos.values())/float(genomeSize)
#=======================================================3.0 Calculate the probability of presence========================================================
	predict=predProb(seqDepth,seqBreadth)
	print basename(file)+'\t'+str(seqDepth)+'\t'+str(seqBreadth)+'\t'+str(1-predict)
