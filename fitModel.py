#!/usr/bin/env python
'''
@name: fitModel.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 03-Mar-2017
@version: 1.0
@license: artistic license 2.0
please type "fitModel.py -h" for usage help
'''


'''1.0 Import modules, define functions, and initialize variables'''
#========================1.1 Import modules=========================
import os, sys, subprocess #Interact with the system
import argparse #Parse the arguments to the program
import gzip #Zip and unzip files
import random #Functions for randomization
import scipy.optimize as opt #Function optimization using scipy
import numpy as np #Numerical operations in python
from numpy import loadtxt, where #Load text data
from Bio import SeqIO #Work with Fasta files
from Bio.Seq import Seq #Work with sequence data
from Bio.Blast.Applications import NcbiblastxCommandline #Local BLAST for alignment
try:
	import screed #Short reads utils
except:
	sys.stderr.write('Module "screed" is missing please install it before running fitModel.py\n')
try:
	import statsmodels.discrete.discrete_model as sm
except:
	sys.stderr.write('Consider installing python module "statsmodels" for faster calculation of model parameters \n')
#=====================1.2 Initialize variables======================
#1.2.1 Parser variables=============================================
parser = argparse.ArgumentParser(description="fitModel: Using simulated metagenomes built with Grinder estimate the training parameters for detection in a metagenomic dataset [jccastrog@gatech.edu]")
group = parser.add_argument_group('Required arguments') #Required
group.add_argument('-t', action='store', dest='target', required=True, help='The target genome to be detected in FASTA format.')
group.add_argument('-sp', action='store', dest='sp', required=True, help = 'The species of the genome written as binomial name in quotations. (e.g. "Escherichia coli")')
group = parser.add_argument_group('Simulated datasets arguments') #Simulated datasets
group.add_argument('-l', action='store', dest='genomes', required=False, help='The genomes list to create the training dataset')
group.add_argument('-s', action='store', dest='train_size', required=False, default=200, help='Size of the training dataset in number of genomes (Incompatible with -l). (default : %(default)s)')
group.add_argument('-e', action='store', dest='training_examples', required=False, default=100, help='Number of training examples (metagenomic datasets) used to train the model by default 200 (100 positive and 100 negative examples. (default : %(default)s)')
group.add_argument('-j', action='store', dest='platform', required=False, default='illumina', help='The sequencing platform used to generate the reads in the datasets. (default : %(default)s)')
group.add_argument('-r', action='store', dest='num_reads', required=False, default=1000000, help='Number of reads per training example (Metagenomic dataset simulated). (default : %(default)s)')
group.add_argument('-d', action='store', dest='read_length', required=False, default=150, help='Average read length for the simulated datasets. (default : %(default)s)')
group = parser.add_argument_group('Mapping arguments') #Mapping
group.add_argument('-p', action='store', dest='prog', required=False, default='blastn', choices=["blastn","blat"], help='Program to be used to align the metagenomes and the reference sequence. (default: %(default)s)')
group.add_argument('-i', action='store', dest='perc_identity', required=False, default='95', help='Percentage of identity to recruit a read to the genome. (default: %(default)s)')
group.add_argument('-a', action= 'store', dest='aln_length', required=False, default='100', help='Alingment length to recruit a read to the genome. (default: %(default)s)')
args = parser.parse_args() #parse all the arguments into the array "args"
#1.2.2 Global variables==============================================
downloadDict=dict()
trainName = 'trainingValues.csv'
refName = 'trainingGenomes.txt'
paramName = 'parameters.txt'
#=======================1.3 Define functions=========================
def fastq_to_fasta(fastqFile,fastaFile):
	fasta = open(fastaFile,'w')
	for n, record in enumerate(screed.open(fastqFile)):
		del record['quality']
		fasta.write('>'+record['name']+'\n'+record['sequence']+'\n')
def get_genome_size(genomeFile):
	genomeSize = 0
	with open(genomeFile) as fastaFile:
		for fastaParse in SeqIO.parse(fastaFile,"fasta"):
			ID = fastaParse.id
			seq = fastaParse.seq
			genomeSize = genomeSize+len(seq)
	return(genomeSize)
def breadth_depth(genomeSize,mappingFile,reqIden,reqLength):
	wholeDepth=0 #The depth of the overall genome
	genPos=dict() #Position directory
	with open(mappingFile) as mapping: #Open Mapping file
		lines = mapping.readlines() #Read the file line by line
		for line in lines: #For each line
			line=line.rstrip('\n') #Chomp
                        fields=line.split('\t') #Split the fields of the blast file
                        alnLength=int(fields[3]) #The alignment length
                        perIden=float(fields[2]) #The percentage of identity
                        if alnLength>=reqLength and perIden>=reqIden:
				subSta=int(fields[8]) #The subject start of the alignment
				subEnd=int(fields[9]) #The subject end of the alignment
				wholeDepth=wholeDepth+alnLength
				keys= range(subSta,subEnd+1) #Create an arrray of the coordinates being mapped
				for key in keys: #For each of the positions ask if it is in the position directory or not
					if key in genPos:
						continue
					else:
						genPos[key]=1
	seqBreadth=sum(genPos.values())/float(genomeSize)
	seqDepth = wholeDepth/float(genomeSize)
	return([seqBreadth,seqDepth])	
def sigmoid(z):
    e = np.exp(1)
    den = 1.0 + e ** (-1.0 * z)
    d = 1.0 / den
    return d
def compute_cost(theta, X, y, n):
    theta.shape = (1, n)
    m = y.size
    h = sigmoid(X.dot(theta.T))
    J = (1.0 / m) * ((-y.T.dot(np.log(h))) - ((1.0 - y.T).dot(np.log(1.0 - h))))
    return -1*J.sum()
def compute_grad(theta, X, y, n):
    theta.shape = (1, n)
    grad = np.zeros(n)
    h = sigmoid(X.dot(theta.T))
    delta = h - y
    l = grad.size
    for i in range(l):
        sumdelta = delta.T.dot(X[:, i])
        grad[i] = (1.0 / m) * sumdelta * - 1
    theta.shape = (n,)
    return  grad
def decorated_cost(it, y, n):
    def f(theta):
        return compute_cost(theta, it, y, n)
    def fprime(theta):
        return compute_grad(theta, it, y, n)
    theta = np.zeros(n)
    return opt.fmin_bfgs(f, theta, fprime, disp=True, maxiter=400)

'''2.0 Download the gneomes from NCBI'''
#======2.1 Initialize the download======
if not os.path.exists("_tempdir"): #Create a temporary directory for the genomes and datasets
	os.makedirs("_tempdir") 
os.system("curl ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -o _tempdir/genomes.txt --silent") #Download the assembly_summary.txt from NCBI
os.system("awk -F '\t' -v OFS='\t' '{if($12==\"Complete Genome\") print $8, $20}' _tempdir/genomes.txt > _tempdir/assembly_summary_complete_genomes.txt") #Filter the complete genomes and get their URLs
os.remove('_tempdir/genomes.txt') #Remove the full genomes file
summFile='_tempdir/assembly_summary_complete_genomes.txt' #The summary file with the names and URL of the genomes
#=======2.2 Download the genomes========
sys.stderr.write('Downloading the training genomes from NCBI...\n')
genomesFile = open('_tempdir/genomes.fna', 'w')
#2.2.1 Download the list genomes========
if args.genomes is not None: #Is a list is provided
	if args.train_size is not None: #If the size is also specified
		sys.stderr.write('Error! Argument -s ignored when using -l\n') #Output anr error message
	genomes = [line.rstrip('\n') for line in open(args.genomes)] #Create a list with the genome names
	with open(summFile) as summary: #With the summary file
		lines = summary.readlines() #Read through the file
		refFile = open(refName, 'w') #The file with the reference files names
		for line in lines: #For each line in the summary file
			line = line.rstrip('\n') #Chomp
			fields = line.split('\t') #Split the fields by tabs
			if fields[0] in genomes: #If the genome name is in the list
				fileName = fields[1].split('/')[-1] #Split the URL by '/'
				fileName = fields[1]+'/'+fileName+'_genomic.fna.gz' #Add the file name and extension to be downloaded
				outName = "_tempdir/"+"_".join(fields[0].split(' '))+".fna.gz" #The output name
				refFile.write(outName+'\n') #Write the file name into the reference files file
				subprocess.call(["curl", fileName, "-o", outName, "--silent"]) #Download the genome
				fastaName = outName.rstrip('.gz') #The name of the fasta file
				zipRef = zip.open(outName, 'rb') #Unzip the file just downloaded
				fastaFile = open(fastaName, 'wb') #The destination fasta
				fastaFile.write( zipRef.read() )
				zipRef.close()
				fastaFile.close()
				os.remove(outName) #Remove the files
				with open(outName.rstrip('.gz')) as inFile: #Add the genome to genomes.fna
					for line in inFile:
				                genomesFile.write(line)
				os.remove(fastaName)
			os.remove("_tempdir/assembly_summary_complete_genomes.txt") #Remove the summary file
		refFile.close()
#2.2.2 Download the default genomes=====
else:
	with open(summFile) as summary: #With the summary file
		lines = random.sample(summary.readlines(),int(int(args.train_size)*1.1)) #Read thorough the file lines
		refName = 'trainingGenomes.txt'
		refFile = open(refName, 'w') #The file with the reference files names
		spTarget = args.sp.split(' ')
		genoCount = 0
		if len(spTarget)==1:
	                for line in lines: #For each line in the file
				line = line.rstrip('\n') #Chomp
                	        fields = line.split('\t') #Split the fields by tabs
	                        ftpName = fields[1].split('/')[-1] #Split the URL by '/'
        	                ftpName = fields[1]+'/'+ftpName+'_genomic.fna.gz' #Add the file name and extension to be downloaded
				spRef = fields[0].split(' ')
				if genoCount==int(args.train_size):
					break
				elif str(spRef[0])==str(spTarget[0]):
					continue
				else :
					outName = "_".join(fields[0].split(' '))+".fna.gz" #The output name
		                        outName = outName.replace("(","") #Remove special characters from the file name
	        	                outName = outName.replace(")","")
	                	        outName = outName.replace(":","")
		                        outName = outName.replace("/","_")
		                        outName = outName.replace("'","")
		                        outName = "_tempdir/"+outName
		                        refFile.write(outName+'\n') #Write the file name into the reference files file
		                        subprocess.call(["curl", ftpName, "-o", outName, "--silent"]) #Download the genome
		                        fastaName = outName.rstrip('.gz') # The name of the fasta file
		                        zipRef = gzip.open(outName, 'rb') #Unzip the file just downloaded
		                        fastaFile = open(fastaName, 'wb') #The destination fasta
		                        fastaFile.write( zipRef.read() )
		                        zipRef.close()
		                        fastaFile.close()
		                        os.remove(outName) #Remove the files
	        	                with open(fastaName) as inFile: #Add the genome to genomes.fna
		                                for line in inFile:
		                                        genomesFile.write(line)
		                        os.remove(outName.rstrip('.gz'))
					genoCount+=1
		elif len(spTarget)==2:
			for line in lines: #For each line in the file
                                line = line.rstrip('\n') #Chomp
                                fields = line.split('\t') #Split the fields by tabs
                                ftpName = fields[1].split('/')[-1] #Split the URL by '/'
                                ftpName = fields[1]+'/'+ftpName+'_genomic.fna.gz' #Add the file name and extension to be downloaded
                                spRef = fields[0].split(' ')
				if genoCount==int(args.train_size):
                                        break
                                elif spRef[0]==spTarget[0] and spRef[1]==spTarget[1]:
                                        continue
                                else :
                                        outName = "_".join(fields[0].split(' '))+".fna.gz" #The output name
                                        outName = outName.replace("(","") #Remove special characters from the file name
                                        outName = outName.replace(")","")
                                        outName = outName.replace(":","")
                                        outName = outName.replace("/","_")
                                        outName = outName.replace("'","")
                                        outName = "_tempdir/"+outName
                                        refFile.write(outName+'\n') #Write the file name into the reference files file
                                        subprocess.call(["curl", ftpName, "-o", outName, "--silent"]) #Download the genome
                                        fastaName = outName.rstrip('.gz') # The name of the fasta file
                                        zipRef = gzip.open(outName, 'rb') #Unzip the file just downloaded
                                        fastaFile = open(fastaName, 'wb') #The destination fasta
                                        fastaFile.write( zipRef.read() )
                                        zipRef.close()
                                        fastaFile.close()
                                        os.remove(outName) #Remove the files
                                        with open(fastaName) as inFile: #Add the genome to genomes.fna
                                                for line in inFile:
                                                        genomesFile.write(line)
                                        os.remove(outName.rstrip('.gz'))
					genoCount+=1
		else:
			sys.stderr.write('Error! Invalid species name')
			sys.exit()
                os.remove("_tempdir/assembly_summary_complete_genomes.txt") #Remove the summary file
                refFile.close()
genomesFile.close()
sys.stderr.write('Download finished a list of the genomes used can be found in "'+trainName+'"...\n')
### Compare the genomes to the ANI DB ###

'''3.0 Create training datasets from random reads'''
sys.stderr.write('Constructing the training dataset...\n')
#=3.1 Create random reads for the selected genomes with grinder=
#3.1.1 Create reads from the background genomes=================
if str(args.platform) == 'illumina':
	for i in range(1,int(args.training_examples)+1):
		numSeqs = int(int(args.num_reads)/int(args.train_size))
		os.system("art_"+str(args.platform)+" -ss HS25 -i _tempdir/genomes.fna -o _tempdir/simulatedPos-"+str(i)+" -c "+str(numSeqs)+" -l "+str(args.read_length)+" > /dev/null")#Positive examples
		os.system("art_"+str(args.platform)+" -ss HS25 -i _tempdir/genomes.fna -o _tempdir/simulatedNeg-"+str(i)+" -c "+str(numSeqs)+" -l "+str(args.read_length)+" > /dev/null")#Negative examples
		os.remove("_tempdir/simulatedPos-"+str(i)+".aln") #Remove the alignment file
		os.remove("_tempdir/simulatedNeg-"+str(i)+".aln") #Remove the alignment file
	os.remove("_tempdir/genomes.fna") #Remove the genomes file
	#3.1.2 Generate reads from the target genome at varying coverage
	for i in range(1,int(args.training_examples)+1): #For each training datasets
		os.system("art_illumina -ss HS25 -i "+str(args.target)+" -o _tempdir/simulatedTarget-"+str(i)+" -f "+str(random.uniform(0.01,0.1))+" -l "+str(args.read_length)+" > /dev/null") #Target datasets
		os.remove("_tempdir/simulatedTarget-"+str(i)+".aln") #Remove the alignment file
	#3.1.3 Spike the negative datasets with varying amounts of the target genome
		destFile = open("_tempdir/simulatedPos-"+str(i)+".fq","w")
		with open("_tempdir/simulatedTarget-"+str(i)+".fq") as inFile:
			for line in inFile:
				destFile.write(line)
		destFile.close()
#================3.2 Convert the reads to fasta================
for i in range(1,int(args.training_examples)+1):
	fastqFileName = "_tempdir/simulatedPos-"+str(i)+".fq"
	fastaFileName = "_tempdir/simulatedNeg-"+str(i)+".fa"
	fastq_to_fasta(fastqFileName,fastaFileName)

'''4.0 Align the datasets to the target genome (Using BLAST of BLAT)'''
sys.stderr.write('Formating BLAST database...\n')
#================4.1 Create a database for the target==================
if not os.path.exists('_tempdb'): #Create a temporary directory for the DB
        os.makedirs('_tempdb')
if args.prog=="blastn": #Check if the program to be used is blast
	os.system("makeblastdb -in "+args.target+" -input_type fasta -dbtype nucl -out _tempdb/targetDB > /dev/null")
#==========4.2 Align the reads of the datasets to the target===========
sys.stderr.write('Recruiting reads to the target...\n')
if not os.path.exists('_tempaln'): #Create a temporary directory for the alignments
        os.makedirs('_tempaln')
if args.prog=="blastn": #Check if the program to be used is blast
	for filename in os.listdir("_tempdir/"):
		if filename.endswith(".fa"):
			blastCMD = NcbiblastxCommandline(cmd='blastn', outfmt=6, query="_tempdir/"+filename, db='_tempdb/targetDB', evalue=0.1, out="_tempaln/"+os.path.splitext(filename)[0]+".tbl")
			blastCMD()
else :
	for filename in os.listdir("_tempdir/"):
		os.system("blat "+args.target+" _tempdir/"+filename+" _tempaln/"+os.path.splitext(filename)[0]+".tbl -t=dna -out=blast8 > /dev/null")

'''5.0 Calculate the sequencing depth and breadth of the training datasetets'''
sys.stderr.write('Calculating the training values...\n')
genomeSize = get_genome_size(args.target) #Calculate the size of the target genome
trainFile = open(trainName, 'w') #Create a new file for the values
for filename in os.listdir("_tempaln/"): #For each alignment file
	if filename.startswith("simulatedPos"): #If the set is positive
		targetPresent = 1
	elif filename.startswith("simulatedNeg"): #If the set is negative
		targetPresent = 0
	trainArr = breadth_depth(genomeSize , "_tempaln/"+filename, int(args.perc_identity), int(args.aln_length)) #Calcuate sequencin depth and breadth
	trainArr.append(targetPresent) #Append the value of presence
	trainStr=str(trainArr[0])+","+str(trainArr[1])+","+str(trainArr[2])+"\n" #A stirng with the comma separated values
	trainFile.write(trainStr) #Write the results as a csv
trainFile.close()

'''6.0 Determine the logistic model parameters based on the training data'''
sys.stderr.write('Training the logistic model...\n')
#===========================6.1 Load the dataset============================
data = loadtxt(trainName, delimiter=',')
X = data[:, 0:2]
y = data[:, 2]
#================6.2 Separate positive and begative examples================
pos = where(y == 1)
neg = where(y == 0)
m, n = X.shape
y.shape = (m, 1)
it = np.ones(shape=(m, n+1))
it[:, 1:n+1] = X
#=======================6.3 Calculated the parameters=======================
try:
	logit = sm.Logit(y, it)
	theta = logit.fit().params
except:
	iniTheta = np.ones(shape=(m, n+1))
	iniTheta[:, 1:n+1] = X
	costGrad = compute_cost(initial_theta, X, y)
	theta=decorated_cost(it, y, n)
paramFile = open(paramName, 'w')
paramFile.write(str(theta[0])+","+str(theta[1])+","+str(theta[2]))
os.system("rm "+trainName)
paramFile.close()
sys.stderr.write('Saved paremeters can be found in parameters.txt\n')
#===========================================================================
