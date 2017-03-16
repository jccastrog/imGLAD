#!/usr/bin/env python
'''
@name: fitModel.py
@author: Juan C. Castro <jccastrog at gatech dot edu>
@update: 10-Mar-2017
@version: 1.0.4
@license: GNU General Public License v3.0.
please type "./fitModel.py -h" for usage help
'''


'''1.0 Import modules, define functions, and initialize variables'''
#========================1.1 Import modules=========================
import os, sys, subprocess 
import argparse 
import gzip 
import random 
import scipy.optimize as opt 
import numpy as np 
from numpy import loadtxt, where 
try :
	from Bio import SeqIO
	from Bio.Seq import Seq
	from Bio.Blast.Applications import NcbiblastxCommandline
except:
	sys.stderr.write('ERROR! Module BioPython (Bio) is missing please install it before running fitModel.py\n')
	sys.exit()
try:
	import screed 
except:
	sys.stderr.write('ERROR! Module "screed" is missing please install it before running fitModel.py\n')
	sys.exit()
try:
	import statsmodels.discrete.discrete_model as sm
except:
	sys.stderr.write('WARNING! Consider installing python module "statsmodels" for faster calculation of model parameters \n')
#=====================1.2 Initialize variables======================
#1.2.1 Parser variables=============================================
parser = argparse.ArgumentParser(description="fitModel: Using simulated metagenomes built with ART estimate the training parameters for detection in a metagenomic dataset [jccastrog@gatech.edu]")
group = parser.add_argument_group('Required arguments') 
group.add_argument('-t', action='store', dest='target', required=True, help='The target genome to be detected in FASTA format.')
group.add_argument('-sp', action='store', dest='sp', required=True, help = 'The species of the genome written as binomial name in quotations. (e.g. "Escherichia coli")')
group = parser.add_argument_group('Simulated datasets arguments') 
group.add_argument('-l', action='store', dest='genomes', required=False, help='The genomes list to create the training dataset')
group.add_argument('-s', action='store', dest='train_size', required=False, default=200, help='Number of genomes included in the training dataset (Incompatible with -l). (default : %(default)s)')
group.add_argument('-e', action='store', dest='training_examples', required=False, default=100, help='Number of training examples (metagenomic datasets) used to train the model by default 200 (100 positive and 100 negative examples. (default : %(default)s)')
group.add_argument('-j', action='store', dest='platform', required=False, default='illumina', help='The sequencing platform used to generate the reads in the datasets. (default : %(default)s)')
group.add_argument('-r', action='store', dest='num_reads', required=False, default=1000000, help='Number of reads per training example (Metagenomic dataset simulated). (default : %(default)s)')
group.add_argument('-d', action='store', dest='read_length', required=False, default=150, help='Average read length for the simulated datasets. (default : %(default)s)')
group = parser.add_argument_group('Mapping arguments') 
group.add_argument('-p', action='store', dest='prog', required=False, default='blat', choices=["blat","blastn"], help='Program to be used to align the metagenomes and the reference sequence. (default: %(default)s)')
group.add_argument('-i', action='store', dest='perc_identity', required=False, default='95', help='Percentage of identity to recruit a read to the genome. (default: %(default)s)')
group.add_argument('-a', action= 'store', dest='aln_length', required=False, default='100', help='Alingment length to recruit a read to the genome. (default: %(default)s)')
args = parser.parse_args() 
#1.2.2 Global variables==============================================
downloadDict=dict()
trainName = 'trainingValues.csv'
refName = 'trainingGenomes.txt'
paramName = 'parameters.txt'
detName = 'detectionLimit.txt'
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
	wholeDepth = 0 
	genPos = dict() 
	with open(mappingFile) as mapping:
		lines = mapping.readlines()
		for line in lines:
			line.rstrip('\n') 
                        fields = line.split('\t')
                        alnLength = int(fields[3])
                        perIden = float(fields[2])
                        if alnLength>=reqLength and perIden>=reqIden:
				subSta = int(fields[8])
				subEnd = int(fields[9])
				wholeDepth = wholeDepth+alnLength
				keys = range(subSta,subEnd+1)
				for key in keys:
					if key in genPos:
						continue
					else:
						genPos[key] = 1
	seqBreadth = sum(genPos.values())/float(genomeSize)
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
if not os.path.exists("_tempdir"):
	os.makedirs("_tempdir") 
os.system("curl ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/bacteria/assembly_summary.txt -o _tempdir/genomes.txt --silent")
os.system("awk -F '\t' -v OFS='\t' '{if($12==\"Complete Genome\") print $1, $8, $20}' _tempdir/genomes.txt > _tempdir/assembly_summary_complete_genomes.txt")
os.remove('_tempdir/genomes.txt')
summFile='_tempdir/assembly_summary_complete_genomes.txt'
#=======2.2 Download the genomes========
sys.stderr.write('Downloading the training genomes from NCBI...\n')
genomesFile = open('_tempdir/genomes.fna', 'w')
#2.2.1 Download the list genomes========
if args.genomes is not None:
	if args.train_size is not None:
		sys.stderr.write('ERROR! Argument -s ignored when using -l\n')
	genomes = [line.rstrip('\n') for line in open(args.genomes)]
	with open(summFile) as summary:
		lines = summary.readlines()
		refFile = open(refName, 'w')
		for line in lines:
			line = line.rstrip('\n')
			fields = line.split('\t')
			if fields[0] in genomes:
				fileName = fields[2].split('/')[-1]
				fileName = fields[2]+'/'+fileName+'_genomic.fna.gz'
				outName = "_tempdir/"+"_".join(fields[1].split(' '))+".fna.gz"
				refFile.write(os.path.basename(outName)+'\n')
				subprocess.call(["curl", fileName, "-o", outName, "--silent"])
				fastaName = outName.rstrip('.gz')
				zipRef = zip.open(outName, 'rb')
				fastaFile = open(fastaName, 'wb')
				fastaFile.write( zipRef.read() )
				zipRef.close()
				fastaFile.close()
				os.remove(outName)
				with open(outName.rstrip('.gz')) as inFile:
					for line in inFile:
				                genomesFile.write(line)
				os.remove(fastaName)
			os.remove("_tempdir/assembly_summary_complete_genomes.txt")
		refFile.close()
#2.2.2 Download the default genomes=====
else:
	with open(summFile) as summary:
		lines = random.sample(summary.readlines(),int(int(args.train_size)*1.1))
		refName = 'trainingGenomes.txt'
		refFile = open(refName, 'w')
		spTarget = args.sp.split(' ')
		genoCount = 0
		if len(spTarget)==1:
	                for line in lines:
				line = line.rstrip('\n')
                	        fields = line.split('\t')
	                        ftpName = fields[2].split('/')[-1]
        	                ftpName = fields[2]+'/'+ftpName+'_genomic.fna.gz'
				spRef = fields[1].split(' ')
				if genoCount==int(args.train_size):
					break
				elif str(spRef[0])==str(spTarget[0]):
					continue
				else :
					outName = "_".join(fields[1].split(' '))+".fna.gz"
		                        outName = outName.replace("(","")
	        	                outName = outName.replace(")","")
	                	        outName = outName.replace(":","")
		                        outName = outName.replace("/","_")
		                        outName = outName.replace("'","")
		                        outName = "_tempdir/"+outName
		                        refFile.write(os.path.basename(outName)+'\n')
		                        subprocess.call(["curl", ftpName, "-o", outName, "--silent"])
		                        fastaName = outName.rstrip('.gz')
		                        zipRef = gzip.open(outName, 'rb')
		                        fastaFile = open(fastaName, 'wb')
		                        fastaFile.write( zipRef.read() )
		                        zipRef.close()
		                        fastaFile.close()
		                        os.remove(outName)
	        	                with open(fastaName) as inFile:
		                                for line in inFile:
		                                        genomesFile.write(line)
		                        os.remove(outName.rstrip('.gz'))
					genoCount+=1
		elif len(spTarget)==2:
			for line in lines:
                                line = line.rstrip('\n')
                                fields = line.split('\t')
                                ftpName = fields[2].split('/')[-1]
                                ftpName = fields[2]+'/'+ftpName+'_genomic.fna.gz'
                                spRef = fields[1].split(' ')
				if genoCount==int(args.train_size):
                                        break
                                elif spRef[0]==spTarget[0] and spRef[1]==spTarget[1]:
                                        continue
                                else :
                                        outName = "_".join(fields[1].split(' '))+".fna.gz"
                                        outName = outName.replace("(","")
                                        outName = outName.replace(")","")
                                        outName = outName.replace(":","")
                                        outName = outName.replace("/","_")
                                        outName = outName.replace("'","")
                                        outName = "_tempdir/"+outName
                                        refFile.write(os.path.basename(outName)+'\n')
                                        subprocess.call(["curl", ftpName, "-o", outName, "--silent"])
                                        fastaName = outName.rstrip('.gz')
                                        zipRef = gzip.open(outName, 'rb')
                                        fastaFile = open(fastaName, 'wb')
                                        fastaFile.write( zipRef.read() ) 
                                        zipRef.close()
                                        fastaFile.close()
                                        os.remove(outName)
                                        with open(fastaName) as inFile:
                                                for line in inFile:
                                                        genomesFile.write(line)
                                        os.remove(outName.rstrip('.gz'))
					genoCount+=1
		else:
			sys.stderr.write('ERROR! Invalid species name')
			sys.exit()
                os.remove("_tempdir/assembly_summary_complete_genomes.txt")
                refFile.close()
genomesFile.close()
sys.stderr.write('Download finished a list of the genomes used can be found in "'+refName+'"...\n')
### Compare the genomes to the ANI DB ###

'''3.0 Create training datasets from random reads'''
sys.stderr.write('Constructing the training datasets...\n')
#===3.1 Create random reads for the selected genomes with ART===
#3.1.1 Create reads from the background genomes=================
if str(args.platform) == 'illumina':
	for i in range(1,int(args.training_examples)+1):
		numSeqs = int(int(args.num_reads)/int(args.train_size))
		os.system("art_"+str(args.platform)+" -ss HS25 -i _tempdir/genomes.fna -o _tempdir/simulatedNeg-"+str(i)+" -c "+str(numSeqs)+" -l "+str(args.read_length)+" > /dev/null")
		os.remove("_tempdir/simulatedNeg-"+str(i)+".aln")
	os.remove("_tempdir/genomes.fna")
	#3.1.2 Generate reads from the target genome at varying coverage
	for i in range(1,int(args.training_examples)+1):
		os.system("art_illumina -ss HS25 -i "+str(args.target)+" -o _tempdir/simulatedTarget-"+str(i)+" -f "+str(random.uniform(0.005,0.05))+" -l "+str(args.read_length)+" > /dev/null")
		os.remove("_tempdir/simulatedTarget-"+str(i)+".aln")
	#3.1.3 Spike the negative datasets with varying amounts of the target genome
		destFile = "_tempdir/simulatedPos-"+str(i)+".fq"
		with open(destFile, 'w') as outFile:
			with open("_tempdir/simulatedNeg-"+str(i)+".fq") as inFile:
				for line in inFile:
					outFile.write(line)
			with open("_tempdir/simulatedTarget-"+str(i)+".fq") as inFile:
				for line in inFile:
					outFile.write(line)
		os.remove("_tempdir/simulatedTarget-"+str(i)+".fq")
#================3.2 Convert the reads to fasta================
for i in range(1,int(args.training_examples)+1):
	fastqFileName = "_tempdir/simulatedPos-"+str(i)+".fq"
	fastaFileName = "_tempdir/simulatedPos-"+str(i)+".fa"
	fastq_to_fasta(fastqFileName,fastaFileName)
	os.remove(fastqFileName)
	fastqFileName = "_tempdir/simulatedNeg-"+str(i)+".fq"
	fastaFileName = "_tempdir/simulatedNeg-"+str(i)+".fa"
	fastq_to_fasta(fastqFileName,fastaFileName)
	os.remove(fastqFileName)

'''4.0 Align the datasets to the target genome (Using BLAST of BLAT)'''
sys.stderr.write('Recruiting reads to the target genome...\n')
if not os.path.exists('_tempaln'):
	os.makedirs('_tempaln')
if args.prog=="blastn":
	sys.stderr.write('Formating BLAST database...\n')
	#============4.1 Create a database for the target==============
	if not os.path.exists('_tempdb'):
		os.makedirs('_tempdb')
	try:
		os.system("makeblastdb -in "+args.target+" -input_type fasta -dbtype nucl -out _tempdb/targetDB > /dev/null")
		#======4.2 Align the reads of the datasets to the target=======
		sys.stderr.write('Recruiting reads to the target...\n')
		for filename in os.listdir("_tempdir/"):
			if filename.endswith(".fa"):
				blastCMD = NcbiblastxCommandline(cmd='blastn', outfmt=6, query="_tempdir/"+filename, db='_tempdb/targetDB', evalue=0.1, out="_tempaln/"+os.path.splitext(filename)[0]+".tbl")
				blastCMD()
		os.system("rm -rf _tempdb _tempdir")
	except:
		sys.stderr.write('ERROR! Could not find "makeblastdb" or "blastn", make sure than the blast binaries are added to your $PATH!\n')
		os.system("rm -rf _tempdb _tempaln _tempdir")
		sys.exit()
elif args.prog=="blat":
	try:
		for filename in os.listdir("_tempdir/"):
			if filename.endswith(".fa"):
				os.system("blat "+args.target+" _tempdir/"+filename+" _tempaln/"+os.path.splitext(filename)[0]+".tbl -t=dna -out=blast8 > /dev/null")
#		os.system("rm -rf _tempdir")
	except:
		sys.stderr.write('ERROR! Could not find "blat" make sure the Blat binaries are added to your $PATH!\n')
		os.system("rm -rf _tempaln _tempdir")
		sys.exit()
else:
	sys.stderr.write('Invalid option -p '+str(args.prog)+' use either blat or blastn')
	os.system("rm -rf _tempaln _tempdir")
	sys.exit()

'''5.0 Calculate the sequencing depth and breadth of the training datasetets'''
sys.stderr.write('Calculating the training values...\n')
genomeSize = get_genome_size(args.target)
trainFile = open(trainName, 'w')
for filename in os.listdir("_tempaln/"):
	if filename.startswith("simulatedPos"):
		targetPresent = 1
	elif filename.startswith("simulatedNeg"):
		targetPresent = 0
	trainArr = breadth_depth(genomeSize , "_tempaln/"+filename, int(args.perc_identity), int(args.aln_length))
	trainArr.append(targetPresent)
	trainStr = str(trainArr[0])+","+str(trainArr[1])+","+str(trainArr[2])+"\n"
	trainFile.write(trainStr)
trainFile.close()
os.system("rm -rf _tempaln")

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
#=======================6.3 Calculate the parameters========================
paramFile = open(paramName, 'w')
#6.3.1 Calculate based on sequencing breadth only===========================
it = np.ones(shape=(m, n+1))
it[:, 1:n+1] = X
try:
	logit = sm.Logit(y, it)
	theta = logit.fit().params
except:
	iniTheta = np.ones(shape=(m, n+1))
	iniTheta[:, 1:n+1] = X
	theta=decorated_cost(it, y, n)
paramFile.write(str(theta[0])+","+str(theta[1])+","+str(theta[2])+"\n")
#6.3.2 Calculate based on sequencing breadth and depth======================
it = np.ones(shape=(m, n))
it[:, n-1] = X[:,0]
try:
	logit = sm.Logit(y, it)
	theta = logit.fit().params
except:
	iniTheta = np.ones(shape=(m, n+1))
	iniTheta[:, 1:n+1] = X
	theta=decorated_cost(it, y, n)
paramFile.write(str(theta[0])+","+str(theta[1])+"\n")
paramFile.close()
os.system("rm "+str(trainName))
sys.stderr.write('Saved paremeters can be found in parameters.txt\n')
#=====================6.4 Establish the detection limit=====================
detLimit = (2.944439 - theta[0])/theta[1]
detFile = open(detName, 'w')
detFile.write('Thanks for using imGLAD\n')
detFile.write('The detection limit for your '+args.sp+' genome is: '+str(detLimit)+' in terms of sequencing breadth of the genome. This means you need to sample at least '++str(detLimit*100)+' of the genome in order to detect it.') 
#===========================================================================
