#!/usr/bin/python
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
'''
Objective: Generate a model NGS seq read file from a model genome with defined
	cytosine methylation positions that has been treated to RRBS - Reduced Representation
	Bisulfite Sequencing.  The resulting seq files can then be used to
	assess algorithm accuracy.
	
	The key is generating differentially methylated CCGG sites among the frag
	copies such that the final DNA population sampled has %MET at CpG(i) ranging
	from 0% to 100%.
	
	The pool of seq reads generated is dependent upon a HpaII digest and a random
	shearing function to open ends for library preparation.
	
	Seq read output file is structured as a fasta file that has already been QC
	filtered so it has the ".qcfasta".
	
	Based off the RE SeqSim script from 2014. 
	
GenPro/AGM-2014
Modified 2015 Y.X
'''
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
import os
import re
import sys
import random     # use choice and randint
import string     # use maketrans and translate


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - -   U S E R    V A R I A B L E S  - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ID 	  		= "03"                      # file prefix for simulated genome output
TAG       	= "0720"                    # unique ID str for folders/files
gMBsize   	= 2                         # genome size (MB)
genCopyNum	= 100                       # number of genome copies (each with different MET patterns)
seqCycles 	= 1                     	# sequencing cycles (depth of frag sampling to generateseqreads)
readLen 	= 76                      	# sequence read length for each tag
runfolder	= "001-RefGenome"           # use if working from project folder with the 00.0-SimQuantPipe.sh script
ReadGenFile = 1				    		# 0 = generate DNA seq; 1 = read DNA seq from an input file
fracBIS		= 99                        # percent efficiency of bisulfite conversion chemistry
fastQ		= 0                         # 0=fasta; 1 = fastQ format with phred scores
REtoss 	    = 100                    	# % probability of a RE frag in the sample population being sequenced
shearToss 	= 0                      	# % probability of a sheared frag being sequenced
loadMET     = 0 						# use a ref file with an existing MET Table; found in ../$ID-$TAG
loadRef 	= "00-Hs37/Hs37-01-2MB.fa"     # location for a reference genome but no existing met table
geneNumber  = 400                       # number of gene seqs in refgen file to be used
sizeSelect  = 200                       # minimum length retained during size selection for library synthesis

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - -  G L O B A L   V A R I A B L E S  - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Set the NT pool frequency from which random characters will be drawn
#        HUMAN promoter domain composition:
#        A = 0.247; C = 0.251; G = 0.254; T = 0.248
freq = [247, 251, 254, 248]
pNT = ''
for nt in ["A", "C", "G", "T"]:
	fq   = freq.pop(0)
	pNT += ("%s " % nt) * fq

if len(sys.argv) > 1:
	# for ag in sys.argv: print ag
	ID          = sys.argv[1]
	TAG         = sys.argv[2]
	gMBsize     = float(sys.argv[3])            
	genCopyNum  = int(sys.argv[4]) 
	seqCycles   = int(sys.argv[5]) 
	readLen     = int(sys.argv[6])  
	runfolder   = sys.argv[7] 
	ReadGenFile = int(sys.argv[8])
	fracBIS     = int(sys.argv[9])
	fastQ       = int(sys.argv[10])
	REtoss    = int(sys.argv[11]) 
	shearToss   = int(sys.argv[12])
	loadMET     = int(sys.argv[13])
	loadRef     = sys.argv[14]
	sizeSelect  = int(sys.argv[15])
	geneNumber  = int(sys.argv[16])
reffolder = runfolder + "/"
runfolder += "/" +ID+"-"+TAG+"/"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Structure of genome if being generated randomly . . . . . . . . . 
genFragLen = int(gMBsize * 10**6)   # number NT positions in each large frag
CCGG = []
cmet = []
# . . . . . . . . . . . . . 
# Establish distribution of differential methylation values:
methylstates = [0 for i in range(5)] + [i*5 for i in range(20)] + [100 for i in range(5)]  # Possible set of fractional (%) CpG site methylation states
#methylstates = [0, 25, 50, 75, 100]
ATGC = [ "A", "A", "A", "T", "T", "T", "G", "G", "G", "G", "G", "G", "C" ]  # Limit random CCGG freq by controlling p(C) at second position
ATCG = [ "A", "A", "A", "T", "T", "T", "C", "C", "C", "C", "C", "C", "G" ]  # Limit random CCGG freq by controlling p(G) at fourth position
phred = "456789:;<=>?@ABCDEFGHIJ"
phredz = list(phred)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Defined Files - - - - - - - - - - - - - - - - - - - -
genfile = runfolder + "01-GenomeSequence.fa"
GenFile = reffolder + loadRef
if loadMET == 1:
	GenFile =  runfolder + "01-GenomeSequence.fa"
refMetTable =  runfolder + "02-MetCountDataTable.txt"
datafile    =  runfolder + "03-BIS-SeqReadData-%s.fa" % (str(readLen))	
bisMetTable =  runfolder + "04-BIS-MetCountDataTable.txt" 

runlog =  runfolder + "00-runlog.txt"

WIPE = open(datafile,"w")
WIPE.close()
WIPE = open(bisMetTable,"w")
WIPE.close()
LOG = open(runlog,'w')
LOG.close()
NT = pNT.split()

if loadMET == 0:
	OUT=open(refMetTable, 'w')
	OUT.write("# CG reference table for %MET (5mC) by position.\n# POS = 0 indexed \
and identifies the cytosine in a CG location on the CHROMOSOME.\n# The 5mC value is a \
percentage score between 0% and 100%\n")
	OUT.write("ID\tCGpos\tpMET-5mC\n")
	OUT.close()


RC    = string.maketrans("ACGTx", "TGCAy")
Mspxy = string.maketrans('xy','CG')
BiSO3 = string.maketrans('C','T')

# - - - - - - - - - - - - - - - - - - - - - -
class Ddict(dict):
    def __init__(self, default=None):
        self.default = default
    def __getitem__(self, key):
        if not self.has_key(key):
            self[key] = self.default()
        return dict.__getitem__(self, key)
refGEN   = Ddict( lambda: 'NA')
refMET   = Ddict( lambda: 0 )
refGAP   = Ddict( lambda: 0 )
CCGGall  = Ddict( lambda: 0 )
METtable = Ddict( lambda: 0 )
# - - - - - - - - - - - - - - - - - - - - - -
def DUMP(mssg):
	LOG = open(runlog,'a')
	LOG.write(mssg)
	LOG.close()
	print mssg
# - - - - - - - - - - - - - - - - - - - - - -
def GeneID(id):
	n = len(id)
	for i in range (1,(11-n)):
		id = "0"+id
	return id
# - - - - - - - - - - - - - - - - - - - - - -
def ReadFasta(file):
	FileFasta = open(file,'r')
	FileData = FileFasta.readlines()
	FileFasta.close()
	count = -1
	for line in FileData:
		if count < geneNumber:
			headerline = re.match('>', line)
			if headerline is not None:
				line = re.sub('>','',line)
				if geneNumber == 1:          # << just a single fasta entry sequence block
					d = line.split('|')
					chrnum = re.sub(r'^>','',d[0])
					pos = "1"                  # << use 1-indexing for position
					id = chrnum + "-" + pos 
				elif loadMET == 0:
					# Parse header line for HS37.fnt files . . . .
					d = line.split('|')
					DUMP("GenFile = %s\n" %(GenFile))
					pos = int(d[7]) - 3000
					headermatch = re.match(r'^Hs\d\d\-', d[0])
					if headermatch is None:
						print "\n\n* * * * * * * * * * * * * * * * * * *"
						print "ERROR: refgen input file headers can not be parsed."
						print "       check head code in ReadFasta function."
						print "* * * * * * * * * * * * * * * * * * *\n\n"
						sys.exit()
					h = d[0].split("-")
					chrnum = re.sub(r'^0*','',h[1])
					id = "chr" + chrnum + "-" + str(pos)
				else:
					# GenFile has already previously been processed:
					#    the correct 'id' string is last element in header line
					d = line.split('|')
					id = d[3].rstrip()
				seq = ''
				count += 1
			else:
				seq = seq + line.rstrip()
				refGEN[id] = seq
				# print id
	print "There were %d gene seqs found in refGenFile." % count
# - - - - - - - - - - - - - - - - - - - - - -

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#
#      '7MMM.     ,MMF'          db          `7MMF'    `7MN.   `7MF'    
#        MMMb    dPMM           ;MM:           MM        MMN.    M      
#        M YM   ,M MM          ,V^MM.          MM        M YMb   M      
#        M  Mb  M' MM         ,M  `MM          MM        M  `MN. M      
#        M  YM.P'  MM         AbmmmqMA         MM        M   `MM.M      
#        M  `YM'   MM        A'     VML        MM        M     YMM      
#      .JML. `'  .JMML.    .AMA.   .AMMA.    .JMML.    .JML.    YM
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

DUMP("\n\nRunning Genome Seq Tag Construction Set: ID = %s\n\n" % (ID))

DUMP("GenFile = %s\n" %(GenFile))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# . . . . . Generate the DNA sequence . . . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fracGC = 0
if (ReadGenFile == 0):
	fragSeq = 'AAAAAAAAAA'
	lastNT = 'A'
	while(len(fragSeq) < genFragLen-10): 			#genFraglen is the length of whole genome: size times 1000000
		ccgg = random.randint(0,100000)
		if ccgg > 99999:
			fragSeq += 'CCGG'
			fracGC += 4                  			#fracGC is the counter for number of G+C in genome seq
			lastNT = 'G'
		else:
			if lastNT == 'C':
				# Need to reduce random CCGG freq in model sequence; reduce p(C) in second pos
				atgc = random.choice(ATGC)
			elif lastNT == 'G':
				atgc = random.choice(ATCG)
			else:
				# Select from random NT distributions defined at start of script 
				atgc = random.choice(NT)      
			fragSeq += atgc
			if atgc == 'C' or atgc == 'G':
				fracGC += 1
			lastNT = atgc
			
	fragSeq += 'TTTTTTTTTT'
	DUMP("\n\nGenome Seq of %.2f MB.\n" % gMBsize)
	OUT = open(GenFile,'a')
	OUT.write("> Genome Model Sequence %s: %.2f MB\n%s\n" % (ID,gMBsize,fragSeq))
	OUT.close()
else:
	# Load starting sequence to work with . . . . . . 
	ReadFasta(GenFile)
	gMBsize = 0
	gcount = 0
	for seq in refGEN.itervalues():
		gMBsize += len(seq)
		gcount += seq.count('G')
		gcount += seq.count('C')
	fracGC = 100*float(gcount)/float(gMBsize)
	gMBsize = float(gMBsize)/10**6
	if loadMET == 0:
		# Copy the genome sequence file to 01-GenomeSequence.fa . . . . . . . .
		OUT = open(genfile,'w')
		for (pos,seq) in refGEN.iteritems():
			OUT.write(">chr1|Hs37|%s-%s|%s\n%s\n" % (ID,TAG,pos,seq))
		OUT.close()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# . . Produce differentially methylated copies in sample population . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
countCCGG = 0

if (loadMET == 0):
	# Generate the refmet table to be used for the population . . . . 
	for (id,fragSeq) in refGEN.iteritems():
		genFragLen = len(fragSeq)
		chromosome = re.match(r'(chr.{1,2})-(\d+)', id)
		chrid = chromosome.group(1)
		chrpos = chromosome.group(2)
		lastCCGG = 0
		index = 1
		while(index > -1):
			index = fragSeq.find('CCGG', index+2)
			if index > -1 and index < genFragLen:
				index += 1   # position marks the internal cytosine in CCGG
				posid = chrpos + "-" + GeneID(str(index))
				# The list "methylstates" sets probability value for the degree a site will be methylated
				# Dict key = "chr1-11363-4078" where chromosome id - gene chr position- CpG position within gene
				# Remember: gene positions are -3000 bp in refgen so they include 5' domains
				CCGGall[chrid+"-"+posid] = random.choice(methylstates)
				countCCGG += 1
				#refGAP[chrid+"-"+posid] = index - lastCCGG
				lastCCGG = index
				#print chrid+"-"+posid
	Mettable=open(refMetTable, 'a')
	for pos, pcnt in sorted(CCGGall.iteritems()):
		x = pos.split("-")
		Mettable.write("%s-%s\t%s\t%0.1f\n" % (x[0],x[1],x[2],pcnt))
	Mettable.close()

else:
	# Load Methylation data table . . . . . . . .
	print "Loading an existing RefMetTable for RefGenome sequences . . . . " 
	IN=open(refMetTable, "r")
	for i in range(4): drop=IN.readline()
	while(1):
		line=IN.readline().rstrip()
		if not line:
			break
		x = line.split("\t")
		# Dict key: x[0] = "chr1-11363" + "-" + "4078"
		# The refMetTable has x[0] = gen in chrom pos, x[1] = CpG in chrom pos
		#  << Need to convert chr pos back to gene pos to meld with CCGGall dict in script >>
		chrposmatch = re.search(r'\-(\d+)$', x[0])     # ex.starting x[0] = chr1-1 , x[1] = 89
		chrpos = int(chrposmatch.group(1))
		posCG = str(int(x[1]) - chrpos)            # when chrpos = 1, then this just becomes a zero-indexed position
		CCGGall[x[0] + "-" + GeneID(posCG)] = float(x[2])    # ex. key = chr1-0000000088
		countCCGG += 1
	if countCCGG == 0:
		print "\n\nERROR: the file %s is empty. A new table needs to be geneerated.\n\n" % refMetTable
		sys.exit()

DUMP("Genome Size = %d MB\nGenome G+C = %0.1f%%\nGenome CCGG = %d\n" % (gMBsize,fracGC,countCCGG))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# . . . Fragment Seq Read Population from multiple genome copies . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

MET   = 0   # uncut
UMT   = 0   # cut
library = 0      # counter for fake gene id numbers in headers

Ncount = len(refGEN)
for (id,fragSeq) in refGEN.iteritems():
	# 150513: Note that the refMET keying scheme only works for ONE chromosome at the moment.
	DUMP("%d. Fragment ID: %s . . . " % (Ncount,id) )
	Ncount -= 1
	for f in range(1,genCopyNum+1):
		# DUMP("               %d" % f)
		total = 0
		genFragLen = len(fragSeq)
		metSeq = list(fragSeq)
		chrPOS = 'x'
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# - - - - - - - - - Methylate - - - - - - - - - - - - - - -
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# Use a random selection from 0-100 to translate that probability to a real distribution
		for cgPosition,fracMET in sorted(CCGGall.iteritems()):
			cgPos = cgPosition.split("-")                    # << "chr1-1-0000000088
			geneFragID = re.search(r'%s$' % cgPos[1], id)    # 0=chrid, 1=genepos on chromosome, 2=CpG pos on chromosome
			if geneFragID is not None:
				chrPOS  = int(cgPos[1])
				ccggPOS = int(cgPos[2])
				Cposition = chrPOS + ccggPOS                 # now back to a 1-indexed position value
				metTHRESH = random.randint(0,99)
				#print "   **", fracMET, metTHRESH
				if (fracMET > metTHRESH):
					#testSeq = "".join(metSeq)
					#print ">>> %s <<<" % testSeq[ccggPOS-14:ccggPOS+10]
					#OLD: metSeq = metSeq[:ccggPOS] + 'x' + metSeq[ccggPOS+1:]
					# - - - - - - - - - - - - - - - - - - - - - - - - - - 
					metSeq[ccggPOS] = 'x'
					MET += 1
					refMET[GeneID(str(Cposition))] += 1
					# - - - - - - - - - - - - - - - - - - - - - - - - - - 
					#testSeq = "".join(metSeq)
					#print ">>> %s <<<" % testSeq[ccggPOS-14:ccggPOS+10]
					##OLD: print ">>> %s <<<" % metSeq[ccggPOS-14:ccggPOS+10]
					#print ">>>",chrPOS,ccggPOS,Cposition, fracMET, "\n\n"
					#sys.exit()
				else:
					refMET[GeneID(str(Cposition))] += 0    # !! need to establish position key value in dictionary
					UMT += 1
		metSeq = "".join(metSeq)
		
		# Reverse Complement of the sequence . . . . . . . . . .
		# Methylation site needs to be reconstructed as CxGG in RC strand
		#   "CxGG": reverse = "GGxC" > complement ="CCyG"
		#    Replace "CCyG"  with "CxGG" to move 5mC marker back to interior C position
		#    For MspI loop, CxCC is also a recognition site for a cut.
		rcseq = metSeq[::-1]
		rcseq = rcseq.translate(RC)
		rcseq = re.sub(r'CCyG', 'CxGG', rcseq)

		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# - - - - - - - - - Msp I  D i g e s t - - - - - - - - - -
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		Msp1 = []
		recogseq = [ 'CCGG', 'CxGG']    # not methyl-sensitive
		frag = 0
		# Process (+) strand . . . . . . . . . 
		position = [m.start(0) for m in re.finditer(r'(CCGG)|(CxGG)', metSeq)]
		if len(position) > 0:
			start = 0
			for i in position:
				fragment = metSeq[start:i+1]  + 'CG' #< end-repair
				start = i + 1 
				if len(fragment) > sizeSelect and len(fragment) < 300:
					Msp1.append(fragment)
			if len(metSeq[position[-1]+1:]) > sizeSelect and len(metSeq[position[-1]+1:]) < 300:
				 Msp1.append(metSeq[position[-1]:])
		# Repeat for RC (-) strand . . . . . . . . . 
		position = [m.start(0) for m in re.finditer(r'(CCGG)|(CxGG)', rcseq)]
		if len(position) > 0:
			start = 0
			for i in position:
				fragment = rcseq[start:i+1] + 'CG' #< end-repair
				start = i + 1
				if len(fragment) > sizeSelect and len(fragment) < 300:
					Msp1.append(fragment)
			if len(rcseq[position[-1]:]) > sizeSelect and len(rcseq[position[-1]+1:]) < 300:
				Msp1.append(rcseq[position[-1]:])
		
		#DUMP( "There are %d frags in Msp1 size-selected digest.\n" % len(Msp1) )
		
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# - - - Sampling Loop to provide coverage depth for Msp cuts 
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		for i in xrange(seqCycles):
			# Shear Msp frags - - - - - - -
			Shear = []
			for msp in Msp1:
				Shear.append(msp[0:readLen])
				# 150708: already size-selected so shearing isn't functional here
				#size = random.randint(2*readLen, 4*readLen)
				#start = size
				#for stop in xrange(size,len(msp)-size-1,size):
				#	Shear.append(msp[start:stop])
				#	start = stop
			# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
			# Build SeqTags - - - - - - - - -
			OUT = open(datafile,'a')
			for seqfrag in Shear:
				odds = random.randint(0,101)
				if (odds <= REtoss):
					seqfrag = seqfrag.translate(BiSO3)

					#index = seqrip.find('C',0)
					#while(index > -1):
					#	# Bisulfite conversion of 5'-OH-Cytosine to Uracil --> Thymidine
					#	# base x base conversion using reaction efficiency specified at runtime
					#	bisrex = random.randint(0,100)
					#	if (fracBIS > bisrex):
					#		seqrip = seqrip[:index] + 'T' + seqrip[index+1:]   # the end-index in python is non-inclusive.
					#	index = seqrip.find('C',index+1)
					#seqtrans = seqrip.translate(Mspxy)       # << replace x,y : C,C			

					seqfrag = seqfrag.translate(Mspxy)
					library += 1
					nid = GeneID("%d" % library)
					fasta = ''
					if fastQ == 1:
						fasta = "@GenProSimSeq-RRBS:%s-%s-%s\n%s\n+\n" % (TAG,ID,nid,seqfrag)
						for i in xrange(0,readLen):
							fasta += "%s" % (random.choice(phredz))
						fasta += "\n"
					else:
						fasta = ">GenProSimSeq-RRBS:%s-%s-%s\n%s\n" % (TAG,ID,nid,seqfrag)
					OUT.write(fasta)
					# print fasta
					# fasta = "> faketag-%s\n%s\n" % (nid,seqtrans[0:readLen])
					# OUT.write(fasta)
			OUT.close()
			
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		# - - - - Sampling loop to provide coverage for random shear sites
		# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		if shearToss > 0:
			interval = 2
			shift = 10
			for i in xrange(seqCycles):
				Shear = []
				# POSTIVE STRAND . . . . . . . 
				# add the ends and tail . . . . . . . .
				lead = re.sub('C', 'T', metSeq[0:readLen+1])
				tail = re.sub('C', 'T', metSeq[-readLen:])
				Shear.append(lead.translate(Mspxy))
				Shear.append(tail.translate(Mspxy))
				# Start shearing . . . . ..  ..  .. . . . 
				L = random.randint(shift,interval*readLen)
				while (L < genFragLen - readLen - 2):
					stop  = L + readLen + 1

					#index = shear.find('C',0)
					#while(index > -1):
					#	bisrex = random.randint(0,100)
					#	if (fracBIS > bisrex):
					#		shear = shear[:index] + 'T' + shear[index+1:]
					#	index = shear.find('C',index+1)
					
					shear = metSeq[L:stop].translate(BiSO3)
					Shear.append(shear.translate(Mspxy))
					L = L + random.randint(shift,interval*readLen)
			
				# NEGATIVE STRAND . . . . . . .
				# add the ends and tail . . . . . . . . 
				Shear.append(rcseq[0:readLen+1].translate(Mspxy))
				Shear.append(rcseq[-readLen:].translate(Mspxy))
				# Start shearing . . . . ..  ..  .. . . . 
				L = random.randint(shift,interval*readLen)
				while (L < genFragLen - readLen - 2):
					stop  = L + readLen + 1 
					shear = rcseq[L:stop].translate(BiSO3)
					Shear.append(shear.translate(Mspxy))
					L = L + random.randint(shift,interval*readLen)
				# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
				# Build SeqTags - - - - - - - - - - - - - - - - - -
				OUT = open(datafile,'a')
				for s in Shear:
					odds = random.randint(0,101)
					if (odds <= shearToss):
						library += 1
						nid = GeneID("%d" % library)
						fasta = ''
						if fastQ == 1:
							fasta = "@GenProSimSeq-RRBS:%s-%s-%s\n%s\n+\n" % (TAG,ID,nid,s[0:readLen])
							for i in xrange(0,readLen):
								fasta += "%s" % (random.choice(phredz))
							fasta += "\n"
						else:
							fasta = ">GenProSimSeq-RRBS:%s-%s-%s\n%s\n" % (TAG,ID,nid,s[0:readLen])
						OUT.write(fasta)
						# print fasta
				OUT.close()
				

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# OUTPUT Methylation Reference Table . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
OUTref=open(bisMetTable, 'w')
for pos, count in sorted(refMET.iteritems()):
	cov = 100*(float(count)/float(genCopyNum))
	OUTref.write("%s\t%0.1f\n" % (pos,cov))
OUTref.close()


DUMP("\n\n------------------------------------\nThere are %d seq tags in the library." % library)
TOT = MET + UMT
DUMP("\nThere are %d CCGG TOTAL sites.\n" % TOT)
DUMP("There are %d CCGG unmethylated sites.\n" % UMT)
DUMP("There are %d CCGG methylated sites.\n\n" % MET)
DUMP("\n\n * * *   D O N E   * * * \n\n")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  - - - -   E O F  - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


		
