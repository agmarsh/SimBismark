#!/usr/bin/python
"""
Analyze the Simulated Genome Methylation Profiles.

Compare the algorithm metric scores from the GenPro pipe to the defined levels
of 5'mC methylation in the simulated gDNA sample with differential methylation.

Input Files:
	1.  "000-RefGen/141020d-MetCountTable.txt"
		This file contains two columns. The first is the C<C>GG position in the
		ref genome sequence. The second is the differential methylation level as
		a percentage across all DNA copies in the sample.
	2. "006-Results/04-S2-RefGen-TABLE-CCGGpos-ALL.txt"
		Second column contains the C<C>GG position in the ref genome sequence.
		The score metrics are in the 10th and 11th columns (0 indexing).  
		
AGM:2014/2015
"""
import os
import re
import sys

# - - - - - - - - - - -   USER RUN TIME VARIABLES - - - - - - - - - - - - - -
id = '004'
tag = "R2MB"
loadMET = 0

# - - - - - - - - - - -   GLOBAL VARIABLES - - - - - - - - - - - - - -
#headers = "Pos\tExpMet\tObsMet\tMETbis\tUMTbis\tSeq\n"
headers = "Pos\tExpMet\tObsMet\tMETbis\tUMTbis\n"

if len(sys.argv) > 1:
	id  = sys.argv[1]
	tag   = sys.argv[2]
	loadMET   = int(sys.argv[3])
	
# REFtable     = "02-MetCountDataTable.txt"
REFtable     = "04-BIS-MetCountDataTable.txt"
bismarktable = "03-BIS-SeqReadData-76.fa_bismark_bt2.bismark.cov"
OUTfile      = "11-BisSimMethylScoreTable-CpG-"+tag+".txt"
OUTlost      = "12-BisSimMethylScoreTable-Lost-"+tag+".txt"
OUTother     = "13-BisSimMethylScoreTable-Other-"+tag+".txt"

results = "002-Bismark/%s-%s/"  % (id,tag)
refdir  = "001-RefGenome/%s-%s/" % (id,tag)

def GeneID(id):
	n = len(id)
	for i in range (1,(11-n)):
		id = "0"+id
	return id

class Ddict(dict):
    def __init__(self, default=None):
        self.default = default

    def __getitem__(self, key):
        if not self.has_key(key):
            self[key] = self.default()
        return dict.__getitem__(self, key)
RefMet   = Ddict(lambda: -1.0)
RefSeq   = Ddict(lambda: 'nada')
FoundMet  = Ddict(lambda: -1.0)
ScoreMet = Ddict(lambda: Ddict(lambda: 'na'))
ConCov   = Ddict(lambda: Ddict(lambda: 'na'))     # Initialize with empty list
CovPos   = []

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - -   M A I N  - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print "\n\nMatching %Methylation to Metric Scores . . . . . "

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Load the RefMet data . . . ..
print "    Loading Reference MET CpG Data . . . . . "
IN=open(refdir+REFtable, 'r')
FILE=IN.readlines()
IN.close()
for i in range(4): FILE.pop(0)
for line in FILE:
	x = line.rstrip().split('\t')
	RefMet[x[0]] = x[1]
	#if len(x) == 4:
	#	RefSeq[x[1]] = x[3]
	#else:
	#	RefSeq[x[1]] = 'na'
FILE=[]

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Match OBS and EXP %MET Data . . . . . . .
print "    Matching OBS and EXP %MET Data . . . .  "
IN=open(results + bismarktable, 'r')
OUT1=open(results + OUTfile, "w")
OUT1.write(headers)
OUT2=open(results + OUTother, "w")
OUT2.write(headers)
found = 0
while(1):
	line = IN.readline().rstrip()
	if not line: break
	# chr1|GenProSim|004-R2MB|2MB	958	959	57.0469798657718	85	64
	# line format: 0=header; 1=5mCpos; 2=1indexpos; 3=%met; 4=metcount; 5=unmetcount
	x = line.split("\t")
	pos = GeneID(str(int(x[1])-1))   # << 5mC position indexing is +1 indexed in bismark cov 	
	if RefMet[pos] > 0:
		met = float(x[3])
		# headers = "Pos\tExpMet\tObsMet\tMETscore\tUMTscore\tSeq\n"
		#OUT1.write("%s\t%s\t%0.1f\t%s\t%s\t%s\n" % (pos, RefMet[pos], met, x[4], x[5], RefSeq[pos]) )
		OUT1.write("%s\t%s\t%0.1f\t%s\t%s\n" % (pos, RefMet[pos], met, x[4], x[5]) )
		FoundMet[pos] = x[3]
		found += 1
	else:
	# Option for dumping all other 5mC calls . . . . . . 
		met = float(x[3])
		# headers = "Pos\tExpMet\tObsMet\tMETscore\tUMTscore\n"
		OUT2.write("%s\t%s\t%0.1f\t%s\t%s\tnada\n" % (pos, RefMet[pos], met, x[4], x[5]) )
		
IN.close()
OUT1.close()
OUT2.close()

OUT3=open(results + OUTlost,'w')
OUT3.write("POS\tpMET\tLOST\n")
lost = 0
for (pos, count) in sorted(RefMet.iteritems()):
	if count > 0:
		if FoundMet[pos] < 0:
			OUT3.write("%s\t%s\t0\n" % (pos, RefMet[pos]))
			lost += 1

OUT3.close()
print "    FOUND = ", found," ;  LOST = ", lost
print "\n\n\n * * * * *   D O N E   * * * * * * \n\n\n"

# EOF ------------------------------------------------------------------------
# AGM2010;2011;2012;2013;2014 . .  .   .     .      .       .        .         .          .

	
	



