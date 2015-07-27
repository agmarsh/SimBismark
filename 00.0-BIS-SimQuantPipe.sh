#!/bin/bash
# ------------------------------------------------------------------
# Pipeline for Analysis of SIMULATED methyl-DNA seq reads.
# ------------------------------------------------------------------
# AGM-2014
# Modified 2015 Y.X
# ------------------------------------------------------------------

# When ReadGenFile is "1" with the goal of processing a standard genome sequence
#   for comparison to a GenPro seqread set, then the runtime variables need to be
#   matched to the comparative set:
#       genCopyNum . . . . epigenetic sequence complexity in sample pool
#       depthFragSampling  . . . . max potential read representation
#       REtoss . . . . . probability of a frag being sequenced

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#                R U N T I M E   V A R S
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
SeqID="10"                      # file prefix for simulated genome output
TAG="0001"                      # unique ID str for folders/files
gMBsize=252                     # genome size (MB)
genCopyNum=30                  # number of genome copies (each with different MET patterns)
depthFragSampling=1                     # sequencing cycles (depth of frag sampling to generate seq reads)
readLen=76                      # sequence read length for each tag
simRefGenFolder="001-RefGenome"   # use if working from project folder with the 00.0-SimQuantPipe.sh script
ReadGenFile=1                   # 0 = generate DNA seq; 1 = read DNA seq from an input file
fracBIS=99                      # percent efficiency of bisulfite conversion chemistry
fastQ=0                         # 0=fasta; 1 = fastQ format with phred scores
REtoss=100                      # % probability of a RE frag in the sample population being sequenced
shearToss=5                     # % probability of a sheared frag being sequenced
loadMET=1 						# use an existing MET Table and 01-GenomeSequence.fa file; found in ..001-RefGenome/$ID-$TAG/
loadRef="00-Hs37/Hs37-01-2MB.fa"   # location for a reference genome/chromosome for the alignment
geneNumber=1                    # number of gene seqs in input refgen file to be used
sizeSelect=150                  # minimum length retained during size selection for library synthesis
currDir=`pwd`					# current directory of this script

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#                F L O W   C O N T R O L
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 0 = NO, do not execute this step . . .
# 1 = YES, execute this step . . . . .
STEPALL=1       # Set all flags to 1 before entering pipe
STEP0=0          # Generate Simulated SeqRead data files
STEP1=0          # File/folder prep
STEP2=0          # BISMARK genome prep
STEP3=0          # BISMARK
STEP4=0          # BISMARK methylation extractor 
STEP5=0          # Prep data tables for R
STEP6=0          # Run R script for OBS vs EXP analysis

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "             ,,                                                        "
echo "\`7MM\"\"\"Yp,   db           .M\"\"\"bgd                              "                     
echo "  MM    Yb               ,MI    \"Y                                    "                    
echo "  MM    dP \`7MM  ,pP\"Ybd \`MMb.      .gP\"Ya    ,dW\"Yvd             "  
echo "  MM\"\"\"bg.   MM  8I   \`\"   \`YMMNq. ,M\`   Yb  ,W\`   MM          "  
echo "  MM    \`Y   MM  \`YMMMa. .     \`MM 8M\"\"\"\"\"\"  8M    MM         "
echo "  MM    ,9   MM  L.   I8 Mb     dM YM.    ,  YA.   MM                  "  
echo ".JMMmmmd9  .JMML.M9mmmP\` P\"Ybmmd\"  \`Mbmmd\`   \`Mbm:dMM            "  
echo "                                                   MM                  "  
echo "                                                 .JMML.                "
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "           Simulated Seq Read Processing and Analysis Pipeline"
echo "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -"
echo "GenPro/2015"

if [ $STEPALL = "1" ]; then
	STEP0=1          # Generate Simulated SeqRead data files . . . script 1.3
	STEP1=1          # 
	STEP2=1          # 
	STEP3=1          # 
	STEP4=1          # 
	STEP5=1          # 
	STEP6=1          # 
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 0. Run the Simulation Seq Read script . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP0 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "0. Running simulated seq read script to generate .qcfasta file . . . "
    if [ ! -d $simRefGenFolder/$SeqID-$TAG ]; then
        mkdir -p $simRefGenFolder/$SeqID-$TAG
    fi
    head -n 40 00.0-BIS-SimQuantPipe.sh > $simRefGenFolder/$SeqID-$TAG/00-RunConfig.txt
    python 10.1-GenerateBisulfiteSeqTagData.py $SeqID $TAG $gMBsize $genCopyNum $depthFragSampling $readLen $simRefGenFolder $ReadGenFile $fracBIS $fastQ $REtoss $shearToss $loadMET $loadRef $sizeSelect $geneNumber
	echo "- - - - - - - - - - - - - - - - - - - - - - - - - -" >> $simRefGenFolder/$SeqID-$TAG/00-RunConfig.txt
	echo "Number of simulated sequences generated"  >> $simRefGenFolder/$SeqID-$TAG/00-RunConfig.txt
	grep -c ">" $simRefGenFolder/$SeqID-$TAG/03-BIS-SeqReadData-76.fa >> $simRefGenFolder/$SeqID-$TAG/00-RunConfig.txt
	echo "- - - - - - - - - - - - - - - - - - - - - - - - - -" >> $simRefGenFolder/$SeqID-$TAG/00-RunConfig.txt

fi

#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 1. Setup folder and file locations . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP1 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "1. Setup files/folder Genome File . . . "
    if [ ! -d $simRefGenFolder/$SeqID-$TAG/01-Genome ]; then
        mkdir -p $simRefGenFolder/$SeqID-$TAG/01-Genome
    fi
    echo "cp $simRefGenFolder/01-GenomeSequence.fa $simRefGenFolder/$SeqID-$TAG/01-Genome/"
    cp $simRefGenFolder/$SeqID-$TAG/01-GenomeSequence.fa $simRefGenFolder/$SeqID-$TAG/01-Genome/      
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 2. Run Bismark Genome_Preparation Script . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP2 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "2. Prepare Genome File . . . "
        $currDir/000-Tools/bismark_v0.13.0/bismark_genome_preparation \
        --path_to_bowtie $currDir/000-Tools/bowtie2-2.2.4/ \
        --bowtie2 \
        $simRefGenFolder/$SeqID-$TAG/01-Genome/  
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 3. Run Bismark and map reads to reference target sequence . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP3 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "3. Running BISMARK . . . . "
    if [ ! -d 002-Bismark/$SeqID-$TAG/ ]; then
        mkdir -p 002-Bismark/$SeqID-$TAG/
    fi
        $currDir/000-Tools/bismark_v0.13.0/bismark \
        --fasta \
        --output_dir 002-Bismark/$SeqID-$TAG/ \
        --path_to_bowtie $currDir/000-Tools/bowtie2-2.2.4/ \
        --bowtie2 \
        $simRefGenFolder/$SeqID-$TAG/01-Genome/ \
        $simRefGenFolder/$SeqID-$TAG/03-BIS-SeqReadData-$readLen.fa
fi

#        --samtools_path $currDir/000-Tools/samtools-1.1/
#        --sam \
#           --fastq
#           --output_dir 002-Bismark/$TAG-$SeqID/ \


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 4. Run Bismark methylation extractor to generate reports . . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP4 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "4. Running BISMARK Methylation Extractor . . . . "
        $currDir/000-Tools/bismark_v0.13.0/bismark_methylation_extractor \
        --single-end \
        --output 002-Bismark/$SeqID-$TAG/ \
        --bedGraph \
        --zero_based \
        002-Bismark/$SeqID-$TAG/03-BIS-SeqReadData-$readLen.fa_bismark_bt2.sam
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 5. Merge the results into an R table for analysis . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP5 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "5. Merge all data sources into R table . . . . . "
    python 17.bis-ScoreSimMethyl.py $SeqID $TAG $loadMET
fi


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# 6. Generate R plots . . . .
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if [ $STEP6 = "1" ]; then
    echo "- - - - - - - - - - - - - - - - - - - - - - - - - -"
    echo "6. Generate R plots . . . . . "
    
    Rscript --vanilla 18.bis-SimMetPlots.R $SeqID $TAG "002-Bismark/"
    
fi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo ''
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo '* * * * *   J O B   R U N   D O N E   * * * * * * * * '
echo '- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
echo ''
echo ''

## EOF------------------------------------------

