####################################################################################################
### SAMPLE CONFIGURATION ###########################################################################
####################################################################################################

# experiment id (analysisID)
exid	test

### SAMPLE DATA
# Species
spec	imaginary
# gender
gend	asexual
# tissue
tiss	some
# library type
type	polya
# description
desc	test dataset for tophat/bowtie
# date (library?)
date	25.11.08/11:00
# read file
read	11.fq

### FILTERING
# required informative bases (above infoqual)
infobase	0
# quality threshold
infoqual	5
# tile dropped if quality goes below threshold befor this position (median)
infodrop	50

### BOWTIE
# alignments to report
bwtk	200
# discard read if more mappings
bwtm	50
# use mismatch count instead of quality (for spliced map only)
bwtv	-1
# use best stratum
bwtb	1
# seed length
bwtl	28
# reads are phred64 (else phred33)
bwts	0
# max mismatch quality score
bwte	70

### ANALYSIS
# use strand information
strands	0
# build new retros
retros	0
# strict filter for best stratum
bestmap	1
# exon duplication model (0:RZ, 1:RC, 2:RD 3:RB)
emod	1
# exon analysis tag
atag	1

### PAIREDEND
# is paired end
paired	0
# minimal pair distance
mindist	0
# maximal pair distance
maxdist	0

####################################################################################################
# REFERENCE DATABASE ###############################################################################
####################################################################################################
# reference genome filename
refs	ref

### SQL
# database
db	rex
# external annotation database (ensembl)
edb	homo_sapiens_core_51_36m
# sql connection parameters
sqlport	3306
sqluser	root
sqlhost	localhost
# protocol tables
proto	ANALYSIS
annotid	ANNOTATION
mapid	MAPPING

### MOCKREADS
# quality baseline for mockreads
quba	67666676777777777777777777777777666666666666666655555555554444333222111000//
# read length for reference and reads
readl	76

### ANNOTATION
# splice juntion file
splc	allJunctions.txt
# exon file
cexo	ConstitutiveExons.txt
# retro table (SQL)
rtab	none

### SPLICESITE HANDLING
# build stringent splice sites (force exonic within gene)
estr	1
# anchor length
ovrh	1
# minimum intron size
mintron	40

####################################################################################################
# OTHER ############################################################################################
####################################################################################################
# threads to use in bowtie
threads	4
# core exon tag (when building from ensembl)
ctag	1
# file prefixes
unmp	unmapped
unrp	repeat