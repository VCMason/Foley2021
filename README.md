# Foley2021
Scripts used for "A genomic timescale for placental mammal evolution" (Foley et al. 2021)


###############################################################
###############################################################
ConcatenateAlignmentBlocksToOneLine_NoMemoryManagement.py was used in publication. However, the _MemoryManagement version should generate the same result while minimizing RAM requirements while increasing run-time. _MemoryManagement write temporary files to the hard drive so it requires more hard drive space.

To run the ConcatenateAlignmentBlocks scripts:
The program needs three things.
1. The block alignment file for one chromosome
2. the species files
3. a specified reference species
To change these open the file and modify the variables in first few lines of code:

alignmentfile = 'chr22_version1.aln'  # test with chr22_version1.aln.n1000000.txt
speciesfile = 'Species.txt'
referencespecies = 'Homo_sapiens'
linelimit = 1000000

alignmentfile: fasta version of HAL alignment
speciesfile: lists all species to concatenate from alignmentfile and orders species in output file. One species name per line.
linelimit: specifies the size of the temporary files. Larger numbers makes them bigger, smaller.. smaller.

To run the program (requires python 3, and two modules "datetime" and "os" (these are likely already installed):
1. put the program file, the alignment file, and the species file in the same directory
2. open a command prompt and navigate to the directory and type:
python ConcatenateAlignmentBlocksToOneLine_MemoryManagement.py
or
python ConcatenateAlignmentBlocksToOneLine_NoMemoryManagement.py

I assumed in the program that humans (or whatever reference species you specify) is the first species in all alignment blocks since i use these human coordinates for all folllowing species (until the next human (reference) sequence).
I assumed that any species present in an alignment block has the same number of characters as the human sequence in that alignment block.
I assumed all sequences were oriented in relation to human (I did not reverse complement anything).


###############################################################
###############################################################
DeletionDetector.py
Script Editor: Andrew Harris

Example hypothesis file MammalHypotheses.txt

Hypothesis file format has three sections per hypothesis delimited by "|" . Indicidual species names are delimited by ", ". File follows a .fasta-like format

>hypothesis_1
All species in the test clade | Required species from the test clade | Outgroup species
All species in test clade may or may not have deletion | Required species from test clade must have deletion | Outgroup species must not have deletion

example:
>hypothesis_1
species1, species2, species3 | species2, species3 | species4


###############################################################
###############################################################

FilterAln.py
Performs sliding window analysis of alignments and masks highly divergent sequences (all sequences are compared to reference sequence) within windows with missing data characters.
