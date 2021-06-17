''' Author: Victor Mason '''
''' Date: April 14, 2015 '''
''' Note: MMU is Macaca mulatta '''
''' Description: Records replaces an individual's sequence within a window that has pairwise deletion distance > cutoff value and replaces nucleotide bases ONLY (not gaps) with missingchar. '''
''' Description: Pairwise deletion is calculated by comparing reference individuals sequence with that of each other individual. Pairwise comparisons that have '-' or missingchar are not included in the calculation.''' 
''' Output:  '''

import re
import os
import glob
import math

fileextension = 'trim.fa'
refname = 'Homo_sapiens' # name of reference sequence for all alignments. We used human.
cutoff = 0.80 # pairwise deletion distance frequency cutoff value to mask all nucleotides in window when pairwise distance is greater than pairwise distance cutoff
Zscore = 6 # Zscore # significance threshold # The Maximum number of standard deviations away from mean for a window's sequence to be kept in analysis. # z = (x - mean) / sigma # 2.575
AAaligncutofflen = 20 # required length of AA alignment to be filtered and passed to output directory. If less than AAaligncutofflen alignment is excluded.
PCcutoff = 0 # minimum accepted percent base coverage. If one idividuals sequence is < PCcutoff exclude entire alignment.
windsize = 50
step = 1
phybindir = 'input' # directory should be in current directory
outdir = 'FilteredAlignments_MedianSig%dMin%.2fWin%d' % (Zscore, cutoff, windsize) # output directory to be made in current working directory
missingchar = 'N' # X for AA and N for nucleotide (make sure it is uppercase)
excludelist = ['Homo_sapiens'] # taxa names to exclude from standard deviation calculation # , 'MMU', 'CJA'

print 'Reference Individual = ' + refname
print 'Zscore = %f' % (Zscore) # how many standard deviations away from mean 
print 'Pairwise deletion distance frequency cutoff value = %.2f' % (cutoff)
print 'Minimum AA alignment length accepted = %d' % (AAaligncutofflen)
print 'Minimum percent base coverage accepted = %d' % (PCcutoff)
print 'Window Size = %d' % (windsize)
print 'Window Step Value = %d' % (step)
print 'Missing data character = ' + missingchar

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()

def FormatDictionaryOfNucleotideSeqsToPhylip(seqs, Ns):
	longestseq = max(seqs.values(), key=len) # returns longest string in list
	longestname = max(seqs.keys(), key=len)
	o = '%d %d\n' % (len(seqs.keys()), len(longestseq))
	for n in Ns:
		nameout = n + (len(longestname)-len(n))*' '
		seqout = seqs[n] + (len(longestseq)-len(seqs[n]))*'N'
		o += '%s %s\n' % (nameout, seqout)
	return(o)

def FormatDictionaryOfNucleotideSeqsToFasta(d, Ns):  
	o = ''
	for k in Ns:
		o += '>%s\n%s\n' % (k, d[k])
	return(o)

def Return1IfValueIsGreaterThanCutoffForAllInDictionary(perccov, PCcutoff, countshit):
	boolean = 1 # gate open
	for n in perccov.keys():
		if perccov[n] < PCcutoff:
			boolean = 0 # if one individuals percent base coverage is below PCcutoff close gate.
			countshit[n] = countshit.get(n, 0) + 1
	return(boolean)

def CalculatePercentAACoverageForEachSequenceInDictionary(newseqs, lenseqs):
	perccov = {}
	for n in newseqs.keys():
		countbases = 0.0
		for char in newseqs[n]:
			if char != '-' and char != missingchar and char.isalpha():	
				countbases += 1.0
		perccov[n] = (countbases/lenseqs)*100.0
	return(perccov)
			
def MaskAlignedSeqsUsingStartAndEndCoordinates(seqs, start, end, missingchar, windsize ,lenseqs):
	for n in seqs.keys():
		try:
			start[n]
		except:
			pass
		else:
			for s,e in map(None, start[n], end[n]):
				mask = ''
				for char in seqs[n][s-1:e]:
					if char == '-':
						mask += char
					else:
						mask += missingchar
				seqs[n] = seqs[n][:s-1] + mask + seqs[n][e:]
	return(seqs)

def RecordValueIfGreaterThanCutoffValue(pdist, cutoff, i, windsize, start, end): # Not used
	for n in pdist.keys():
		if pdist[n] > cutoff: # if greater mask nucleotides in sequence of individual
			start[n] = start.get(n, []) + [i+1]
			end[n] = end.get(n, []) + [i+windsize]
	return(start, end)		

def RecordValueIfGreaterThanTwoStandardDevsAwayFromMean(mean, stddev, Zscore, cutoff, pdist, i, windsize, start, end):
	for n in pdist.keys():
		if pdist[n] > (mean+(Zscore*stddev)) and pdist[n] > cutoff: # value must be significantly different according to Zscore and also greater than pairwise distance cutoff.
			start[n] = start.get(n, []) + [i+1]
			end[n] = end.get(n, []) + [i+windsize]
	return(start, end)

def MakeNewReducedListAndExcludeKeyValues(d, excludelist):
	newdict = {}
	for key in d.keys():
		if key not in excludelist:
			newdict[key] = d[key]
	return(newdict)

def StandardDeviationUsingMedian(mylist):
	mn = Median(mylist)
	variance = sum([(e-mn)**2 for e in mylist]) / (len(mylist)-1) # average of the squared differences of values from mylist from their mean value. # minus 1 in denominator because window is sample of larger population (Bessel's correction)
	return(math.sqrt(variance)) # stddev = sigma = sqrt of variance

def StandardDeviation(mylist):
	mn = Mean(mylist)
	variance = sum([(e-mn)**2 for e in mylist]) / (len(mylist)-1) # average of the squared differences of values from mylist from their mean value. # minus 1 in denominator because window is sample of larger population (Bessel's correction)
	return(math.sqrt(variance)) # stddev = sigma = sqrt of variance

def Median(mylist):
	sorts = sorted(mylist)
	length = len(sorts)
	if not length % 2:
		return ((sorts[length / 2] + sorts[length / 2 - 1]) / 2.0)
	return sorts[length / 2]

def Mean(mylist):
	#print mylist
	return((sum(mylist)*1.0)/len(mylist)+1)

def CalculatePairwiseDeletionDistanceForAlignedSeqsComparedToRefSeq(seqs, refname):
	pdist = {}
	xlist = excludelist
	for n in seqs.keys():
		countgap = 0.0
		countmissing = 0.0
		countmismatch = 0.0
		denominator = 0.0
		for b,r in map(None, seqs[n], seqs[refname]): # b is query base, r is reference base
			if b == '-':
				countgap += 1.0
			elif b.upper() == missingchar: # X for AA and N for nucleotide
				countmissing += 1.0
			elif b.upper() != r.upper():
				countmismatch += 1.0
				denominator += 1.0
			else:
				denominator += 1.0
		if denominator != 0:
			pdist[n] = countmismatch / denominator
		else:
			pdist[n] = 0.0000001 # initiate very small Non-zero frequency value
			xlist = xlist + [n]
	return(pdist, xlist)

def SplitAlignedSeqsIntoWindows(seqs, numindivs, lenseqs, windsize, step, missingchar):
	start = {} # if window pdist is > cutoff record start and end positions of alignment for indiv
	end = {}
	for i in range(0, lenseqs, step):
		window = {}
		for n in seqs.keys(): # cut one window for all aligned sequences in seqs.
			window[n] = seqs[n][i:(i+windsize)]
		pdist, xlist = CalculatePairwiseDeletionDistanceForAlignedSeqsComparedToRefSeq(window, refname)
		redpdist = MakeNewReducedListAndExcludeKeyValues(pdist, xlist)
		if redpdist == {}:
			pass  # since this program is not supposed to mask gaps or missing data chars (with missingchar) then skip this window, i.e. do not mask anything if the window is ONLY gaps and missing data in all other species besides humans (i.e. ref sequence)
		else:
			#mean = Mean(redpdist.values())
			median = Median(redpdist.values())
			#stddev = StandardDeviationUsingMedian(redpdist.values())
			stddev = StandardDeviation(redpdist.values())
			start, end = RecordValueIfGreaterThanTwoStandardDevsAwayFromMean(median, stddev, Zscore, cutoff, redpdist, i, windsize, start, end) # we use median instead of mean because less influence from outliers
			#start, end = RecordValueIfGreaterThanCutoffValue(pdist, cutoff, i, windsize, start, end)
	newseqs = MaskAlignedSeqsUsingStartAndEndCoordinates(seqs, start, end, missingchar, windsize ,lenseqs)
	return(newseqs)

def ParsePhylip(fle): ### NOT USED ###
	seqs = {} # names keys, sequences value
	Ns = []	# names
	Ss = [] # Sequences
	FILE = open(fle, 'r')
	first = FILE.readline()
	numindivs, lenseqs = int(first.strip().split()[0]), int(first.strip().split()[1])
	if numindivs == 0 or lenseqs < AAaligncutofflen:
		return('empty', 'empty', 'empty', 'empty', 'empty')
	else:
		line = FILE.readline()
		while line:
			n, s = line.strip().split()[0], line.strip().split()[1]
			seqs[n] = s
			Ns.append(n)
			Ss.append(s)
			line = FILE.readline()
		FILE.close()
		return(seqs, Ns, Ss, numindivs, lenseqs)

def ParseFasta(fle):
	seqs = {}
	Ns = []
	Ss = []
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip()
			seqs[n] = ''
			Ns.append(n)
		else:
			seqs[n] += line.strip()
		line = FILE.readline()
	FILE.close()
	if len(Ns) == 0:
		return('empty', 'empty', 'empty', 'empty', 'empty')
	Ss = [seqs[n] for n in Ns]
	if len(Ss[0]) < AAaligncutofflen:
		return('empty', 'empty', 'empty', 'empty', 'empty')
	else:	
		return(seqs, Ns, Ss, len(Ns), len(Ss[0]))

def main():
	os.system('mkdir %s' % (outdir))
	
	cwd = os.getcwd()
	path = '%s/%s/*%s' % (cwd, phybindir, fileextension)
	filelist = glob.glob(path)

	countshit = {}
	outphyliplist = []
	countperccovcutoff = 0
	countemptyalignments = 0
	for f in filelist: # f is complete path to file with filename
		stump = f.split('/')[-1][:-len(fileextension)]
		seqs, Ns, Ss, numindivs, lenseqs = ParseFasta(f)
		if seqs != 'empty':
			newseqs = SplitAlignedSeqsIntoWindows(seqs, numindivs, lenseqs, windsize, step, missingchar)
			perccov = CalculatePercentAACoverageForEachSequenceInDictionary(newseqs, lenseqs)
			boolean = Return1IfValueIsGreaterThanCutoffForAllInDictionary(perccov, PCcutoff, countshit)
			if boolean == 1:
				output = FormatDictionaryOfNucleotideSeqsToFasta(newseqs, Ns)
				o = stump + '_aligned.filter.fa'
				outfile = '%s/%s/%s' % (cwd, outdir, o)
				WriteOUT(outfile, output)
				outphyliplist.append(outfile)
			elif boolean == 0:
				countperccovcutoff += 1
		elif seqs == 'empty':
			countemptyalignments += 1

	print 'Number of left after excluding empty alignments: %d' % (len(outphyliplist))
	print 'Number of excluded empty alignments (i.e. gblocks removed all bases) or Alignment was < %d AA long: %d' % (AAaligncutofflen, countemptyalignments)
	print 'Number of excluded alignments because alignment had < %d percentage of bases covered: %d' % (PCcutoff, countperccovcutoff)
	for n in countshit.keys():
		print '%s: %d' % (n, countshit[n])

main()
