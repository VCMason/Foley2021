''' Author: Victor C Mason '''
''' Date: May 14, 2015 '''
'''Editor: Andrew Harris'''
'''Date Modified: 6/22/2020'''
''' Output: '''
''' Summary: v14.1 only requires the required individuals from the hypothesis file to be present to report the INDEL. '''
''' Summary: v13.7 requires that the number of individuals with the INDEL must be >= number of individuals from hypothesis minus number of individuals with missing data from hypothesis that overlaps the INDEL. '''
# example hypothesis file
# before | is possible species to have deletion, after | are REQUIRED species to have deletion
# once the program sees a '#' at the beginning of a line it will stop reading additional hypotheses, breaks at '#####'
# >EleManHyr
# Loxodonta_africana, Trichechus_manatus, Heterohyrax_brucei | Loxodonta_africana, Trichechus_manatus, Heterohyrax_brucei
# >EleMan
# Loxodonta_africana, Trichechus_manatus | Loxodonta_africana, Trichechus_manatus
# >EleHyr
# Loxodonta_africana, Heterohyrax_brucei | Loxodonta_africana, Heterohyrax_brucei
# >ManHyr
# Trichechus_manatus, Heterohyrax_brucei | Trichechus_manatus, Heterohyrax_brucei
#####
# >Euarchonta
# Homo, Pan, Gorilla, Pongo, Nomascus, Papio, Macaca, Callithrix, Tarsius, Otolemur, Microcebus, GVA, Tupaia, TCH | GVA, TCH, Homo
#

import os
import re
import timeit
import shutil
import logging
import pprint

start = timeit.default_timer()

# Inititate logging file
logging.basicConfig(filename='./Victors_script_info.log', level=logging.DEBUG, filemode='w', format='')

hypotheses = 'Chiroptera_HypothesesList.txt'
speciesfile = 'species241.txt'
min = 2  # lowest number of taxa allowed to grant support to a hypothesis
minlen = 10 # minimum length of INDELs counted or reported
missminlen = 0  # minimum length of missing data 1character stretches considered by program
missdatachar = '-'  # use capital letter
# taxamisscutoff = the number of taxa included in the hypothesis allowed to have missing data overlapping a shared deletion
taxamisscutoff = 0
fileextension = '.fas'  # '.aln.reformat.Concat251Files.fa'  # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

# ----------------------------------------------------------------------
# New parameters (AJH - 6/22/2020)
N_CHECK_LENGTH = 100 # Number of bases to check around deletions for N's
N_COUNT_THRESH = 0



# NOTE: !!! ONLY TAXA IN outputtaxa LIST WILL BE logging.debugED OUT TO HIGHLIGHTED ALIGNED FASTA FILES FOR HIGHLIGHTED INDELs
# read in species names
with open(speciesfile, 'r') as FILE:
    outputtaxa = [line.strip() for line in FILE]


### not used ###
# notallowedalone = []  # not allowed to represent a shared indel with only these two taxa  # 'TCH', 'Tupaia'


def WriteOUT(outfile, output):
    OUT = open(outfile, 'w')
    OUT.write(output)
    OUT.close()


def FormatDictionaryOfNucleotideSeqsToFasta(d, orderoutputtaxa):
    o = ''
    for k in orderoutputtaxa:  # d.keys()
        try:
            d[k]
        except:
            pass
        else:
            o += '>%s\n%s\n' % (k, d[k])
    return o


def FormatOutput(d, f):  # list of dictionaries, output, fasta(ish) name (in this case gene name)
    o = ''
    o += '>%s\n' % (f)
    for key in sorted(d.keys(), key=len)[::-1]:  # sort keys by length and reverse so decending in length
        o += '%d\t%s\n' % (d[key], key)
    return o


def FormatOutputLists(listdict, listnames):  # dictionary, output, fasta(ish) name (in this case gene name)
    o = ''
    count = 0
    for d, g in zip(listdict, listnames):  # python 2.7 map(None, listdict, listnames)
        o += '>%s\n' % (g)
        for key in sorted(d.keys(), key=len)[::-1]:
            o += '%s\t%s\n' % (key, ','.join(d[key]))
            count += 1
    return o


def N_flank_check(seq, start, end, outgroup=False):
    # These variables are the flanking sequences of given length by "n_check_length"
    start_flank_seq = seq[start-N_CHECK_LENGTH:start]
    end_flank_seq = seq[end: end + N_CHECK_LENGTH]
    if outgroup:
        # First check for N's within the indel span
        deletion_seq = seq[start:end]
        if deletion_seq.count("N") > N_COUNT_THRESH:
            # print("Deletion seq(Outgroup) Fail:", deletion_seq)
            return False
        elif start_flank_seq.count("N") > N_COUNT_THRESH:
            # print("Start flank(Outgroup) Fail:",start_flank_seq)
            return False
        elif end_flank_seq.count("N") > N_COUNT_THRESH:
            # print("End flank(Outgroup) Fail:",end_flank_seq)
            return False
        else:
            # print("deletion seq:", deletion_seq)
            return True
    else:
        if start_flank_seq.count("N") > N_COUNT_THRESH:
            # print("Start flank Fail:",start_flank_seq)
            return False
        elif end_flank_seq.count("N") > N_COUNT_THRESH:
            # print("End flank Fail:",end_flank_seq)
            return False
        else:
            return True


def FormatGeneFastasToHighlightINDELs(seqs, indelspergene, outgrouptaxa):
    listoutseqs = []
    starts = []
    ends = []
    """
    AJH - 6/25/2020
    Changes:
        - New variable "rejected_deletions".
        - "rejected_deletions" is new Return for function.
    """
    rejected_deletions = []
    # logging.debug indelspergene.values()
    for listBS in indelspergene.values():
        for span in listBS:
            # logging.debug span
            # logging.debug type(span)
            starts.append(int(span.split('_')[0]))
            ends.append(int(span.split('_')[1]))

    ssort = sorted(starts)
    esort = [x for (y, x) in sorted(zip(starts, ends))]  # sort ends based on starts

    count = 0
    for s, e in zip(ssort, esort):  # access sorted list of start and end coordinates for this gene
        outseqs = {}  # only make this shit around one INDEL at a time, so if there are three indels in a gene then you
        # will have three files of the same gene with one of the three indels highlighted in each gene file.
        # print("start:", s, "End:", e)  
        for n in seqs.keys():
            # trim the output sequences to those flanking the indel 
            outseqs[n] = seqs[n][s - 1000:e + 1000]
        # output the whole gene file with indel highlighted with missdatachar, but maybe change to unique
        # outseqs[n] = seqs[n][:s-10] + missdatachar*10 + seqs[n][s-10:e+10] + missdatachar*10 + seqs[n][e+10:]
        listoutseqs.append(outseqs)
        count += 1

    # ssort = sorted(s)[::-1]
    # esort = [x for (y,x) in sorted(zip(s,e))][::-1] # sort e based on s, then reverse list
    return listoutseqs, rejected_deletions



def RecordINDELsThatSupportEachHypothesisAccountForMissData(htaxa, hreqtaxa, outreqtaxa, sharedspans, overlapmiss, seqs, taxamisscutoff=0):
    indelspergene = {}  # key is taxa combinations, value is spans (start_end) coordinates of INDELs
    indelcount = 0
    indelcountexcluded = 0
    for s in sharedspans.keys():
        # this is basically the gateway to hypothesis heaven.
        # already checked len(sharedspans[s]) >= min so don't need that now
        # requires that there are fewer species in hypothesis than in the alignment len(sharedspans[s]) < len(seqs)
        # requires that species required for hyp are present in shared span (deletion)
        # requires that species  with shared span (deletion) are allowed to have deletion (h[hyp])
        # nummiss will equal the number of taxa from the hypothesis file (not required taxa)
        # that have missing data overlapping with the shared span (deletion)

        ############ requires that the shared DEL be present in less than numtaxa - 1.
        ############ and len(list(set(taxlist) & set(sharedspans[s]))) >= 1
        ############ removed: and len(list(set(sharedspans[s]) ^ set(notallowedalone))) >= 1
        ############ removed: len(sharedspans[s]) >= (len(hreq[hyp]))

        nummiss = 0
        try:
            overlapmiss[s]
        except:
            pass
        else:
            nummiss = len(list(set(overlapmiss[s]) & set(outreqtaxa)))

        start = int(s.split("_")[0])
        end = int(s.split("_")[1])
        

        if len(sharedspans[s]) < len(seqs) and set(hreqtaxa) <= set(sharedspans[s]) and \
                set(sharedspans[s]) <= set(htaxa) and set(hreqtaxa) <= set(sharedspans[s]) \
                and nummiss <= taxamisscutoff:
            # INDEL supports hypothesis so record it.
            """
            First check that all DEL sequences pass N_check, ignore span N_flank_check fails
            """
            break_check = False
            for samp in htaxa:
                if not N_flank_check(seqs[samp], start, end):
                    break_check = True
                    continue
            for samp2 in hreqtaxa:
                if not N_flank_check(seqs[samp2], start, end):
                    break_check = True
                    continue

            for samp3 in outreqtaxa:
                if not N_flank_check(seqs[samp3], start, end, outgroup=True):
                    break_check = True
                    continue
            if break_check:
                indelcountexcluded += 1
                continue
            indelcount += 1
            # change our set back into list to sort alphabetically before
            k = '_'.join(sorted(list(sharedspans[s])))
            indelspergene[k] = indelspergene.setdefault(k, []) + [s]
        elif len(sharedspans[s]) < len(seqs) and set(hreqtaxa) <= set(sharedspans[s]) and \
                set(sharedspans[s]) <= set(htaxa) and set(hreqtaxa) <= set(sharedspans[s]) \
                and nummiss > taxamisscutoff:
            # INDEL supports hypothesis so record it, however, there are taxa from the hypothesis file that
            # have missing data overlapping the deletion
            indelcountexcluded += 1

    return indelcount, indelcountexcluded, indelspergene


def RecordINDELsThatSupportEachHypothesis(htaxa, hreqtaxa, sharedspans, seqs):
    indelspergene = {}  # key is taxa combinations, value is spans (start_end) coordinates of INDELs
    indelcount = 0
    for s in sharedspans.keys():
        # this is basically the gateway to hypothesis heaven.
        # already checked len(sharedspans[s]) >= min so don't need that now
        # requires that there are fewer species in hypothesis than in the alignment
        # requires the shared DEL be present more than or equal to number of required taxa from hyp (hreq).
        # requires that the shared DEL be present in less than numtaxa - 1.
        # and len(list(set(taxlist) & set(sharedspans[s]))) >= 1
        # removed: and len(list(set(sharedspans[s]) ^ set(notallowedalone))) >= 1
        if len(sharedspans[s]) < len(seqs) and len(sharedspans[s]) >= (len(hreqtaxa)) and \
                set(sharedspans[s]) <= set(htaxa) and set(hreqtaxa) <= set(sharedspans[s]):
            indelcount += 1  # INDEL supports hypothesis so record it.
            # change our set back into list to sort alphabetically before
            k = '_'.join(sorted(list(sharedspans[s])))
            indelspergene[k] = indelspergene.setdefault(k, []) + [s]

    return indelcount, indelspergene


def TrimSharedSpansToHypothesisTaxa(htaxa, hreqtaxa, sharedspans, seqs):
    trimsharedspans = {}  # is a container for dictionary indelspergene, key is hyp, value is dictionary indelspergene
    for s in sharedspans.keys():
        # this is basically the gateway to hypothesis heaven.
        # already checked len(sharedspans[s]) >= min so don't need that now
        # requires that there are fewer species in hypothesis than in the alignment
        # requires the shared DEL be present more than or equal to number of required taxa from hyp (hreq).
        # requires that the shared DEL be present in less than numtaxa - 1.
        # and len(list(set(taxlist) & set(sharedspans[s]))) >= 1
        # removed: and len(list(set(sharedspans[s]) ^ set(notallowedalone))) >= 1

        if len(sharedspans[s]) < len(seqs) and len(sharedspans[s]) >= (len(hreqtaxa)) and \
                set(sharedspans[s]) <= set(htaxa) and set(hreqtaxa) <= set(sharedspans[s]):
            trimsharedspans[s] = sharedspans[s]

    return trimsharedspans


def CollectAllContinuousStringsOfCharFromDict_Span(char, d, minlen=1, useminlen=False):
    spans = {}
    for n in d.keys():  # record spans of all indels from all taxa
        matchObjs = [m for m in re.finditer(r'%s+' % (char), d[n])]  # record all indel events
        for m in matchObjs:
            indellength = m.span()[1] - m.span()[0]
            span_start = m.span()[0]
            span_end = m.span()[1]
            if useminlen == True:
                if indellength >= minlen:  # minlen is an integer = minimum length of INDEL
                    s = '%d_%d' % (m.span()[0], m.span()[1])
                    spans[s] = spans.get(s, []) + [n]
            else:
                s = '%d_%d' % (m.span()[0], m.span()[1])
                spans[s] = spans.get(s, []) + [n]
    return (spans)


def RecordMissDataOverlapingINDELsFromSequenceDict(d, sharedspans, missdatachar, outgrouptaxa, missminlen):
    # outgrouptaxa = list of outgroup taxa from hyp file. I will find miss data character spans overlapping shared spans
    # sharedpsans == all spans(INDELs) ('start_end') as key and all taxa that have that span(INDEL) as list in value.
    # records all spans(of missdatachar) ('start_end') as key
    # and all taxa that have that span(of missdatachar) as list in value.
    # here i limit all sequences to only those specified in the outgroup section of the hypothesis file
    # this is require outgroup taxa to have sequence overlapping deletions (and not have missing data)
    reducedseqs = {species: d[species] for species in outgrouptaxa}
    # print(reducedseqs)
    spanmiss = CollectAllContinuousStringsOfCharFromDict_Span(missdatachar, reducedseqs, missminlen, useminlen=True)
    logging.debug('Number of missing character spans (%s): %d' % (missdatachar, len(spanmiss)))
    # indel span as key, and value is set()tuple of all taxa that have a missingdatachar inside of indel span.
    overlapmiss = {}
    # identify which taxa have missing data within the range(span) of the INDELs identified above.
    for mspan in spanmiss.keys():
        ms, me = int(mspan.split('_')[0]), int(mspan.split('_')[1])
        # print(ms, me)
        for spanindel in sharedspans.keys():
            ins, ine = int(spanindel.split('_')[0]), int(spanindel.split('_')[1])
            # if we assume that  ranges are well-formed (so that x1 <= x2 and y1 <= y2) then
            # x1 <= y2 && y1 <= x2
            # meaning the aligned sequence overlaps in some way to the gff3 feature
            if (ins <= me-1) and (ms <= ine-1):
                overlapmiss[spanindel] = spanmiss[mspan]  # record list of all taxa with missing data for this span
    for k in overlapmiss.keys():
        overlapmiss[k] = set(overlapmiss[k])  # just making sure that the list of taxa is unique.
    logging.debug('Number sharedspans: %d. Number of shared spans that overlap missing data character: %d' % (len(sharedspans), len(overlapmiss)))
    return overlapmiss


def RecordIndelsFromDict(d, minlen):
    # records all spans(INDELs) ('start_end') as key and all taxa that have that span(INDEL) as list in value.
    spans = CollectAllContinuousStringsOfCharFromDict_Span('-', d, minlen, useminlen=True)
    logging.debug('Number of gap \'-\' spans: %d' % (len(spans)))
    sharedspans = {}
    for s in spans.keys():
        spans[s] = set(spans[s])  # make list of taxa a set of unique taxa (as value) that has INDEL (span) as key
        if len(spans[s]) >= min:  # then it means the span(INDEL) is shared with at least two taxa.
            # save set of taxa to sharedspans if there is more than one taxa in value spans[s]
            sharedspans[s] = spans[s]

    logging.debug('Number sharedspans: %d.' % (len(sharedspans)))
    return sharedspans


### NOT USED, but correct ###
def ReplaceGapsInFrontAndBackOfEachSequenceWithMissDataCharFromDict(seqs):  ### NOT USED, but correct ###
    Xseqs = {}
    for key in seqs.keys():

        Obj = re.search(r'(^-+[^-]).+', seqs[key])
        try:
            Obj.group(1)
        except:
            frontfix = seqs[key]
        else:
            front = len(Obj.group(1)) - 1
            # minus 1 because matching NON dash character too in regex.
            frontfix = re.sub(r'^-+', missdatachar * front, seqs[key])

        Obj = re.search(r'.+([^-]-+$)', seqs[key])
        try:
            Obj.group(1)
        except:
            Xseqs[key] = frontfix
        else:
            back = len(Obj.group(1)) - 1
            Xseqs[key] = re.sub(r'-+$', missdatachar * back, frontfix)
    return Xseqs  ### NOT USED, but correct ###


def ParseFasta(fle):
    # returns dictionary of sequences, names as key (without >), seqs as values
    seqs = {}
    names = []
    FILE = open(fle, 'r')
    line = FILE.readline()
    while line:
        if line[0] == '>':
            n = line[1:].strip()
            names.append(n)
            seqs[n] = ''
        elif line.strip() == '#####':  # for hypotheses file
            break
        else:
            seqs[n] += line.strip()
        line = FILE.readline()
    FILE.close()
    return seqs, names


def main(taxamisscutoff, minlen, missminlen):
    filelist = []
    files = filter(os.path.isfile, os.listdir('.'))
    for filename in files:
        if extension.match(filename) != None:
            filelist.append(filename)

    logging.debug('Reading in fasta file: %s' % (hypotheses))
    hyps, hypnames = ParseFasta(hypotheses)
    logging.debug('Finished reading in fasta file: %s' % (hypotheses))
    h = {}  # will taxa that represent various hypotheses
    hreq = {}
    outreq = {}  # dictionary, each key is hypothesis name, each value is list of species
    outdirs = []
    for key in hyps.keys():
        # Reads in hypotheses from .fasta(ish) file and has all taxa as list. i.e. >Sundatheria\nGVA, TCH, Tupaia\n
        h[key] = ''.join(hyps[key].split('|')[0].split()).split(',')
        # Reads in Taxa defined to be required for that hypothesis to be valid.
        hreq[key] = ''.join(hyps[key].split('|')[1].split()).split(',')
        # Reads in outgroup taxa that is required to have no missing data under specified conditions
        outreq[key] = ''.join(hyps[key].split('|')[2].split()).split(',')
        outdir = '_'.join([key, 'DEL', 'MinLen%d' % minlen, 'MinTaxa%d' % min, 'MissMinLen%d' % missminlen])
        logging.debug('Hypothesis: %s' % key)
        logging.debug('Outgroup taxa: %s' % outreq[key])
        logging.debug('Required taxa to have shared deletion: %s' % hreq[key])
        logging.debug('Possible taxa to have shared deletion: %s' % h[key])
        logging.debug('Output directory: %s' % outdir)
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        else:
            shutil.rmtree(outdir)
            os.makedirs(outdir)
        outdirs.append(outdir)

    oneindeldict = {}  # counts all INDELs that support each hypothesis
    hypallgene = {}  # dictionary container of all genes, key is hypothesis name
    hypindelspergene = {}  # is a container for dictionary indelspergene, key is hyp, value is dictionary indelspergene
    hypallindelspergene = {}  # dictionary container of all indels per gene, key is hypothesis name
    allgene = []
    allindelspergene = []
    # sort files by filename because of the number prefix to keep them in order [00001_*, 00002_*, ... 07118_*]
    for fle in sorted(filelist):
        gene = fle  # makes chromosome name from filename
        allgene.append(gene)
        logging.debug('Reading in fasta file: %s' % fle)
        seqs, seqnames = ParseFasta(fle)
        logging.debug('Finished reading in fasta file: %s. Number of sequences: %d' % (fle, len(seqs)))
        sharedspans = RecordIndelsFromDict(seqs, minlen)  # do this once, it's kinda slow
        for hyp in h.keys():
            #outgrouptaxa = outreq[hyp]
            logging.debug('******************************************************************\n')
            logging.debug('Hypothesis: %s' % hyp)
            logging.debug('Outgroup taxa: %s' % outreq[hyp])
            logging.debug('Required taxa to have shared deletion: %s' % hreq[hyp])
            logging.debug('Possible taxa to have shared deletion: %s' % h[hyp])
            if outreq[hyp] == ['NA']:
                # there are no outgroup taxa specified, therefore we do not need to account for missing data
                logging.debug('Recording indels from sequences that support hypothesis: %s' % hyp)
                indelcount, indelspergene = RecordINDELsThatSupportEachHypothesis(h[hyp], hreq[hyp], sharedspans, seqs)
                oneindeldict[hyp] = oneindeldict.setdefault(hyp, 0) + indelcount
                hypindelspergene[hyp] = indelspergene
                oneindeldict[hyp + '_OverlappingMissingDataCharacters'] = 0
            elif outreq[hyp] != ['NA']:  # if there are outgroup taxa specified
                # overlapmiss contains streches of missdatachar (from outgroup taxa) that overlap deletions in shared spans
                logging.debug('Trim all shared spans to only those possible for hypothesis')
                trimsharedspans = TrimSharedSpansToHypothesisTaxa(h[hyp], hreq[hyp], sharedspans, seqs)
                logging.debug('Finding missing data in outgroup taxa that overlap shared deletions')
                overlapmiss = RecordMissDataOverlapingINDELsFromSequenceDict(seqs, trimsharedspans, missdatachar, outreq[hyp], missminlen)
                logging.debug('Recording indels from sequences that support hypothesis: %s and account for missing data in outgroup taxa' % hyp)
                indelcount, indelcountexcluded, indelspergene = RecordINDELsThatSupportEachHypothesisAccountForMissData(h[hyp], hreq[hyp], outreq[hyp], trimsharedspans, overlapmiss, seqs, taxamisscutoff)
                # print(oneindeldict.setdefault(hyp, 0) + indelcount)
                oneindeldict[hyp] = oneindeldict.setdefault(hyp, 0) + indelcount
                oneindeldict[hyp + '_OverlappingMissingDataCharacters_OrRejected'] = oneindeldict.setdefault(hyp + '_OverlappingMissingDataCharacters_OrRejected', 0) + indelcountexcluded
                hypindelspergene[hyp] = indelspergene
            logging.debug(oneindeldict)
            # hypindelspergene[hyp] == indelspergene
            hypallindelspergene.setdefault(hyp, []).append(hypindelspergene[hyp])
            logging.debug('Highlighting DELs')
            """
            This next portion of code will go through each sequence and will format the output fasta sequences
            and will remove anything not being used.
            """
 # ------------------------------------------------------------------------------------------------------------------------         
            # if gene has an indel that supports a hypothesis then copy it to directory outdir
            if hypindelspergene[hyp] != {}:
                logging.debug('Found shared DELETIONS in alignment %s for hyp %s' % (fle, hyp))
                listofDicts, rejected_dels = FormatGeneFastasToHighlightINDELs(seqs, hypindelspergene[hyp], outreq[hyp])
                # remove all possible taxa from output list
                trimoutputtaxa = [ele for ele in outputtaxa if (ele not in h[hyp]) and (ele not in outreq[hyp])]
                # remove required taxa from possible taxa
                allowedbutnotrequired = [ele for ele in h[hyp] if ele not in hreq[hyp]]
                # order the taxa for output
                orderoutputtaxa = hreq[hyp] + allowedbutnotrequired + outreq[hyp] + trimoutputtaxa
                #logging.debug(orderoutputtaxa)
                logging.debug('Writing highlighted INDELs to files')
                for count, d in enumerate(listofDicts):
                    outdir = '_'.join([hyp, 'DEL', 'MinLen%d' % minlen, 'MinTaxa%d' % min, 'MissMinLen%d' % missminlen])
                    output = FormatDictionaryOfNucleotideSeqsToFasta(d, orderoutputtaxa)
                    dest = os.path.join(os.getcwd(), outdir, '%s_%d.fa' % (fle[:-len(fileextension)], count))
                    WriteOUT(dest, output)
                    logging.debug('Files written to: %s' % dest)
            else:
                logging.debug('No shared DELETIONS found in alignment %s for hyp %s' % (fle, hyp))

    logging.debug('Number Alignment files Queried: %d' % (len(filelist)))

    for hyp in hypindelspergene.keys():
        outdir = '_'.join([hyp, 'DEL', 'MinLen%d' % minlen, 'MinTaxa%d' % min, 'MissMinLen%d' % missminlen])

        outfile = '_'.join(['NumberOfDELsThatSupportEachHypothesis', 'MinLen%d' % minlen, 'MinTaxa%d' % min])
        output = FormatOutput(oneindeldict, outfile)

        dest = os.path.join(os.getcwd(), outdir, 'DELsSupporting%sCounts.v17.0.2.DelMinLen%dMinTaxa%d.txt' % (hyp, minlen, min))
        WriteOUT(dest, output)

        output = FormatOutputLists(hypallindelspergene[hyp], allgene)  # (allindelspergene, allgene)
        dest = os.path.join(os.getcwd(), outdir, 'DELsSupporting%s_PerGene.v17.0.2.DelMinLen%dMinTaxa%d.txt' % (hyp, minlen, min))
        WriteOUT(dest, output)


main(taxamisscutoff, minlen, missminlen)

stop = timeit.default_timer()

logging.debug(stop - start)
