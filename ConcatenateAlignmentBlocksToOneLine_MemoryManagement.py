import datetime
import os

alignmentfile = 'chr22_version1.aln.n1000000.v2.txt'  # test with chr22_version1.aln.n1000000.txt
speciesfile = 'Species251.txt'
referencespecies = 'Homo_sapiens'
linelimit = 10000


def concatenate_files(filelist, specieslist):
    # assuming non-interleaved fasta file
    # assuming sequences have same number of species and same length already
    print('Concatenating species from all %d files, one species at a time' % (len(filelist)))
    outfile = '.'.join(filelist[0].split('.')[:-2] + ['Concat%dFiles' % (len(filelist)), 'fa'])
    with open(outfile, 'w') as OUT:  # clear file before appending
        OUT.write('')
    with open(outfile, 'a') as OUT:
        print('Concatenating sequences to file: %s' % outfile)
        for species in specieslist:
            print('Starting to concatenate species: %s' % species)
            print(datetime.datetime.now())
            species_seq = []  # will contain all alignment fragments for one species # and free memory
            for f in filelist:
                gate = 'open'  # open gate at beginning of each file
                with open(f, 'r') as FILE:
                    line = FILE.readline()
                    while line:
                        if line.strip() == '':  # if the last line is empty skip it
                            pass
                        elif gate == 'closed':  # stop reading file when already read species of interest
                            break
                        elif line[0] == '>':
                            if line.strip()[1:] == species:
                                species_seq.append(FILE.readline().strip())
                                gate = 'closed'  # close gate when species sequence for this file has been appended
                        line = FILE.readline()

            OUT.write('>%s\n%s\n' % (species, ''.join(species_seq)))
            print('Finished concatenating species: %s, Sequence length = %d' % (species, len(''.join(species_seq))))
            print(datetime.datetime.now())
    for f in filelist:
        os.remove(f)


def main(alignmentfile, speciesfile, referencespecies):
    ''' filename (full path) '''

    # linelimit = 5000000

    print(alignmentfile)
    print(speciesfile)
    print(referencespecies)
    print(datetime.datetime.now())
    # read in all species from species file, one species name per line
    with open(speciesfile, 'r') as FILE:
        specieslist = [line.strip() for line in FILE]

    filelist = []
    # d is a dictionary number of keys = number of species and values (will) = all concatenated seq blocks
    d = {species: [] for species in specieslist}  # initiate keys for all species in species list file
    d_coords = {species: 1 for species in specieslist}  # tracks the last long N extension when missing from prev block
    count = 1
    filecount = 1
    gate = 'closed'  # controls when data should be written to temp files
    print('Starting to reformat alignment blocks.')
    with open(alignmentfile, 'r') as FILE:
        line = FILE.readline()
        while line:
            if (line[0] == '>') and (line.split('-')[0] == '>' + referencespecies):
                # >Homo_sapiens-chr22(+)/10509923-10509970
                # for humans only extract the start and end coordinates
                start, end = [int(i) for i in line.strip().split('/')[1].split('-')]
                # extend all sequences by some number of N's if the human start coordinate for the present block is
                # some distance from the previous block (start - end) 1001 - 1001 = 0, 1005 - 1001 = 4
                # the end coordinates are +1 from the actual end base for the block (seq is 1000 bases, end is 1001)
                species = line.split('-')[0][1:]  # species name
                line = FILE.readline()  # read next line
            elif line[0] == '>':
                # for all other species
                species = line.split('-')[0][1:]  # species name
                line = FILE.readline()  # read next line
            # here you add the N extension to the beginning of the seq block (if needed as a spacer)
            seq = 'N'*(start - d_coords[species])  # also corrects for if species was missing from all prev blocks
            d_coords[species] = end
            if line[0] != '>':
                # should only enter on line immediately after name lines
                # gather the sequence (which could be interleaved) to non-interleaved
                while '>' not in line:  # can not use line[0] != '>' if the last line is empty..
                    seq += line.strip()
                    line = FILE.readline()  # it will break when it is at a name line where line[0] == '>'
                    if line.strip() == '':
                        break
                # add the sequence for each species
                # creates keys = species name (if not already made)
                # concatenates sequences
                # if concatenating strings (slow) d[species] = d.setdefault(species, '') + seq
                d.setdefault(species, []).append(seq)
            # chrom 22 is about 51 million bps long which is about 71 Gb in size
            # if end > 10509971:
                # break
            if count % linelimit == 0:
                gate = 'open'
            if (gate == 'open') and (line.split('-')[0] == '>' + referencespecies):
                for species in specieslist:
                    if d_coords[species] < d_coords[referencespecies]:
                        d.setdefault(species, []).append('N' * (d_coords[referencespecies] - d_coords[species]))
                with open('%s.reformat.%d.fa' % (alignmentfile, filecount), 'w') as OUT:
                    # + '\n' just add a new line to the end of the file
                    OUT.write('\n'.join(['>%s\n%s' % (species, ''.join(d[species])) for species in specieslist]) + '\n')
                d = {species: [] for species in specieslist}  # initiate keys for all species in species list file
                d_coords = {species: end for species in specieslist}
                print('%s.reformat.%d.fa' % (alignmentfile, filecount))
                filelist.append('%s.reformat.%d.fa' % (alignmentfile, filecount))
                filecount += 1
                gate = 'closed'
            count += 1

    for species in specieslist:
        if d_coords[species] < d_coords[referencespecies]:
            d.setdefault(species, []).append('N' * (d_coords[referencespecies] - d_coords[species]))

    with open('%s.reformat.%d.fa' % (alignmentfile, filecount), 'w') as OUT:
        # + '\n' just add a new line to the end of the file
        OUT.write('\n'.join(['>%s\n%s' % (species, ''.join(d[species])) for species in specieslist]) + '\n')
    print('%s.reformat.%d.fa' % (alignmentfile, filecount))
    filelist.append('%s.reformat.%d.fa' % (alignmentfile, filecount))

    d = None  # clear memory

    print('Finished reformatting alignment blocks.')

    concatenate_files(filelist, specieslist)


main(alignmentfile, speciesfile, referencespecies)  # test with align file chr22_version1.aln.n1000000.txt



