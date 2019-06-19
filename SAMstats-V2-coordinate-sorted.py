##################################
#                                #
# Last modified 2019/06/18       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import gc
import pysam
import math
import string
import os
from sets import Set

# FLAG field meaning
# 0x0001 1 the read is paired in sequencing, no matter whether it is mapped in a pair
# 0x0002 2 the read is mapped in a proper pair (depends on the protocol, normally inferred during alignment) 1
# 0x0004 4 the query sequence itself is unmapped
# 0x0008 8 the mate is unmapped 1
# 0x0010 16 strand of the query (0 for forward; 1 for reverse strand)
# 0x0020 32 strand of the mate 1
# 0x0040 64 the read is the first read in a pair 1,2
# 0x0080 128 the read is the second read in a pair 1,2
# 0x0100 256 the alignment is not primary (a read having split hits may have multiple primary alignment records)
# 0x0200 512 the read fails platform/vendor quality checks
# 0x0400 1024 the read is either a PCR duplicate or an optical duplicate
# 0x0800 2048 supplementary alignment

def FLAG(FLAG):

    Numbers = [0,1,2,4,8,16,32,64,128,256,512,1024,2048]

    FLAGDict = {}

    MaxNumberList=[]
    for i in Numbers:
        if i <= FLAG:
            MaxNumberList.append(i)

    Residual=FLAG
    maxPos = len(MaxNumberList)-1

    while Residual > 0:
        if MaxNumberList[maxPos] <= Residual:
            Residual = Residual - MaxNumberList[maxPos]
            FLAGDict[MaxNumberList[maxPos]] = 1
            maxPos-=1
        else:
            maxPos-=1
  
    return FLAGDict

def run():

    if len(sys.argv) < 3:
        print 'usage: python %s BAMfilename chrom.sizes stats_out_filename' % sys.argv[0]
        print '\tthe script assumes the presence of NH tags'
        print '\t!!!do not run this script on a subset of the chromosomes, it assumes that it will see all alignments in the BAM file!!!'
        sys.exit(1)

    BAM = sys.argv[1]
    chrominfo = sys.argv[2]
    outfilename = sys.argv[3]

    chromInfoList = []
    chromInfoDict = {}
    linelist=open(chrominfo)
    for line in linelist:
        fields = line.strip().split('\t')
        chr = fields[0]
        start = 0
        end = int(fields[1])
        chromInfoList.append((chr,start,end))
        chromInfoDict[chr] = end

    samfile = pysam.Samfile(BAM, "rb" )

    try:
        print 'testing for NH tags presence'
        for alignedread in samfile.fetch():
            multiplicity = alignedread.opt('NH')
            print 'file has NH tags'
            break
    except:
        print 'no NH: tags in BAM file, exiting'
        sys.exit(1)

    ReadStatsDict = {}
    ReadStatsDict['proper_pairs'] = {}
    ReadStatsDict['proper_pairs']['all'] = {}
    ReadStatsDict['proper_pairs']['all'][1] = 0
#    ReadStatsDict['proper_pairs']['after_filtering'] = {}
#    ReadStatsDict['proper_pairs']['after_filtering'][1] = 0
    ReadStatsDict['not_proper_pairs'] = {}
    ReadStatsDict['not_proper_pairs']['all'] = {}
    ReadStatsDict['not_proper_pairs']['all'][1] = 0
#    ReadStatsDict['not_proper_pairs']['after_filtering'] = {}
#    ReadStatsDict['not_proper_pairs']['after_filtering'][1] = 0
    ReadStatsDict['unpaired'] = {}
    ReadStatsDict['unpaired']['all'] = {}
    ReadStatsDict['unpaired']['all'][1] = 0
#    ReadStatsDict['unpaired']['after_filtering'] = {}
#    ReadStatsDict['unpaired']['after_filtering'][1] = 0
    ReadStatsDict['paired_unaligned'] = 0
    ReadStatsDict['unpaired_unaligned'] = 0
#    ReadStatsDict['paired_unaligned_after_filtering'] = 0
#    ReadStatsDict['unpaired_unaligned_after_filtering'] = 0
    ReadStatsDict['spliced_reads'] = {}
    ReadStatsDict['spliced_reads']['all'] = {}
    ReadStatsDict['spliced_reads']['all'][1] = 0
#    ReadStatsDict['spliced_reads']['after_filtering'] = {}
#    ReadStatsDict['spliced_reads']['after_filtering'][1] = 0
    ReadStatsDict['unspliced_reads'] = {}
    ReadStatsDict['unspliced_reads']['all'] = {}
    ReadStatsDict['unspliced_reads']['all'][1] = 0
#    ReadStatsDict['unspliced_reads']['after_filtering'] = {}
#    ReadStatsDict['unspliced_reads']['after_filtering'][1] = 0

    ReadLengthDict = {}

    ComplexityDict = {}
    UP = 0.0
    UR = 0.0
    M0 = 0.0
    M1 = 0.0
    M2 = 0.0

    SeenDict = {}

    i=0
    samfile = pysam.Samfile(BAM, "rb" )
    for (chr,start,end) in chromInfoList:
        ComplexityChromosomesDict = {}
        try:
            jj=0
            for alignedread in samfile.fetch(chr, start, end):
                jj+=1
                if jj==1:
                       break
        except:
            print 'problem with region:', chr, start, end, 'skipping'
            continue
        for alignedread in samfile.fetch(chr, start, end):
            i+=1
            if i % 5000000 == 0:
                print str(i/1000000) + 'M alignments processed', chr,start,alignedread.pos,end
            fields=str(alignedread).split('\t')
            NH = alignedread.opt('NH')
            FLAGfields = FLAG(int(fields[1]))
            RL = 0
            for (m,bp) in alignedread.cigar:
                if m == 0:
                    RL += bp
            if ReadLengthDict.has_key(RL):
                pass
            else:
                ReadLengthDict[RL] = 0
            ReadLengthDict[RL] += 1./NH
            if FLAGfields.has_key(4):
                if FLAGfields.has_key(1):
                    if FLAGfields.has_key(8):
                        ReadStatsDict['paired_unaligned'] += 0.5
                else:        
                    ReadStatsDict['unpaired_unaligned'] += 1
                continue
            CIGAR = fields[5]
            if 'N' in CIGAR:
                if ReadStatsDict['spliced_reads']['all'].has_key(NH):
                    pass
                else:
                    ReadStatsDict['spliced_reads']['all'][NH] = 0
                ReadStatsDict['spliced_reads']['all'][NH] += 1./NH
            else:
                if ReadStatsDict['unspliced_reads']['all'].has_key(NH):
                    pass
                else:
                    ReadStatsDict['unspliced_reads']['all'][NH] = 0
                ReadStatsDict['unspliced_reads']['all'][NH] += 1./NH
            if FLAGfields.has_key(1) or FLAGfields.has_key(64) or FLAGfields.has_key(128):
                IsPairedEnd = True
            if IsPairedEnd:
                ProperPairs = True
                if FLAGfields.has_key(8):
                    ProperPairs = False
                if FLAGfields.has_key(64):
                    if NH == 1:
                        pos1 = int(fields[3])
                        pos2 = int(fields[7])
                        if FLAGfields.has_key(16):
                            pos1 = -pos1
                        else:
                            pos2 = -pos2
                if ProperPairs:
                    if ReadStatsDict['proper_pairs']['all'].has_key(NH):
                        pass
                    else:
                        ReadStatsDict['proper_pairs']['all'][NH] = 0
                    ReadStatsDict['proper_pairs']['all'][NH] += 0.5/NH
                    if NH == 1:
                        if FLAGfields.has_key(64):
                            if ComplexityChromosomesDict.has_key(pos1):
                                pass
                            else:
                                ComplexityChromosomesDict[pos1] = {}
                            if ComplexityChromosomesDict[pos1].has_key(pos2):
                                pass
                            else:
                                ComplexityChromosomesDict[pos1][pos2] = 0
                            ComplexityChromosomesDict[pos1][pos2] += 1
                            if ComplexityChromosomesDict[pos1][pos2] == 1:
                                UP += 1
                                UR += 1
                                M0 += 1
                                M1 += 1
                            elif ComplexityChromosomesDict[pos1][pos2] == 2:
                                UR += 1
                                M1 -= 1
                                M2 += 1
                            elif ComplexityChromosomesDict[pos1][pos2] == 3:
                                UR += 1
                                M2 -= 1
                            else:
                                UR += 1
                else:
                    if ReadStatsDict['not_proper_pairs']['all'].has_key(NH):
                        pass
                    else:
                        ReadStatsDict['not_proper_pairs']['all'][NH] = 0
                    ReadStatsDict['not_proper_pairs']['all'][NH] += 1./NH
            else:
                if NH == 1:
                    pos1 = int(fields[3])
                    if FLAGfields.has_key(16):
                        pos1 = -(pos1 + len(fields[9]))
                    pos2 = pos1
                    if ComplexityChromosomesDict.has_key(pos1):
                        pass
                    else:
                        ComplexityChromosomesDict[pos1] = {}
                    if ComplexityChromosomesDict[pos1].has_key(pos2):
                        pass
                    else:
                        ComplexityChromosomesDict[pos1][pos2] = 0
                    ComplexityChromosomesDict[pos1][pos2] += 1
                    if ComplexityChromosomesDict[pos1][pos2] == 1:
                        UP += 1
                        UR += 1
                        M0 += 1
                        M1 += 1
                    elif ComplexityChromosomesDict[pos1][pos2] == 2:
                        UR += 1
                        M1 -= 1
                        M2 += 1
                    elif ComplexityChromosomesDict[pos1][pos2] == 3:
                        UR += 1
                        M2 -= 1
                    else:
                        UR += 1
                if ReadStatsDict['unpaired']['all'].has_key(NH):
                    pass
                else:
                    ReadStatsDict['unpaired']['all'][NH] = 0
                ReadStatsDict['unpaired']['all'][NH] += 1./NH

    outfile=open(outfilename, 'w')

    outline='unique pairs, proper:\t' + str(int(ReadStatsDict['proper_pairs']['all'][1]))
    outfile.write(outline+'\n')
    outline='unique pairs, not proper:\t' + str(int(ReadStatsDict['not_proper_pairs']['all'][1]))
    outfile.write(outline+'\n')
    outline='unique unpaired reads:\t' + str(int(ReadStatsDict['unpaired']['all'][1]))
    outfile.write(outline+'\n')

    M = 0
    for RM in ReadStatsDict['proper_pairs']['all'].keys():
        if RM != 1:
            M += ReadStatsDict['proper_pairs']['all'][RM]
    outline='multiread pairs, proper:\t' + str(int(M))
    outfile.write(outline+'\n')

    M = 0
    for RM in ReadStatsDict['not_proper_pairs']['all'].keys():
        if RM != 1:
            M += ReadStatsDict['not_proper_pairs']['all'][RM]
    outline='multiread pairs, not proper:\t' + str(int(M))
    outfile.write(outline+'\n')

    M = 0
    for RM in ReadStatsDict['unpaired']['all'].keys():
        if RM != 1:
            M += ReadStatsDict['unpaired']['all'][RM]
    outline='multiread unpaired:\t' + str(int(M))
    outfile.write(outline+'\n')

    outline='unaligned pairs:\t' + str(int(ReadStatsDict['paired_unaligned']))
    outfile.write(outline+'\n')
    outline='unaligned unpaired reads:\t' + str(int(ReadStatsDict['unpaired_unaligned']))
    outfile.write(outline+'\n')
    outline='spliced unique reads:\t' + str(int(ReadStatsDict['spliced_reads']['all'][1]))
    outfile.write(outline+'\n')

    M = 0
    for RM in ReadStatsDict['spliced_reads']['all'].keys():
        if RM != 1:
            M += ReadStatsDict['spliced_reads']['all'][RM]
    outline='spliced multireads:\t' + str(int(M))
    outfile.write(outline+'\n')

    outline='unspliced unique reads:\t' + str(int(ReadStatsDict['unspliced_reads']['all'][1]))
    outfile.write(outline+'\n')

    M = 0
    for RM in ReadStatsDict['unspliced_reads']['all'].keys():
        if RM != 1:
            M += ReadStatsDict['unspliced_reads']['all'][RM]
    outline='unspliced multireads:\t' + str(int(M))
    outfile.write(outline+'\n')
    
    outline='Read Length, Minimum:\t'+str(min(ReadLengthDict.keys()))
    outfile.write(outline+'\n')
    outline='Read Length, Maximum:\t'+str(max(ReadLengthDict.keys()))
    outfile.write(outline+'\n')

    TotalReads=0.0
    TotalLength=0.0
    for length in ReadLengthDict.keys():
        TotalReads += ReadLengthDict[length]
        TotalLength += length*ReadLengthDict[length]
    outline = 'Read Length, Average:\t'+str(TotalLength/TotalReads)
    outfile.write(outline+'\n')

    outline='#Library complexity metrics (unique reads/pairs only):'
    outfile.write(outline+'\n')
    outline='U_P\t' + str(int(UP))
    outfile.write(outline+'\n')
    outline='U_R\t' + str(int(UR))
    outfile.write(outline+'\n')
    outline='M_0\t' + str(int(M0))
    outfile.write(outline+'\n')
    outline='M_1\t' + str(int(M1))
    outfile.write(outline+'\n')
    outline='M_2\t' + str(int(M2))
    outfile.write(outline+'\n')
    outline='PBC1\t' + str(M1/M0)
    outfile.write(outline+'\n')
    outline='PBC2\t' + str(M1/M2)
    outfile.write(outline+'\n')
    outline='NRF\t' + str(UP/UR)
    outfile.write(outline+'\n')

    outfile.close()

run()
