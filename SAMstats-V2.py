##################################
#                                #
# Last modified 2019/06/18       # 
#                                #
# Georgi Marinov                 #
#                                # 
##################################

import sys
import argparse 
import gc
import pysam
import math
import string
import os
#from sets import Set
import pandas as pd 

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

def parse_args():
    parser=argparse.ArgumentParser(description="This script computes mapping statistics for an input aligned bam file.")
    parser.add_argument("-SAMFileName",type=str,help="Unsorted BAM file, direct output of alignment algorithm. Use '-' to read from stdin")
    parser.add_argument("-batchSize",type=int,help="Number of reads for which alignments will be stored in memory. Set this value to 1 if it is certain that the alignments are sorted by read ID ")
    parser.add_argument("-statsOutFileName",default=None,type=str,help="File name to write summary statistics; Default is to print them to stdout")
    parser.add_argument("--addNH",action="store_true",default=False,help="adds NH tags to the output alignments")
    parser.add_argument("--doNotPrintAlignments",action="store_true",default=False)
    parser.add_argument("--minMAPQ",type=float,default=None,help="Minimum mapping quality. No quality filtering is performed if this argument is not provided.")
    parser.add_argument("--Ffilter",nargs="*",type=str,default=None,help="List of flags to filter out")
    return parser.parse_args()

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

def initializeOutputVars():
    """
    Initializes dictionaries to store read statistics
    Returns: ReadStatsDict,ReadLengthDict,ComplexityDict,ComplexityChromosomesDict,UP,UR,M0,M1,M2
    """
    ReadStatsDict = {}
    ReadStatsDict['proper_pairs'] = {}
    ReadStatsDict['proper_pairs']['all'] = {}
    ReadStatsDict['proper_pairs']['all'][1] = 0
    ReadStatsDict['proper_pairs']['after_filtering'] = {}
    ReadStatsDict['proper_pairs']['after_filtering'][1] = 0
    ReadStatsDict['not_proper_pairs'] = {}
    ReadStatsDict['not_proper_pairs']['all'] = {}
    ReadStatsDict['not_proper_pairs']['all'][1] = 0
    ReadStatsDict['not_proper_pairs']['after_filtering'] = {}
    ReadStatsDict['not_proper_pairs']['after_filtering'][1] = 0
    ReadStatsDict['unpaired'] = {}
    ReadStatsDict['unpaired']['all'] = {}
    ReadStatsDict['unpaired']['all'][1] = 0
    ReadStatsDict['unpaired']['after_filtering'] = {}
    ReadStatsDict['unpaired']['after_filtering'][1] = 0
    ReadStatsDict['paired_unaligned'] = 0
    ReadStatsDict['unpaired_unaligned'] = 0
    ReadStatsDict['paired_unaligned_after_filtering'] = 0
    ReadStatsDict['unpaired_unaligned_after_filtering'] = 0
    ReadStatsDict['spliced_reads'] = {}
    ReadStatsDict['spliced_reads']['all'] = {}
    ReadStatsDict['spliced_reads']['all'][1] = 0
    ReadStatsDict['spliced_reads']['after_filtering'] = {}
    ReadStatsDict['spliced_reads']['after_filtering'][1] = 0
    ReadStatsDict['unspliced_reads'] = {}
    ReadStatsDict['unspliced_reads']['all'] = {}
    ReadStatsDict['unspliced_reads']['all'][1] = 0
    ReadStatsDict['unspliced_reads']['after_filtering'] = {}
    ReadStatsDict['unspliced_reads']['after_filtering'][1] = 0
    
    ReadLengthDict = {}
    ComplexityDict = {}
    ComplexityChromosomesDict = {}
    UP = 0.0
    UR = 0.0
    M0 = 0.0
    M1 = 0.0
    M2 = 0.0    
    return ReadStatsDict,ReadLengthDict,ComplexityDict,ComplexityChromosomesDict,UP,UR,M0,M1,M2

def main():
    
    #Read in the arguments to the program from the argument parser
    args=parse_args()

    #determine whether the input comes from a SAM/BAM file or whether it should be read from stdin 
    SAM=args.SAMFileName
    if SAM=="-":
        SAM=sys.stdin
        
    BS=args.batchSize
    outfilename=args.statsOutFileName
    doAddNH=args.addNH
    #set optional quality and flag filtering 
    if args.minMAPQ is None:
        doMF=False
    else:
        doMF=True
        MAPQ=args.minMAPQ
        
    if args.Ffilter is None:
        doFF=False
    else:
        doFF=True
        FF=args.Ffilter
        
    doPrint= not args.doNotPrintAlignments
    
    #Initialize dictionaries and variables to store read summary statistics
    ReadStatsDict,ReadLengthDict,ComplexityDict,ComplexityChromosomesDict,UP,UR,M0,M1,M2=initializeOutputVars() 

    #parse through the input SAM file (or stdin), grabbing one batch of size BS at a time 
    for batch in pd.read_csv(SAM, header=None, chunksize=BS,sep='\t'):
        for index,row in batch.iterrows():
            if row[0].startswith('@'):
                if doPrint:
                    print(row)
                continue
            ID = row[0]
            if AlignmentDict.has_key(ID):
                pass
            else:
                AlignmentDict[ID] = [row]
            AlignmentDict[ID].append(line.strip())
        if len(AlignmentDict.keys()) == BS or EOF:
#            print NR, 'alignments processed'
            for ID in AlignmentDict.keys():
                IsPairedEnd = False
                fields = AlignmentDict[ID][0].split('\t')
                FLAGfields = FLAG(int(fields[1]))
                if FLAGfields.has_key(64) or FLAGfields.has_key(128):
                    IsPairedEnd = True
                if IsPairedEnd:
                    NH1 = 0
                    NH2 = 0
                    ProperPairs = True
                    for lline in AlignmentDict[ID]:
                        fields = lline.split('\t')
                        FLAGfields = FLAG(int(fields[1]))
                        if FLAGfields.has_key(64):
                            if not FLAGfields.has_key(4):
                                NH1 += 1
                            if NH1 == 1:
                                pos1 = int(fields[3])
                                if FLAGfields.has_key(16):
                                    pos1 = -(pos1 + len(fields[9]))
                        elif FLAGfields.has_key(128):
                            if not FLAGfields.has_key(4):
                                NH2 += 1
                            if NH2 == 1:
                                pos2 = int(fields[3])
                                if FLAGfields.has_key(16):
                                    pos2 = -(pos2 + len(fields[9]))
                        if FLAGfields.has_key(8) and not FLAGfields.has_key(4):
                            ProperPairs = False
                    if NH1 != NH2:
                        ProperPairs = False
                    if ProperPairs:
                        if NH1 == 0:
                            ReadStatsDict['paired_unaligned'] += 1
                        else:
                            if ReadStatsDict['proper_pairs']['all'].has_key(NH1):
                                pass
                            else:
                                ReadStatsDict['proper_pairs']['all'][NH1] = 0
                            ReadStatsDict['proper_pairs']['all'][NH1] += 1
                        if NH1 == 1:
                            chr = fields[2]
                            if ComplexityChromosomesDict.has_key(chr):
                                pass
                            else:
                                ComplexityChromosomesDict[chr] = {}
                            if ComplexityChromosomesDict[chr].has_key(pos1):
                                pass
                            else:
                                ComplexityChromosomesDict[chr][pos1] = {}
                            if ComplexityChromosomesDict[chr][pos1].has_key(pos2):
                                pass
                            else:
                                ComplexityChromosomesDict[chr][pos1][pos2] = 0
                            ComplexityChromosomesDict[chr][pos1][pos2] += 1
                            if ComplexityChromosomesDict[chr][pos1][pos2] == 1:
                                UP += 1
                                UR += 1
                                M0 += 1
                                M1 += 1
                            elif ComplexityChromosomesDict[chr][pos1][pos2] == 2:
                                UR += 1
                                M1 -= 1
                                M2 += 1
                            elif ComplexityChromosomesDict[chr][pos1][pos2] == 3:
                                UR += 1
                                M2 -= 1
                            else:
                                UR += 1
                    else:
                        if ReadStatsDict['not_proper_pairs']['all'].has_key(max(NH1,NH2)):
                            pass
                        else:
                            ReadStatsDict['not_proper_pairs']['all'][max(NH1,NH2)] = 0
                        ReadStatsDict['not_proper_pairs']['all'][max(NH1,NH2)] += 1
                    for lline in AlignmentDict[ID]:
                        fields = lline.split('\t')
                        FLAGfields = FLAG(int(fields[1]))
                        CIGAR = fields[5]
                        if FLAGfields.has_key(64) and NH1 > 0:
                            if 'N' in CIGAR:
                                if ReadStatsDict['spliced_reads']['all'].has_key(NH1):
                                    pass
                                else:
                                    ReadStatsDict['spliced_reads']['all'][NH1] = 0
                                ReadStatsDict['spliced_reads']['all'][NH1] += 1./NH1
                            else:
                                if ReadStatsDict['unspliced_reads']['all'].has_key(NH1):
                                    pass
                                else:
                                    ReadStatsDict['unspliced_reads']['all'][NH1] = 0
                                ReadStatsDict['unspliced_reads']['all'][NH1] += 1./NH1
                        if FLAGfields.has_key(128) and NH2 > 0:
                            if 'N' in CIGAR:
                                if ReadStatsDict['spliced_reads']['all'].has_key(NH2):
                                    pass
                                else:
                                    ReadStatsDict['spliced_reads']['all'][NH2] = 0
                                ReadStatsDict['spliced_reads']['all'][NH2] += 1./NH2
                            else:
                                if ReadStatsDict['unspliced_reads']['all'].has_key(NH2):
                                    pass
                                else:
                                    ReadStatsDict['unspliced_reads']['all'][NH2] = 0
                                ReadStatsDict['unspliced_reads']['all'][NH2] += 1./NH2
                        RL = len(fields[9])
                        if ReadLengthDict.has_key(RL):
                            pass
                        else:
                            ReadLengthDict[RL] = 0
                        if not FLAGfields.has_key(4):
                            if FLAGfields.has_key(64):
                                if NH1 != 0:
                                    ReadLengthDict[RL] += 1./NH1
                                    if doAddNH:
                                        if doPrint:
                                            print lline + '\tNH:i:' + str(NH1)
                                else:
                                    if doPrint:
                                        print lline
                            elif FLAGfields.has_key(128 ):
                                if NH2 != 0:
                                    ReadLengthDict[RL] += 1./NH2
                                    if doAddNH:
                                        if doPrint:
                                            print lline + '\tNH:i:' + str(NH2)
                                else:
                                    if doPrint:
                                        print lline
                            else:
                                if doPrint:
                                    print lline
                else:
                    for lline in AlignmentDict[ID]:
                        fields = lline.split('\t')
                        FLAGfields = FLAG(int(fields[1]))
                        CIGAR = fields[5]
                        if FLAGfields.has_key(4):
                            if len(AlignmentDict[ID]) > 1:
                                print 'read is unmapped but multiple alignments are detected for the read ID, exiting'
                                print AlignmentDict[ID]
                                sys.exit(1)
                            else:
                                ReadStatsDict['unpaired_unaligned'] += 1
#                            if doPrint:
#                                print lline
                        else:
                            RL = len(fields[9])
                            NH = len(AlignmentDict[ID]) 
                            if NH == 1:
                                chr = fields[2]
                                pos1 = int(fields[3])
                                if FLAGfields.has_key(16):
                                    pos1 = -(pos1 + len(fields[9]))
                                pos2 = pos1
                                if ComplexityChromosomesDict.has_key(chr):
                                    pass
                                else:
                                    ComplexityChromosomesDict[chr] = {}
                                if ComplexityChromosomesDict[chr].has_key(pos1):
                                    pass
                                else:
                                    ComplexityChromosomesDict[chr][pos1] = {}
                                if ComplexityChromosomesDict[chr][pos1].has_key(pos2):
                                    pass
                                else:
                                    ComplexityChromosomesDict[chr][pos1][pos2] = 0
                                ComplexityChromosomesDict[chr][pos1][pos2] += 1
                                if ComplexityChromosomesDict[chr][pos1][pos2] == 1:
                                    UP += 1
                                    UR += 1
                                    M0 += 1
                                    M1 += 1
                                elif ComplexityChromosomesDict[chr][pos1][pos2] == 2:
                                    UR += 1
                                    M1 -= 1
                                    M2 += 1
                                elif ComplexityChromosomesDict[chr][pos1][pos2] == 3:
                                    UR += 1
                                    M2 -= 1
                                else:
                                    UR += 1
                            if ReadLengthDict.has_key(RL):
                                pass
                            else:
                                ReadLengthDict[RL] = 0
                            ReadLengthDict[RL] += 1./NH
                            if ReadStatsDict['unpaired']['all'].has_key(NH):
                                pass
                            else:
                                ReadStatsDict['unpaired']['all'][NH] = 0
                            ReadStatsDict['unpaired']['all'][NH] += 1
                            if 'N' in CIGAR:
                                if ReadStatsDict['spliced_reads']['all'].has_key(NH):
                                    pass
                                else:
                                    ReadStatsDict['spliced_reads']['all'][NH] = 0
                                ReadStatsDict['spliced_reads']['all'][NH] += 1
                            else:
                                if ReadStatsDict['unspliced_reads']['all'].has_key(NH):
                                    pass
                                else:
                                    ReadStatsDict['unspliced_reads']['all'][NH] = 0
                                ReadStatsDict['unspliced_reads']['all'][NH] += 1
                            if doAddNH:
                                if doPrint:
                                    print lline + '\tNH:i:' + str(NH)
                            else:
                                if doPrint:
                                    print lline
            AlignmentDict = {}
        if EOF:
            break

    outfile=open(outfilename, 'w')

    outline='unique pairs, proper:\t' + str(ReadStatsDict['proper_pairs']['all'][1])
    outfile.write(outline+'\n')
    outline='unique pairs, not proper:\t' + str(ReadStatsDict['not_proper_pairs']['all'][1])
    outfile.write(outline+'\n')
    outline='unique unpaired reads:\t' + str(ReadStatsDict['unpaired']['all'][1])
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
    outline='multiread pairs, not proper:\t' + str(M)
    outfile.write(outline+'\n')

    M = 0
    for RM in ReadStatsDict['unpaired']['all'].keys():
        if RM != 1:
            M += ReadStatsDict['unpaired']['all'][RM]
    outline='multiread unpaired:\t' + str(int(M))
    outfile.write(outline+'\n')

    outline='unaligned pairs:\t' + str(ReadStatsDict['paired_unaligned'])
    outfile.write(outline+'\n')
    outline='unaligned unpaired reads:\t' + str(ReadStatsDict['unpaired_unaligned'])
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

if __name__=="__main__":
    main()
    