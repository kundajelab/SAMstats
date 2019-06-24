import sys 
import argparse 
import pandas as pd 

def parse_args(): 
    parser=argparse.ArgumentParser(description="Compute SAM file mapping statistics for a SAM file sorted by read name")
    parser.add_argument("--inputSortedSamFile",help="SAM file that has been sorted by read name")
    parser.add_argument("--outputStatsFileName") 
    return parser.parse_args()

def main(): 
    argse=parse_args() 
    bam=args.inputSortedBamFile 
    outf=args.outputStatsFileName
    if bam=="-":
        bam=sys.stdin
    flagstat=dict() 
    cur_flag_stat=dict() 
    cur_ID=None
    for line in pd.read_csv(bam, header=None, chunksize=BS,sep='\t'):
        
if __name__=="__main__": 
    main() 

