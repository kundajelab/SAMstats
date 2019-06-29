## Computes samtools flagstat mapping statistics for an input SAM file sorted by read name (i.e. SO= queryname). Mapping statistics are computed per *READ* rather than per alignment 
## (i.e. reads with multiple alignments are counted once)
## stat --> filter (use --filter_before_stat set to False) 
## filter --> stat (use --filter_before_stat set to True) 

import sys 
import argparse 
import pdb 
import line_profiler

def parse_args(): 
    parser=argparse.ArgumentParser(description="Compute SAM file mapping statistics for a SAM file sorted by read name")
    parser.add_argument("--sorted_sam_file",help="Input SAM file. Use '-' if input is being piped from stdin. File must be sorted by read name.")
    parser.add_argument("--outf",default=None,help="Output file name to store alignment statistics. The statistics will be printed to stdout if no file is provided") 
    parser.add_argument("--chunk_size",type=int,default=100000,help="Number of lines to read a time from sortedSamFile")
    parser.add_argument("--threads",type=int,default=1,help="number of threads to use")
    return parser.parse_args()


def write_output_file(outf,
                      stats,
                      global_flagstat,
                      percent_mapped,
                      percent_properly_paired,
                      percent_singletons):
    if outf!=None: 
        outf=open(outf,'w') 
    num_fields=len(stats) 
    for i in range(len(stats)): 
        outstring=' '.join([str(global_flagstat[0][i]),"+",str(global_flagstat[1][i]),stats[i]])
        percent_string=""
        if stats[i]=="mapped": 
            percent_string="".join(["(",
                            str(percent_mapped[0]),
                            ":",
                            str(percent_mapped[1]),
                            ")"])
        elif stats[i]=="properly paired": 
            percent_string="".join(["(",
                            str(percent_properly_paired[0]),
                            ":",
                            str(percent_properly_paired[1]),
                            ")"])
        elif stats[i]=="singletons": 
            percent_string="".join(["(",
                            str(percent_singletons[0]),
                            ":",
                            str(percent_singletons[1]),
                            ")"])
        outstring=' '.join([outstring,percent_string])
        if outf==None: 
            print(outstring) 
        else: 
            outf.write(outstring+'\n') 
    
def calculate_percent(field_index,global_flagstat):
    '''
    field_index = index (0 - 12) in the global_flagstat array indicating the number of reads for the field of interest. Refer to http://www.htslib.org/doc/samtools.html for field index ordering 
    returns % of reads meeting the criteria for the specified field among the reads that pass qc and among the reads that fail qc 
    '''
    qc_passed_primary=global_flagstat[0][13] 
    qc_failed_primary=global_flagstat[1][13]
    qc_passed_field=global_flagstat[0][field_index] 
    qc_failed_field=global_flagstat[1][field_index] 
    field_percent_qc_passed="NA" 
    if qc_passed_primary >0: 
        field_percent_qc_passed=round(qc_passed_field/qc_passed_primary,3) 
    field_percent_qc_failed="NA" 
    if qc_failed_primary > 0: 
        field_percent_qc_failed=round(qc_failed_field/qc_failed_primary,3) 
    return field_percent_qc_passed, field_percent_qc_failed 


#@profile
def add_read_stats(flag,mapq,rnext,cur_flagstat): 
    '''
    implements flagstat logic from http://www.htslib.org/doc/samtools.html to calculate each of the 13 flags
    also keep track of number of primary reads for calculating fraction of mapped reads in the global summary stats 
    '''
    temp_flagstat=[0]*14
    qc_passed = flag & 0x200 == 0    
    secondary =  0x100 & flag == 0x100
    supplementary = flag & 0x800 == 0x800
    primary=flag & 0x800 == 0 
    duplicates=flag & 0x400 == 0x400
    mapped=flag & 0x4==0 

    #note: it is not stated in the samtools documentation, but the paired read statistics are *only* computed for primary reads, so we add a check that primary==1 for all remaining stats 
    paired_in_sequencing= (flag & 0x1 == 0x1) & (primary==1) 
    read1=(flag & 0x1 == 0x1) & (flag & 0x40 == 0x40) & (primary==1)
    read2=(flag & 0x1 == 0x1) & (flag & 0x80 == 0x80) & (primary==1)
    properly_paired=(flag & 0x1==0x1) & (flag & 0x2==0x2) & (flag & 0x4==0) & (primary==1) 
    with_itself_and_mate_mapped=(flag & 0x1==0x1) & (flag & 0x4==0) & (flag & 0x8==0) & (primary==1) 
    singleton=(flag & 0x1==0x1) & (flag & 0x8==0x8) & (flag & 0x4==0) & (primary==1) 
    with_mate_mapped_to_different_chrom=(flag & 0x1==0x1) & (flag & 0x4==0) & (flag & 0x8 ==0) & (primary==1) & (rnext !="=") 
    with_mate_mapped_to_different_chrom_q5=with_mate_mapped_to_different_chrom & (mapq>=5)
    
    #populate the temp_flagstat array with flag values
    temp_flagstat[0]=qc_passed
    temp_flagstat[1]=secondary 
    temp_flagstat[2]=supplementary 
    temp_flagstat[3]=duplicates 
    temp_flagstat[4]=mapped 
    temp_flagstat[5]=paired_in_sequencing 
    temp_flagstat[6]=read1
    temp_flagstat[7]=read2
    temp_flagstat[8]=properly_paired
    temp_flagstat[9]=with_itself_and_mate_mapped
    temp_flagstat[10]=singleton 
    temp_flagstat[11]=with_mate_mapped_to_different_chrom
    temp_flagstat[12]=with_mate_mapped_to_different_chrom_q5
    temp_flagstat[13]=primary 
    
    #compute bitwise or of cur_flagstat and temp_flagstat to see if a flag is set for any of the alignments for the given read 
    #the first entry in cur_flagstat is the numpy array with stats for reads that pass qc; the second entry is the numpy array with stats for reads that fail qc 
    if qc_passed==True:
        cur_flagstat[0]=max(cur_flagstat[0],temp_flagstat)
    else:
        cur_flagstat[1]=max(cur_flagstat[1],temp_flagstat)
    return cur_flagstat 
#@profile
def initialize_flagstat(stat): 
    '''
    numpy array of length 14 (for the 13 metrics in flagstat + one entry to keep track of primary reads) keeps track of number of reads that meet the criteria for each of the 13 statistics 
    numpy array initialized for qc_passed and qc_failed read subsets 
    returns [qc_passed, qc_failed] list ofnumpy arrays 
    initial read counts for each stat are set to 0 
    '''
    qc_passed=[0]*14
    qc_failed=[0]*14
    flagstat=[qc_passed,qc_failed] 
    return flagstat 
#@profile
def update_flagstat_for_readname(global_flagstat, cur_flagstat): 
    '''
    global_flagstat -- flagstat statistics for the full dataset 
    cur_flagstat -- flagstat statistics for a fixed QNAME+SEQ
    This function  updates the global_flagstat arrays with the count of flag stats for a given QNAME +SEQ 
    '''
    global_flagstat[0]=[i+j for i,j in zip(global_flagstat[0],cur_flagstat[0])]
    global_flagstat[1]=[i+j for i,j in zip(global_flagstat[1],cur_flagstat[1])]
    return global_flagstat

#@profile
def main(): 
    #read in the arguments 
    args=parse_args() 
    outf=args.outf
    if args.sorted_sam_file=="-":
        sam=sys.stdin
    else:
        sam=open(args.sorted_sam_file,'r') 
    #13 stats that flagstats reports
    # the order of the stat in the global_flagstat and cur_flagstat arrays matches the order of the stats list 
    stats=['total',
           'secondary',
           'supplementary',
           'duplicates',
           'mapped',
           'paired in sequencing',
           'read1',
           'read2',
           'properly paired',
           'with itself and mate mapped',
           'singletons',
           'with mate mapped to a different chr',
           'with mate mapped to a different chr q5']

    global_flagstat=initialize_flagstat(stats) 

    #since sorted_sam_file is sorted with SO=queryname, we keep track of statistics for a given read name and merge with the larger dict when no further reads 
    #with that name are encountered 
    cur_flagstat=initialize_flagstat(stats) 
    cur_seqid=None
    
    #we read the input SAM file in chunks of size chunksize, iterating one line at a time 
    #skip comment lines in the header that start with @ 
    # use columns : 
    #    0 = QNAME, 
    #    1 = FLAG,
    #    4 = MAPQ, 
    #    6 = RNEXT (reference name of the mate/next read) 
    #    9 = SEQ 
    # ignore all other columns 
    print("starting flag calculation...") 
    line_number=0
    while True:
        cur_lines=sam.readlines(args.chunk_size)
        if len(cur_lines)==0:
            break 
        for line in cur_lines: 
            if line_number % args.chunk_size==0:
                print(line_number)
            line_number+=1
            if line.startswith('@'):
                #this is a comment, we skip
                continue
            tokens=line.split('\t')
            flag=int(tokens[1]) 
            mapq=int(tokens[4])
            rnext=tokens[6] 
            seq=tokens[9] 
            new_seqid=tokens[0]+seq
            #read1=str(flag & 0x40 == 0x40)
            #new_seqid=tokens[0]+read1

            if new_seqid!=cur_seqid: 
                #we are finished processing the readname "cur_ID", updated the global flag statistics for full dataset with statistics for this read 
                global_flagstat=update_flagstat_for_readname(global_flagstat,cur_flagstat) 
                #reinitialize read-specific flagstat for current read name 
                cur_flagstat=initialize_flagstat(stats) 

            #identify flagstat values for current read 
            cur_flagstat=add_read_stats(flag,mapq,rnext,cur_flagstat)
            cur_seqid=new_seqid
        
        
    #we have parsed all the reads in the file
    global_flagstat=update_flagstat_for_readname(global_flagstat,cur_flagstat) 
    #calculate % mapped, properly paired, singletons from total primary reads 
    #note: 4, 8, 10 are the numpy array indices in global_flagstat where the counts for these fields are stored 
    print("finished parsing lines, summarizing...") 
    percent_mapped=calculate_percent(4,global_flagstat)
    percent_properly_paired=calculate_percent(8,global_flagstat) 
    percent_singletons=calculate_percent(10,global_flagstat) 
    
    #write the output file 
    write_output_file(args.outf,
                      stats,
                      global_flagstat,
                      percent_mapped,
                      percent_properly_paired,
                      percent_singletons)

if __name__=="__main__": 
    main() 
