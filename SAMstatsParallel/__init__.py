## Computes samtools flagstat mapping statistics for an input SAM file sorted by read name (i.e. SO= queryname). Mapping statistics are computed per *READ* rather than per alignment 
## (i.e. reads with multiple alignments are counted once)

import sys 
import argparse 
import multiprocess as mp
import threading
#from multiprocessing.pool import ThreadPool

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


def parse_args(): 
    parser=argparse.ArgumentParser(description="Compute SAM file mapping statistics for a SAM file sorted by read name")
    parser.add_argument("--sorted_sam_file",help="Input SAM file. Use '-' if input is being piped from stdin. File must be sorted by read name.",required=True)
    parser.add_argument("--outf",default=None,help="Output file name to store alignment statistics. The statistics will be printed to stdout if no file is provided") 
    parser.add_argument("--chunk_size",type=int,default=100000,help="Number of lines to read a time from sortedSamFile")
    parser.add_argument("--threads",type=int,default=1,help="number of threads to use. Note: the default is ")
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
        field_percent_qc_passed=round(float(qc_passed_field)/qc_passed_primary,3) 
    field_percent_qc_failed="NA" 
    if qc_failed_primary > 0: 
        field_percent_qc_failed=round(float(qc_failed_field)/qc_failed_primary,3) 
    return field_percent_qc_passed, field_percent_qc_failed 

def add_read_stats(flag,mapq,rnext,cur_flagstat): 
    '''
    implements flagstat logic from http://www.htslib.org/doc/samtools.html to calculate each of the 13 flags
    also keep track of number of primary reads for calculating fraction of mapped reads in the global summary stats 
    '''
    temp_flagstat=[]
    qc_passed=flag & 0x200 == 0
    temp_flagstat.append(qc_passed) #qc_passed
    temp_flagstat.append(flag & 0x100 == 0x100) #secondary 
    temp_flagstat.append(flag & 0x800 == 0x800) #supplementary 
    temp_flagstat.append(flag & 0x400 == 0x400) #duplicates
    temp_flagstat.append(flag & 0x4 == 0) #mapped

    #note: it is not stated in the samtools documentation, but the paired read statistics are *only* computed for primary reads, so we add a check that primary==1 for all remaining stats
    primary= flag & 0x800 ==0  #primary
    paired_in_sequencing= (flag & 0x1 == 0x1) & (primary==1)
    temp_flagstat.append(paired_in_sequencing) #paired_in_sequencing 
    temp_flagstat.append(paired_in_sequencing & (flag & 0x40 == 0x40)) #read1
    temp_flagstat.append(paired_in_sequencing & (flag & 0x80 == 0x80)) #read2
    temp_flagstat.append(paired_in_sequencing & (flag & 0x2==0x2) & (flag & 0x4==0)) #properly aligned 
    with_itself_and_mate_mapped=(paired_in_sequencing & (flag & 0x4==0) & (flag & 0x8==0)) #with_itself_and_mate_mapped
    temp_flagstat.append(with_itself_and_mate_mapped)
    temp_flagstat.append(paired_in_sequencing & (flag & 0x8==0x8) & (flag & 0x4==0)) #singleton
    with_mate_mapped_to_different_chrom=(with_itself_and_mate_mapped & (rnext !="="))
    temp_flagstat.append(with_mate_mapped_to_different_chrom) #with_mate_mapped_to_different_chrom
    temp_flagstat.append(with_mate_mapped_to_different_chrom & (mapq>=5)) #with_mate_mapped_to_different_chrom_q5
    temp_flagstat.append(primary)
    #compute bitwise or of cur_flagstat and temp_flagstat to see if a flag is set for any of the alignments for the given read 
    #the first entry in cur_flagstat is the numpy array with stats for reads that pass qc; the second entry is the numpy array with stats for reads that fail qc 
    if qc_passed==True:
        cur_flagstat[0]=[i|j for i,j in zip(cur_flagstat[0],temp_flagstat)]
    else:
        cur_flagstat[1]=[i|j for i,j in zip(cur_flagstat[1],temp_flagstat)]
    return cur_flagstat


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


def update_flagstat_for_readname(global_flagstat, cur_flagstat): 
    '''
    global_flagstat -- flagstat statistics for the full dataset 
    cur_flagstat -- flagstat statistics for a fixed QNAME+SEQ
    This function  updates the global_flagstat arrays with the count of flag stats for a given QNAME +SEQ 
    '''
    global_flagstat[0]=[i+j for i,j in zip(global_flagstat[0],cur_flagstat[0])]
    global_flagstat[1]=[i+j for i,j in zip(global_flagstat[1],cur_flagstat[1])]
    return global_flagstat


def read_group_stats(input_q,output_q,i):
    '''
    pulls a read group from input_q with same seq_id, calculates flagstat for reads in the group, 
    flagstat statistics are added to output_q
    '''
    print(i)
    while True:
        seqs=input_q.get()
        if seqs is None:
            return 
        else:            
            cur_flagstat=initialize_flagstat(stats)
            for seq in seqs:
                flag=seq[0]
                mapq=seq[1]
                rnext=seq[2]
                cur_flagstat=add_read_stats(flag,mapq,rnext,cur_flagstat)
            output_q.put(cur_flagstat)


def aggregate_read_groups(output_q,outf):
    '''
    aggregates flag stats for unique reads to a single list 
    writes output file 
    '''
    global_flagstat=initialize_flagstat(stats)
    while True:
        cur_flagstat=output_q.get()
        if cur_flagstat is None:                
            #calculate % mapped, properly paired, singletons from total primary reads 
            #note: 4, 8, 10 are the numpy array indices in global_flagstat where the counts for these fields are stored 
            print("finished parsing lines, summarizing...") 
            percent_mapped=calculate_percent(4,global_flagstat)
            percent_properly_paired=calculate_percent(8,global_flagstat) 
            percent_singletons=calculate_percent(10,global_flagstat) 
    
            #write the output file 
            write_output_file(outf,
                              stats,
                              global_flagstat,
                              percent_mapped,
                              percent_properly_paired,
                              percent_singletons)
            return 
        else:
            #we are finished processing the readname "cur_ID", updated the global flag statistics for full dataset with statistics for this read 
            global_flagstat=update_flagstat_for_readname(global_flagstat,cur_flagstat) 


def main(): 
    #read in the arguments 
    args=parse_args() 
    outf=args.outf
    if args.sorted_sam_file=="-":
        sam=sys.stdin
    else:
        sam=open(args.sorted_sam_file,'r') 
    
    #create the input/output queues for multi-threading
    input_q=mp.Queue()
    output_q=mp.Queue() 

    #initialize the worker threads that calculate flagstat statistics for each read group (i.e. alignments with the same readname + sequence)
    workers=[]
    for i in range(args.threads):
        worker=threading.Thread(target=read_group_stats,args=(input_q, output_q,i))
        worker.start() 
        workers.append(worker)
        
    #start the aggregator -- sums flagstat statistics for unique read groups 
    aggregator=threading.Thread(target=aggregate_read_groups,args=(output_q,args.outf))
    aggregator.start() 
    
    #since sorted_sam_file is sorted with SO=queryname, we keep track of statistics for a given read name and merge with the larger dict when no further reads 
    #with that name are encountered
    cur_seq_id=None
    cur_seq_tokens=dict() 
    
    #we read the input SAM file in chunks of size chunksize, iterating one line at a time 
    #skip comment lines in the header that start with @ 
    # use columns : 
    #    0 = QNAME, 
    #    1 = FLAG,
    #    4 = MAPQ, 
    #    6 = RNEXT (reference name of the mate/next read) 
    # ignore all other columns 
    print("starting flag calculation...") 
    line_number=0
    while True:
        #read in a user-specified batch size of reads at a time 
        cur_lines=sam.readlines(args.chunk_size)
        if len(cur_lines)==0:
            #we're at the end of the file 
            break        
        for line in cur_lines:
            #update the log indicating how many lines have been processed 
            if line_number % args.chunk_size==0:
                print(line_number)
            line_number+=1

            #this is a comment, we skip
            if line.startswith('@'):
                continue
            
            #get the fields of interest from the current read
            tokens=line.split()
            flag=int(tokens[1]) 
            mapq=int(tokens[4])
            rnext=tokens[6]
            new_seq_id=tokens[0] 
            read2=str(flag & 0x80==0x80) #check to see if this is the 2nd read in a pair, append 0/1 to the QNAME to
            #generate unique name for read1 & read2 
            new_seq_id_expanded=tokens[0]+read2    
            if new_seq_id!=cur_seq_id:
                #We have seen all the reads with a given name, we now process the flags for that read group, separating by read pair (each read in pair is key) 
                for key in cur_seq_tokens: 
                    input_q.put(cur_seq_tokens[key])
                #reinitialize token list for the new seq_id 
                cur_seq_tokens=dict() 
            if new_seq_id_expanded not in cur_seq_tokens:
                cur_seq_tokens[new_seq_id_expanded]=[[flag,mapq,rnext]]
            else:
                cur_seq_tokens[new_seq_id_expanded].append([flag,mapq,rnext])
            #update cur_seq_id to reflect the new_seq_id we have observed 
            cur_seq_id=new_seq_id
            
    #don't forget to process the final read group!
    for key in cur_seq_tokens:
        input_q.put(cur_seq_tokens[key])
    
    # When the worker sees None in the queue, it knows to terminate
    for i in range(args.threads): 
        input_q.put(None)
        
    #we are finished, join the worker and aggregator processes
    for worker in workers: 
        worker.join()
    output_q.put(None) 
    aggregator.join()


if __name__=="__main__": 
    main() 
