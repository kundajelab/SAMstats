#small example file 
#time python SAMstats.sort.stat.filter.parallelized.thread.py --sorted_sam_file examples/wgs_bam_NA12878_20k_b37_NA12878.sam
#time python SAMstats.sort.stat.filterx.py --sorted_sam_file examples/wgs_bam_NA12878_20k_b37_NA12878.sam

#100 million reads 
#time python SAMstats.sort.stat.filter.parallelized.thread.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.sam --chunk_size 1000000 --threads 1
#time python SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.sam --chunk_size 1000000

#100 million reads after dup filtering
#time python SAMstats.sort.stat.filter.parallelized.thread.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam --threads 1 --chunk_size 1000000
#time python SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam --chunk_size 1000000 

