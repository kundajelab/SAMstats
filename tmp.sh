#time python SAMstats.sort.stat.filter.py --sorted_sam_file examples/wgs_bam_NA12878_20k_b37_NA12878.sam
time python SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.sam --chunk_size 1000000
#time python SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam

#kernprof -v -l SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam
