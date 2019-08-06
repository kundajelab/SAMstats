#small example file 
time SAMstats --sorted_sam_file wgs_bam_NA12878_20k_b37_NA12878.sam
time SAMstatsParallel --sorted_sam_file wgs_bam_NA12878_20k_b37_NA12878.sam

#100 million reads 
#time SAMstatsParallel --sorted_sam_file 2D_MCF10A_1_CAGAGAGG_R1.merged.sam --chunk_size 1000000 --threads 1
#time SAMstats --sorted_sam_file 2D_MCF10A_1_CAGAGAGG_R1.merged.sam --chunk_size 1000000

#100 million reads after dup filtering
#time SAMstatsParallel --sorted_sam_file 2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam --threads 1 --chunk_size 1000000
#time SAMstats --sorted_sam_file 2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam --chunk_size 1000000 

