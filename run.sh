#time python SAMstats.sort.stat.filter.parallelized.thread.py --sorted_sam_file examples/wgs_bam_NA12878_20k_b37_NA12878.sam --threads 20
#time python SAMstats.sort.stat.filter.parallelized.thread.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.sam --chunk_size 1000000 --threads 40
#python SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.sam --chunk_size 1000000 --threads 40  --outf pe.raw.stats
#time python SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam --threads 40 --outf pe.dedup.stats --chunk_size 1000000
#kernprof -v -l SAMstats.sort.stat.filter.py --sorted_sam_file examples/2D_MCF10A_1_CAGAGAGG_R1.merged.nodup.sam
