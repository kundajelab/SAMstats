# SAMstats

Scripts that implement samtools flagstat functionality, but provide statistics for individual reads rather than individual alignments

## small example file

### Script can be run in multi-threaded mode (actually turns out slower due to thread lock) or in non-threaded mode:


#### multithreaded version: 
```
python SAMstats.sort.stat.filter.parallelized.thread.py \
  --sorted_sam_file examples/wgs_bam_NA12878_20k_b37_NA12878.sam \
  --outf stats.parallel.txt \
  --chunk_size 1000  \
  --threads 2
```
#### singlethreaded version: 
```
python SAMstats.sort.stat.filter.py \
  --sorted_sam_file examples/wgs_bam_NA12878_20k_b37_NA12878.sam \
  --outf stats.txt \
  --chunk_size 1000 
```

# Benchmarks 
https://docs.google.com/presentation/d/1QyyDMHlYEJdyl-1rH6Z_YoqwT41zsB9MdH-LXjMM8pw/edit#slide=id.g5c152d9ebd_0_40
