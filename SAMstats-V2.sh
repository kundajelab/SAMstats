#!/bin/bash

#input file name (bam file) 
inputf=$1

#output file name to start stats 
outf=$2 

#threads for sorting the bam file 
threads=$3 

#sort the BAM file by name 
samtools sort -n $inputf --threads $threads  -o $inputf.sorted

#index the sorted bam file 
samtools index $inputf.sorted 

#The first row of output gives the total number of reads that are QC pass and fail (according to flag bit 0x200). For example:
qcPassed=`samtools view -F 0x200 $inputf.sorted | cut -f1 | uniq | wc -l` 
echo "qcPassed:$qcPassed"
qcFailed=`samtools view -f 0x200 $inputf.sorted | cut -f1 | uniq | wc -l` 
echo "qcFailed:$qcFailed"

#secondary, 0x100 bit set 
secondary_qcPassed=`samtools view -F 0x200 -f 0x100 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "secondary_qcPassed:$secondary_qcPassed"
secondary_qcFailed=`samtools view -f 0x200 -f 0x100 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "secondary_qcFailed:$secondary_qcFailed"

#supplementary, 0x800 bit set 
supplementary_qcPassed=`samtools view -F 0x200 -f 0x800 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "supplementary_qcPassed:$supplementary_qcPassed"
supplementary_qcFailed=`samtools view -f 0x200 -f 0x800 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "supplementary_qcFailed:$supplementary_qcFailed"


#duplicates,0x400 bit set
duplicates_qcPassed=`samtools view -F 0x200 -f 0x400 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "duplicates_qcPassed:$duplicates_qcPassed"
duplicates_qcFailed=`samtools view -f 0x200 -f 0x400 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "duplicates_qcFailed:$duplicates_qcFailed"

#mapped, 0x4 bit not set
mapped_qcPassed=`samtools view -F 0x200 -F 0x4 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "mapped_qcPassed:$mapped_qcPassed"
mapped_qcFailed=`samtools view -f 0x200 -F 0x4 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "mapped_qcFailed:$mapped_qcFailed"

#paired in sequencing, 0x1 bit set
pairedInsequencing_qcPassed=`samtools view -F 0x200 -f 0x1 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "pairedInsequencing_qcPassed:$pairedInsequencing_qcPassed"
pairedInsequencing_qcFailed=`samtools view -f 0x200 -f 0x1 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "pairedInsequencing_qcFailed:$pairedInsequencing_qcFailed"

#read1, both 0x1 and 0x40 bits set
read1_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x40 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "read1_qcPassed:$read1_qcPassed"
read1_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x40 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "read1_qcFailed:$read1_qcFailed"

#read2, both 0x1 and 0x80 bits set
read2_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x80 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "read2_qcPassed:$read2_qcPassed"
read2_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x80 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "read2_qcFailed:$read2_qcFailed"

#properly paired, both 0x1 and 0x2 bits set and 0x4 bit not set
properlyPaired_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x2 -F 0x4 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "properlyPaired_qcPassed:$properlyPaired_qcPassed"
properlyPaired_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x2 -F 0x4 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "properlyPaired_qcFailed:$properlyPaired_qcFailed"

#with itself and mate mapped, 0x1 bit set and neither 0x4 nor 0x8 bits set
withItselfAndMateMapped_qcPassed=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "withItselfAndMateMapped_qcPassed:$withItselfAndMateMapped_qcPassed"
withItselfAndMateMapped_qcFailed=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "withItselfAndMateMapped_qcFailed:$withItselfAndMateMapped_qcFailed"

#singletons, both 0x1 and 0x8 bits set and bit 0x4 not set
singletons_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x8 -F 0x4 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "singletons_qcPassed:$singletons_qcPassed"
singletons_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x8 -F 0x4  $inputf.sorted | cut -f1 | uniq | wc -l `
echo "singletons_qcFailed:$singletons_qcFailed"

#And finally, two rows are given that additionally filter on the reference name (RNAME), mate reference name (MRNM), and mapping quality (MAPQ) fields:
#with mate mapped to a different chr, 0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME
withMateMappedToDiffChrom_qcPassed=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "withMateMappedToDiffChrom_qcPassed:$withMateMappedToDiffChrom_qcPassed"
withMateMappedToDiffChrom_qcFailed=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "withMateMappedToDiffChrom_qcFailed:$withMateMappedToDiffChrom_qcFailed"

#with mate mapped to a different chr (mapQ>=5), 0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME and MAPQ >= 5
withMateMappedToDifferentChromQC5_qcPassed=`samtools view -F 0x200 -F 0x4 -q5 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "withMateMappedToDifferentChromQC5_qcPassed:$withMateMappedToDifferentChromQC5_qcPassed"
withMateMappedToDifferentChromQC5_qcFailed=`samtools view -f 0x200 -F 0x4 -q5 $inputf.sorted | cut -f1 | uniq | wc -l `
echo "withMateMappedToDifferentChromQC5_qcFailed:$withMateMappedToDifferentChromQC5_qcFailed"

#write the output file 


