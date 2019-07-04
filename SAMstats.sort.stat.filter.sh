#!/bin/bash
usage() {
    echo "Usage:"
    echo "SAMstats.sort.stat.filter.sh InputBamFile OutputStatsFileName NthreadsForSortingBam"
    exit 0
}

[ -z $1 ] && { usage; }

#input file name (bam file) 
inputf=$1

#output file name to save flagstats 
outf=$2 

#threads for sorting the bam file (if it is not currently sorted) 
threads=$3 


#sort the BAM file by read name, if it's not sorted already 
is_sorted=`samtools view -H $inputf | grep "SO:queryname"`
if [ "$is_sorted" == "" ]
then 
    echo "Sorting by read name!" 
    samtools sort -n  --threads $threads -O BAM -o $inputf.sorted $inputf
    sorted_bam=$inputf.sorted
else
    echo "Input file is already sorted by read name"
    sorted_bam=$inputf
fi

#The first row of output gives the total number of reads that are QC pass and fail (according to flag bit 0x200). For example:
#read2
qcPassed_read2=`samtools view -F 0x200 -f 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l` 
echo "qcPassed_read2:$qcPassed_read2"
qcFailed_read2=`samtools view -f 0x200 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l` 
echo "qcFailed_read2:$qcFailed_read2"
qcPassed_read1=`samtools view -F 0x200 -F 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l` 
echo "qcPassed_read1:$qcPassed_read1"
qcFailed_read1=`samtools view -f 0x200 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l` 
echo "qcFailed_read1:$qcFailed_read1"
qcPassed=$(( $qcPassed_read2 + $qcPassed_read1 ))
echo $qcPassed
qcFailed=$(( $qcFailed_read2 + $qcFailed_read1 ))
echo $qcFailed 

#secondary, 0x100 bit set 
secondary_qcPassed_read2=`samtools view -F 0x200 -f 0x100 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "secondary_qcPassed_read2:$secondary_qcPassed_read2"
secondary_qcFailed_read2=`samtools view -f 0x200 -f 0x100 -f 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "secondary_qcFailed_read2:$secondary_qcFailed_read2"
secondary_qcPassed_read1=`samtools view -F 0x200 -f 0x100 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "secondary_qcPassed_read1:$secondary_qcPassed_read1"
secondary_qcFailed_read1=`samtools view -f 0x200 -f 0x100 -F 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "secondary_qcFailed_read1:$secondary_qcFailed_read1"
secondary_qcPassed=$(( $secondary_qcPassed_read2 + $secondary_qcPassed_read1))
echo $secondary_qcPassed
secondary_qcFailed=$(( $secondary_qcFailed_read2 + $secondary_qcFailed_read1))
echo $secondary_qcFailed 


#supplementary, 0x800 bit set 
supplementary_qcPassed_read2=`samtools view -F 0x200 -f 0x800 -f 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "supplementary_qcPassed_read2:$supplementary_qcPassed_read2"
supplementary_qcFailed_read2=`samtools view -f 0x200 -f 0x800  -f 0x80 $sorted_bam | cut -f1| sort |uniq | wc -l `
echo "supplementary_qcFailed_read2:$supplementary_qcFailed_read2"
supplementary_qcPassed_read1=`samtools view -F 0x200 -f 0x800 -F 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "supplementary_qcPassed_read1:$supplementary_qcPassed_read1"
supplementary_qcFailed_read1=`samtools view -f 0x200 -f 0x800 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "supplementary_qcFailed_read1:$supplementary_qcFailed_read1"
supplementary_qcPassed=$(( $supplementary_qcPassed_read2 + $supplementary_qcPassed_read1))
echo $supplementary_qcPassed 
supplementary_qcFailed=$(( $supplementary_qcFailed_read2 + $supplementary_qcFailed_read1))
echo $supplementary_qcFailed 

#for fraction calculations
#primary, 0x800 bit NOT set 
primary_qcPassed_read2=`samtools view -F 0x200 -F 0x800 -f 0x80  $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "primary_qcPassed_read2:$primary_qcPassed_read2"
primary_qcFailed_read2=`samtools view -f 0x200 -F 0x800 -f 0x80  $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "primary_qcFailed_read2:$primary_qcFailed_read2"
primary_qcPassed_read1=`samtools view -F 0x200 -F 0x800 -F 0x80  $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "primary_qcPassed_read1:$primary_qcPassed_read1"
primary_qcFailed_read1=`samtools view -f 0x200 -F 0x800 -F 0x80  $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "primary_qcFailed_read1:$primary_qcFailed_read1"
primary_qcPassed=$(( $primary_qcPassed_read2 + $primary_qcPassed_read1))
echo $primary_qcPassed 
primary_qcFailed=$(( $primary_qcFailed_read2 + $primary_qcFailed_read1))
echo $primary_qcFailed

#duplicates,0x400 bit set
duplicates_qcPassed_read2=`samtools view -F 0x200 -f 0x400 -f 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "duplicates_qcPassed_read2:$duplicates_qcPassed_read2"
duplicates_qcFailed_read2=`samtools view -f 0x200 -f 0x400 -f 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "duplicates_qcFailed_read2:$duplicates_qcFailed_read2"
duplicates_qcPassed_read1=`samtools view -F 0x200 -f 0x400 -F 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "duplicates_qcPassed_read1:$duplicates_qcPassed_read1"
duplicates_qcFailed_read1=`samtools view -f 0x200 -f 0x400 -F 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "duplicates_qcFailed_read1:$duplicates_qcFailed_read1"
duplicates_qcPassed=$(( $duplicates_qcPassed_read2 + $duplicates_qcPassed_read1))
echo $duplicates_qcPassed
duplicates_qcFailed=$(( $duplicates_qcFailed_read2 + $duplicates_qcFailed_read1))
echo $duplicates_qcFailed

#mapped, 0x4 bit not set
mapped_qcPassed_read2=`samtools view -F 0x200 -F 0x4 -f 0x80 $sorted_bam | cut -f1 | sort | uniq |  wc -l `
echo "mapped_qcPassed_read2:$mapped_qcPassed_read2"
mapped_qcFailed_read2=`samtools view -f 0x200 -F 0x4 -f 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "mapped_qcFailed_read2:$mapped_qcFailed_read2"
mapped_qcPassed_read1=`samtools view -F 0x200 -F 0x4 -F 0x80 $sorted_bam | cut -f1 | sort | uniq |  wc -l `
echo "mapped_qcPassed_read1:$mapped_qcPassed_read1"
mapped_qcFailed_read1=`samtools view -f 0x200 -F 0x4 -F 0x80 $sorted_bam | cut -f1 | sort | uniq | wc -l `
echo "mapped_qcFailed_read1:$mapped_qcFailed_read1"
mapped_qcPassed=$(( $mapped_qcPassed_read2 + $mapped_qcPassed_read1))
echo $mapped_qcPassed
mapped_qcFailed=$(( $mapped_qcFailed_read2 + $mapped_qcFailed_read1))
echo $mapped_qcFailed

#paired in sequencing, 0x1 bit set
pairedInSequencing_qcPassed_read2=`samtools view -F 0x200 -F 0x800 -f 0x1 -f 0x80 $sorted_bam | cut -f1 |sort | uniq | wc -l `
echo "pairedInSequencing_qcPassed_read2:$pairedInSequencing_qcPassed_read2"
pairedInSequencing_qcFailed_read2=`samtools view -f 0x200 -F 0x800 -f 0x1 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "pairedInSequencing_qcFailed_read2:$pairedInSequencing_qcFailed_read2"
pairedInSequencing_qcPassed_read1=`samtools view -F 0x200 -F 0x800 -f 0x1 -F 0x80 $sorted_bam | cut -f1 |sort | uniq | wc -l `
echo "pairedInSequencing_qcPassed_read1:$pairedInSequencing_qcPassed_read1"
pairedInSequencing_qcFailed_read1=`samtools view -f 0x200 -F 0x800 -f 0x1 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "pairedInSequencing_qcFailed_read1:$pairedInSequencing_qcFailed_read1"
pairedInSequencing_qcPassed=$(( $pairedInSequencing_qcPassed_read2 + $pairedInSequencing_qcPassed_read1))
echo $pairedInSequencing_qcPassed
pairedInSequencing_qcFailed=$(( $pairedInSequencing_qcFailed_read2 + $pairedInSequencing_qcFailed_read1))
echo $pairedInSequencing_qcFailed

#read1, both 0x1 and 0x40 bits set
read1_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x40 -F 0x800 $sorted_bam | cut -f1 |sort| uniq | wc -l `
echo "read1_qcPassed:$read1_qcPassed"
read1_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x40 -F 0x800 $sorted_bam | cut -f1 | sort|uniq | wc -l `
echo "read1_qcFailed:$read1_qcFailed"
#read2, both 0x1 and 0x80 bits set
read2_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x80 -F 0x800 $sorted_bam | cut -f1 | sort|uniq | wc -l `
echo "read2_qcPassed:$read2_qcPassed"
read2_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x80 -F 0x800 $sorted_bam | cut -f1 | sort|uniq | wc -l `
echo "read2_qcFailed:$read2_qcFailed"

#properly paired, both 0x1 and 0x2 bits set and 0x4 bit not set
properlyPaired_qcPassed_read2=`samtools view -F 0x200 -f 0x1 -f 0x2 -F 0x4 -F 0x800 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "properlyPaired_qcPassed_read2:$properlyPaired_qcPassed_read2"
properlyPaired_qcFailed_read2=`samtools view -f 0x200 -f 0x1 -f 0x2 -F 0x4 -F 0x800 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "properlyPaired_qcFailed_read2:$properlyPaired_qcFailed_read2"
properlyPaired_qcPassed_read1=`samtools view -F 0x200 -f 0x1 -f 0x2 -F 0x4 -F 0x800 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "properlyPaired_qcPassed_read1:$properlyPaired_qcPassed_read1"
properlyPaired_qcFailed_read1=`samtools view -f 0x200 -f 0x1 -f 0x2 -F 0x4 -F 0x800 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "properlyPaired_qcFailed_read1:$properlyPaired_qcFailed_read1"
properlyPaired_qcPassed=$(( $properlyPaired_qcPassed_read2 + $properlyPaired_qcPassed_read1))
echo $properlyPaired_qcPassed
properlyPaired_qcFailed=$(( $properlyPaired_qcFailed_read2 + $properlyPaired_qcFailed_read1))
echo $properlyPaired_qcFailed


#with itself and mate mapped, 0x1 bit set and neither 0x4 nor 0x8 bits set
withItselfAndMateMapped_qcPassed_read2=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "withItselfAndMateMapped_qcPassed_read2:$withItselfAndMateMapped_qcPassed_read2"
withItselfAndMateMapped_qcFailed_read2=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "withItselfAndMateMapped_qcFailed_read2:$withItselfAndMateMapped_qcFailed_read2"
withItselfAndMateMapped_qcPassed_read1=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "withItselfAndMateMapped_qcPassed_read1:$withItselfAndMateMapped_qcPassed_read1"
withItselfAndMateMapped_qcFailed_read1=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "withItselfAndMateMapped_qcFailed_read1:$withItselfAndMateMapped_qcFailed_read1"
withItselfAndMateMapped_qcPassed=$(( $withItselfAndMateMapped_qcPassed_read2 + $withItselfAndMateMapped_qcPassed_read1))
echo $withItselfAndMateMapped_qcPassed
withItselfAndMateMapped_qcFailed=$(( $withItselfAndMateMapped_qcFailed_read2 + $withItselfAndMateMapped_qcFailed_read1))
echo $withItselfAndMateMapped_qcFailed 

#singletons, both 0x1 and 0x8 bits set and bit 0x4 not set
singletons_qcPassed_read2=`samtools view -F 0x200 -f 0x1 -f 0x8 -F 0x4 -F 0x800 -f 0x80 $sorted_bam | cut -f1 | sort |uniq  | wc -l `
echo "singletons_qcPassed_read2:$singletons_qcPassed_read2"
singletons_qcFailed_read2=`samtools view -f 0x200 -f 0x1 -f 0x8 -F 0x4 -F 0x800 -f 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "singletons_qcFailed_read2:$singletons_qcFailed_read2"
singletons_qcPassed_read1=`samtools view -F 0x200 -f 0x1 -f 0x8 -F 0x4 -F 0x800 -F 0x80 $sorted_bam | cut -f1 | sort |uniq  | wc -l `
echo "singletons_qcPassed_read1:$singletons_qcPassed_read1"
singletons_qcFailed_read1=`samtools view -f 0x200 -f 0x1 -f 0x8 -F 0x4 -F 0x800 -F 0x80 $sorted_bam | cut -f1 | sort |uniq | wc -l `
echo "singletons_qcFailed_read1:$singletons_qcFailed_read1"
singletons_qcPassed=$(( $singletons_qcPassed_read2 + $singletons_qcPassed_read1))
echo $singletons_qcPassed
singletons_qcFailed=$(( $singletons_qcFailed_read2 + $singletons_qcFailed_read1))
echo $singletons_qcFailed 

#And finally, two rows are given that additionally filter on the reference name (RNAME), mate reference name (MRNM), and mapping quality (MAPQ) fields:
#with mate mapped to a different chr, 0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME
withMateMappedToDiffChrom_qcPassed_read2=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -f 0x80 $sorted_bam |awk '{ if ($7 != "=") { print } }'| cut -f1 | sort |uniq |  wc -l `
echo "withMateMappedToDiffChrom_qcPassed_read2:$withMateMappedToDiffChrom_qcPassed_read2"
withMateMappedToDiffChrom_qcFailed_read2=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -f 0x80 $sorted_bam |awk '{ if ($7 != "=") { print } }'| cut -f1 | sort |uniq  | wc -l `
echo "withMateMappedToDiffChrom_qcFailed_read2:$withMateMappedToDiffChrom_qcFailed_read2"
withMateMappedToDiffChrom_qcPassed_read1=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -f 0x80 $sorted_bam |awk '{ if ($7 != "=") { print } }'| cut -f1 | sort |uniq |  wc -l `
echo "withMateMappedToDiffChrom_qcPassed_read1:$withMateMappedToDiffChrom_qcPassed_read1"
withMateMappedToDiffChrom_qcFailed_read1=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 -f 0x80 $sorted_bam |awk '{ if ($7 != "=") { print } }'| cut -f1 | sort |uniq  | wc -l `
echo "withMateMappedToDiffChrom_qcFailed_read1:$withMateMappedToDiffChrom_qcFailed_read1"
withMateMappedToDiffChrom_qcPassed=$(( $withMateMappedToDiffChrom_qcPassed_read2 + $withMateMappedToDiffChrom_qcPassed_read1))
echo $withMateMappedToDiffChrom_qcPassed
withMateMappedToDiffChrom_qcFailed=$(( $withMateMappedToDiffChrom_qcFailed_read2 + $withMateMappedToDiffChrom_qcFailed_read1))
echo $withMateMappedToDiffChrom_qcFailed

#with mate mapped to a different chr (mapQ>=5), 0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME and MAPQ >= 5
withMateMappedToDiffChromQC5_qcPassed_read2=`samtools view -F 0x200 -F 0x4 -F 0x800 -f 0x80 -q5 $sorted_bam | awk '{ if ($7 != "=") { print } }'|cut -f1 | sort |uniq  |wc -l `
echo "withMateMappedToDiffChromQC5_qcPassed_read2:$withMateMappedToDiffChromQC5_qcPassed_read2"
withMateMappedToDiffChromQC5_qcFailed_read2=`samtools view -f 0x200 -F 0x4 -F 0x800 -f 0x80 -q5 $sorted_bam | awk '{ if ($7 != "=") { print } }'|cut -f1 | sort |uniq  | wc -l `
echo "withMateMappedToDiffChromQC5_qcFailed_read2:$withMateMappedToDiffChromQC5_qcFailed_read2"
withMateMappedToDiffChromQC5_qcPassed_read1=`samtools view -F 0x200 -F 0x4 -F 0x800 -F 0x80 -q5 $sorted_bam | awk '{ if ($7 != "=") { print } }'|cut -f1 | sort |uniq  |wc -l `
echo "withMateMappedToDiffChromQC5_qcPassed_read1:$withMateMappedToDiffChromQC5_qcPassed_read1"
withMateMappedToDiffChromQC5_qcFailed_read1=`samtools view -f 0x200 -F 0x4 -F 0x800 -F 0x80 -q5 $sorted_bam | awk '{ if ($7 != "=") { print } }'|cut -f1 | sort |uniq  | wc -l `
echo "withMateMappedToDiffChromQC5_qcFailed_read1:$withMateMappedToDiffChromQC5_qcFailed_read1"
withMateMappedToDiffChromQC5_qcPassed=$(( $withMateMappedToDiffChromQC5_qcPassed_read2 + $withMateMappedToDiffChromQC5_qcPassed_read1))
echo $withMateMappedToDiffChromQC5_qcPassed
withMateMappedToDiffChromQC5_qcFailed=$(( $withMateMappedToDiffChromQC5_qcFailed_read2 + $withMateMappedToDiffChromQC5_qcFailed_read1))
echo $withMateMappedToDiffChromQC5_qcFailed

#Get fractions of mapped, paired, singletons relative to the total (this is provided in the original FlagStat output, so we provide it here)
multiplier=100
if (($qcPassed == 0))
then 
    percent_mapped_qcPassed="NaN"
    percent_properlyPaired_qcPassed="NaN"
    percent_singletons_qcPassed="NaN"
else
    percent_mapped_qcPassed=`echo $multiplier*$mapped_qcPassed/$qcPassed|bc -l`
    percent_properlyPaired_qcPassed=`echo $multiplier*$properlyPaired_qcPassed/$primary_qcPassed|bc -l`
    percent_singletons_qcPassed=`echo $multiplier*$singletons_qcPassed/$primary_qcPassed|bc -l`
fi
if (($qcFailed == 0))
then 
    percent_mapped_qcFailed="NaN"
    percent_properlyPaired_qcFailed="NaN"
    percent_singletons_qcFailed="NaN"
else
    percent_mapped_qcFailed=`echo $multiplier*$mapped_qcFailed/$qcFailed|bc -l`
    percent_properlyPaired_qcFailed=`echo $multiplier*$properlyPaired_qcFailed/$primary_qcFailed|bc -l`
    percent_singletons_qcFailed=`echo $multiplier*$singletons_qcFailed/$primary_qcFailed|bc -l`
fi
#write the output file 
echo "$qcPassed + $qcFailed in total (QC-passed reads + QC-failed reads)" >$outf 
echo "$secondary_qcPassed + $secondary_qcFailed secondary">> $outf 
echo "$supplementary_qcPassed + $supplementary_qcFailed supplementary">> $outf 
echo "$duplicates_qcPassed + $duplicates_qcFailed duplicates">> $outf 
echo "$mapped_qcPassed + $mapped_qcFailed mapped ($percent_mapped_qcPassed % : $percent_mapped_qcFailed % )">> $outf 
echo "$pairedInSequencing_qcPassed + $pairedInSequencing_qcFailed paired in sequencing">> $outf 
echo "$read1_qcPassed + $read1_qcFailed read1">> $outf 
echo "$read2_qcPassed + $read1_qcFailed read2">> $outf 
echo "$properlyPaired_qcPassed + $properlyPaired_qcFailed properly paired ($percent_properlyPaired_qcPassed % : $percent_properlyPaired_qcFailed % )">> $outf 
echo "$withItselfAndMateMapped_qcPassed + $withItselfAndMateMapped_qcFailed with itself and mate mapped">> $outf
echo "$singletons_qcPassed + $singletons_qcFailed singletons ( $percent_singletons_qcPassed % : $percent_singletons_qcFailed % )">> $outf
echo "$withMateMappedToDiffChrom_qcPassed + $withMateMappedToDiffChrom_qcFailed with mate mapped to a different chr">> $outf
echo "$withMateMappedToDiffChromQC5_qcPassed + $withMateMappedToDiffChromQC5_qcFailed with mate mapped to a different chr (mapQ>=5)">> $outf


