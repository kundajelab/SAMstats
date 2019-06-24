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
qcPassed=`samtools view -F 0x200 $sorted_bam | cut -f1 | wc -l` 
echo "qcPassed:$qcPassed"
qcFailed=`samtools view -f 0x200 $sorted_bam | cut -f1 | wc -l` 
echo "qcFailed:$qcFailed"

#secondary, 0x100 bit set 
secondary_qcPassed=`samtools view -F 0x200 -f 0x100 $sorted_bam | cut -f1 | wc -l `
echo "secondary_qcPassed:$secondary_qcPassed"
secondary_qcFailed=`samtools view -f 0x200 -f 0x100 $sorted_bam | cut -f1 | wc -l `
echo "secondary_qcFailed:$secondary_qcFailed"

#supplementary, 0x800 bit set 
supplementary_qcPassed=`samtools view -F 0x200 -f 0x800 $sorted_bam | cut -f1 | wc -l `
echo "supplementary_qcPassed:$supplementary_qcPassed"
supplementary_qcFailed=`samtools view -f 0x200 -f 0x800 $sorted_bam | cut -f1 | wc -l `
echo "supplementary_qcFailed:$supplementary_qcFailed"

#for fraction calculations
#primary, 0x800 bit NOT set 
primary_qcPassed=`samtools view -F 0x200 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "primary_qcPassed:$primary_qcPassed"
primary_qcFailed=`samtools view -f 0x200 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "primary_qcFailed:$primary_qcFailed"


#duplicates,0x400 bit set
duplicates_qcPassed=`samtools view -F 0x200 -f 0x400  $sorted_bam | cut -f1 | wc -l `
echo "duplicates_qcPassed:$duplicates_qcPassed"
duplicates_qcFailed=`samtools view -f 0x200 -f 0x400  $sorted_bam | cut -f1 | wc -l `
echo "duplicates_qcFailed:$duplicates_qcFailed"

#mapped, 0x4 bit not set
mapped_qcPassed=`samtools view -F 0x200 -F 0x4 $sorted_bam | cut -f1 | wc -l `
echo "mapped_qcPassed:$mapped_qcPassed"
mapped_qcFailed=`samtools view -f 0x200 -F 0x4 $sorted_bam | cut -f1 | wc -l `
echo "mapped_qcFailed:$mapped_qcFailed"

#paired in sequencing, 0x1 bit set
pairedInSequencing_qcPassed=`samtools view -F 0x200 -F 0x800 -f 0x1 $sorted_bam | cut -f1 | wc -l `
echo "pairedInSequencing_qcPassed:$pairedInSequencing_qcPassed"
pairedInSequencing_qcFailed=`samtools view -f 0x200 -F 0x800 -f 0x1 $sorted_bam | cut -f1 | wc -l `
echo "pairedInSequencing_qcFailed:$pairedInSequencing_qcFailed"

#read1, both 0x1 and 0x40 bits set
read1_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x40 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "read1_qcPassed:$read1_qcPassed"
read1_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x40 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "read1_qcFailed:$read1_qcFailed"

#read2, both 0x1 and 0x80 bits set
read2_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x80 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "read2_qcPassed:$read2_qcPassed"
read2_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x80 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "read2_qcFailed:$read2_qcFailed"

#properly paired, both 0x1 and 0x2 bits set and 0x4 bit not set
properlyPaired_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x2 -F 0x4 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "properlyPaired_qcPassed:$properlyPaired_qcPassed"
properlyPaired_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x2 -F 0x4 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "properlyPaired_qcFailed:$properlyPaired_qcFailed"

#with itself and mate mapped, 0x1 bit set and neither 0x4 nor 0x8 bits set
withItselfAndMateMapped_qcPassed=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "withItselfAndMateMapped_qcPassed:$withItselfAndMateMapped_qcPassed"
withItselfAndMateMapped_qcFailed=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "withItselfAndMateMapped_qcFailed:$withItselfAndMateMapped_qcFailed"

#singletons, both 0x1 and 0x8 bits set and bit 0x4 not set
singletons_qcPassed=`samtools view -F 0x200 -f 0x1 -f 0x8 -F 0x4 -F 0x800 $sorted_bam | cut -f1 | wc -l `
echo "singletons_qcPassed:$singletons_qcPassed"
singletons_qcFailed=`samtools view -f 0x200 -f 0x1 -f 0x8 -F 0x4 -F 0x800  $sorted_bam | cut -f1 | wc -l `
echo "singletons_qcFailed:$singletons_qcFailed"

#And finally, two rows are given that additionally filter on the reference name (RNAME), mate reference name (MRNM), and mapping quality (MAPQ) fields:
#with mate mapped to a different chr, 0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME
withMateMappedToDiffChrom_qcPassed=`samtools view -F 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 $sorted_bam |awk '{ if ($7 != "=") { print } }'| cut -f1 | wc -l `
echo "withMateMappedToDiffChrom_qcPassed:$withMateMappedToDiffChrom_qcPassed"
withMateMappedToDiffChrom_qcFailed=`samtools view -f 0x200 -f 0x1 -F 0x4 -F 0x8 -F 0x800 $sorted_bam |awk '{ if ($7 != "=") { print } }'| cut -f1 | wc -l `
echo "withMateMappedToDiffChrom_qcFailed:$withMateMappedToDiffChrom_qcFailed"

#with mate mapped to a different chr (mapQ>=5), 0x1 bit set and neither 0x4 nor 0x8 bits set and MRNM not equal to RNAME and MAPQ >= 5
withMateMappedToDiffChromQC5_qcPassed=`samtools view -F 0x200 -F 0x4 -F 0x800 -q5 $sorted_bam | awk '{ if ($7 != "=") { print } }'|cut -f1 | wc -l `
echo "withMateMappedToDiffChromQC5_qcPassed:$withMateMappedToDiffChromQC5_qcPassed"
withMateMappedToDiffChromQC5_qcFailed=`samtools view -f 0x200 -F 0x4 -F 0x800 -q5 $sorted_bam | awk '{ if ($7 != "=") { print } }'|cut -f1 | wc -l `
echo "withMateMappedToDiffChromQC5_qcFailed:$withMateMappedToDiffChromQC5_qcFailed"

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


