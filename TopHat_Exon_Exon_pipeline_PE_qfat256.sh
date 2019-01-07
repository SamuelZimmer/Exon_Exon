#!/bin/bash


####TopHat######

while getopts "t:g:1:2:r:n:o:" option
do
	case "${option}"
	in
	o) OUT=${OPTARG};;
	g) GTF=${OPTARG};;
	t) THREADS=${OPTARG};;
	1) FASTQ1=${OPTARG};;
	2) FASTQ2=${OPTARG};;
	r) REF=${OPTARG};;
	s) SAMPLE=${OPTARG};;
	n) NAME=${OPTARG};;
	esac
done

mkdir $OUT
cd $OUT

export JOB="tophat -p $THREADS -G $GTF --no-novel-juncs --keep-tmp --min-intron-length 40 -o $OUT $REF $FASTQ1 $FASTQ2"

echo $JOB > job_time.txt
echo "start time" >> job_time.txt
date >> job_time.txt

ml bioinformatics/tophat/2.1.1 bioinformatics/bowtie2/2.2.8
$JOB



######Intersects#####

cat junctions.bed | sort -n -r -k 5 -o junctions.sorted.bed
sort -k1,1 -k2,2n junctions.sorted.bed > junctions.2.sorted.bed
bedtools intersect -a  accepted_hits.bam -b junctions.2.sorted.bed -F 0.5 -sorted | samtools view - > abam.F.0.5.intersects




####Filtre#####


###Filtre sur mapq

cat abam.F.0.5.intersects | awk '$5 >30' > abam.F.0.5.intersects.Filter.1


###Filtre sur CIGAR ( M >= 20 ) et Filtre sur taille de l'intron ####

cat abam.F.0.5.intersects.Filter.1 | grep -P "MD:Z:([2-9]\d+|\d{3,})\^[A-Z]([2-9]\d+|\d{3,})" | grep -P "[0-9]*[0-9][0-9]N" > abam.F.0.5.intersects.Filter.2





###intersect to bam###


samtools view -H accepted_hits.bam header > header.txt
#sed header.txt -i -e 's/SN:/SN:chr/'
sed header.txt -i -e 's/Mito/M/'
for i in `ls abam.*`; do echo $i; mv $i ${NAME}.$i; done
#cat ${NAME}.abam.F.0.5.intersects.Filter.2 | awk '{ $3="chr"$3; print $0 }' | tr ' ' '\t' > ${NAME}.abam.F.0.5.intersects.chr
#cat header.txt ${NAME}.abam.F.0.5.intersects.chr > ${NAME}.abam.F.0.5.intersects.sam
cat header.txt ${NAME}.abam.F.0.5.intersects.Filter.2 > ${NAME}.abam.F.0.5.intersects.sam
sed ${NAME}.abam.F.0.5.intersects.sam -i -e 's/Mito/M/g'
samtools view -bS ${NAME}.abam.F.0.5.intersects.sam -o ${NAME}.abam.F.0.5.intersects.bam
samtools sort ${NAME}.abam.F.0.5.intersects.bam > ${NAME}.abam.F.0.5.intersects.sorted.bam
samtools index ${NAME}.abam.F.0.5.intersects.sorted.bam


####Filtre sur nombre de reads par region (>20)

ml mugqic/bedtools/2.26.0
bedtools genomecov -ibam ${NAME}.abam.F.0.5.intersects.sorted.bam -bg > ${NAME}.abam.F.0.5.intersects.bedgraph

bedtools merge -d 1 -c 4,4 -o max,sum -i ${NAME}.abam.F.0.5.intersects.bedgraph > ${NAME}.peaks.bedgraph

######



echo "end time" >> job_time.txt
date >> job_time.txt
