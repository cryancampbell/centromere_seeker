#!/bin/bash

#centromere_seeker

#takes two inputs:
#1) the location (or alias) of your local copy of tandem repeat finder (trf)
#2) the sequence file to analyze (fasta or fastq)

#this script looks for long repeats in pacbio data with trf
#and plots them visually to isolate potential centromeric
#repeats

#must have both tandem repeat finder (trf, https://tandem.bu.edu/trf/trf.html)
#and R (https://www.r-project.org/) installed

#centromere seeker has 2 stages and 2 outputs
#
#stage 1 - determines underlying repeat patterns in the sequencing data
#			this stage runs trf, then converts that data to csv and plots
#			the resulting data in R. the R plot output is "cent_seek_graph.pdf"
#			on screen output will inform you when this has been created, and
#			once it has, you should assess it for patterns that are consistent
#			with a broadly present centromere-sized repeat (see 
#			example_cent_seek_graph.pdf). if you observe a pattern input the
#			centromere size that pattern indicates into the terminal prompt to
#			enter stage 2
#
#stage 2 - the program now switches and gathers the sequence data with pattern
#			lengths that match your hypotesized centromere size (n, 2n, and 4n)
#			there are two outputs:
#
#			centromere_pattern_length_matches.fasta	- includes the n-length 
#														pattern that trf found
#
#			centromere_full_length_matches.fasta	- includes the full match
#														that trf found, which 
#														may include slight 
#														variations in the 
#														pattern
#
#with these fasta outputs, use your favorite alignment program (e.g. AliView)
#to determine if there is a large and common connection across these patterns. 
#if so, that repeat may be very common in the genome of interest which may be 
#indicative of a centromeric repeat.

#INPUTS -- change to flags
#location of trf
TRF=$1
#seq file (fasta or fastq)
seqinput=$2

#check character to determine fasta/fastq
charcheck=`head -n3 $seqinput | tail -n1 | cut -c1`

#if statement writes fastq to fasta
#if [ $charcheck = "+" ]; then
if [  $charcheck != ">"  ]; then 
	if [[ $char != A && $char != G && $char != C && $char != T ]]; then
	echo input is fastq
	echo converting to fasta
	cat $seqinput | grep "@" -A1 | sed 's,@,>,g' | grep [">"ACTG] > input.fasta
	seqinput="input.fasta"
	else
	echo input is multi-line fasta
	fi
else
	cat $seqinput > input.fasta
fi

#output the input file so the user can see
head input.fasta | cut -c1-80

echo above is the fasta sequence from your input "file"
sleep 5 

rm *dat

#Run tandem repeat finder

$TRF input.fasta 2 5 7 80 10 50 2000 -h

echo "=========================TANDEM=REPEAT=FINDER=COMPLETE=================="

##flag to start here only (e.g. -d lists "dat" file, and if "dat" file exists start with it)

trfoutput=`ls *dat`

#count the trf hits
grep -n Sequence $trfoutput | cut -d: -f1 > contig_line.count

rm trf_all_hits.csv

m=2

#number of trf hits to loop through
totContigs=`wc -l contig_line.count | cut -d" " -f1`

#show this number to the user, so they know how many to expect
for z in `seq 1 $totContigs`
do
	echo -n "_"
	sleep .01
done
echo

#go through the trf hits and write them to csv, echoing a "." to show progress
for start in `cat contig_line.count`
do
	endplus=`cat contig_line.count | awk 'NR=='$m''`
	end=$((endplus - 1))
	contig=`cat $trfoutput | awk 'NR=='$start'' | cut -d: -f2 | sed 's, ,,g'`
	#echo $contig,$start,$end
	echo -n "."
	#cat $trfoutput | awk 'NR=='$start',NR=='$end'' | grep [ACTG] | grep -v Sequence | cut -d" " -f1,2,3,4,5,14,15 | sed 's/ /,/g' | head > trf_hits.tmp
	cat $trfoutput | awk 'NR=='$start',NR=='$end'' | grep [ACTG] | 		\
		grep -v Sequence | cut -d" " -f1,2,3,4,5,14,15 | sed 's/ /,/g'  \
		| head > trf_hits.tmp
	for h in `cat trf_hits.tmp`
	do 
		echo $contig,$h >> trf_all_hits.csv
	done
	m=$((m+1)) 
done
echo

#get rid of the converted fasta
rm input.fasta

#Use R to plot the results, and print to file
R --slave -f ../centromere_seeker/snitch.R

#ask user to assess the R plot and estimate centromere length
echo "=======================PLOTTING=======COMPLETE=========================="
echo Inspect the results, cent_seek_graph.pdf, "for" a pattern of repeats
echo the occuring "in" multiples of decreasing height from left to right
echo
echo If you see a pattern, what is the x-axis value of the shortest and tallest
echo peak"?"
echo This may be the length of your centromere, enter it here to proceed:

read centLength

echo So your centromere is $centLength "?"
echo
echo I will gather repeats of that length from the trf output"!"

sleep 3

#from user input, gather all the hits from trf into one csv

#nOneBEG=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((centLength - 1)), | head -n1 | cut -d: -f1`
nOneBEG=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | \ 
	grep -n ,$((centLength - 1)), | head -n1 | cut -d: -f1`
#echo $nOneBEG
#nOneEND=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((centLength + 1)), | tail -n1 | cut -d: -f1`
nOneEND=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | \	
	grep -n ,$((centLength + 1)), | tail -n1 | cut -d: -f1`
#echo $nOneEND

if [ -z "$nOneBEG" ]; then 
	if [ -z "$nOneEND"]; then
		echo No Hits at 1x
	else
		#cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nOneEND'' >> cent_matches.csv
		cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | \
			awk 'NR=='$nOneEND'' >> cent_matches.csv
	fi
elif [ -z "$nOneEND" ]; then
	#cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nOneBEG'' >> cent_matches.csv
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | 	\
		awk 'NR=='$nOneBEG'' >> cent_matches.csv
else
	#cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nOneBEG',NR=='$nOneEND'' >> cent_matches.csv
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | 	\
	awk 'NR=='$nOneBEG',NR=='$nOneEND'' >> cent_matches.csv
fi

nTwoBEG=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((2 * centLength - 2)), | head -n1 | cut -d: -f1`
#echo $nTwoBEG
nTwoEND=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((2 * centLength + 2)), | tail -n1 | cut -d: -f1`
#echo $nTwoEND

if [ -z "$nTwoBEG" ]; then 
	if [ -z "$nTwoEND"]; then
		echo No Hits at 2x
	else
		cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nTwoEND'' >> cent_matches.csv
	fi
elif [ -z "$nTwoEND" ]; then
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nTwoBEG'' >> cent_matches.csv
else
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nTwoBEG',NR=='$nTwoEND'' >> cent_matches.csv
fi


nFourBEG=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((4 * centLength - 4)), | head -n1 | cut -d: -f1`
#echo $nFourBEG
nFourEND=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((4 * centLength + 4)), | tail -n1 | cut -d: -f1`
#echo $nFourEND

if [ -z "$nFourBEG" ]; then 
	if [ -z "$nFourEND"]; then
		echo No Hits at 4x
	else
		cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nFourEND'' >> cent_matches.csv
	fi
elif [ -z "$nFourEND" ]; then
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nFourBEG'' >> cent_matches.csv
else
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nFourBEG',NR=='$nFourEND'' >> cent_matches.csv
fi

#finally, convert the trf centromere matches to two fasta files for aligning
cat cent_matches.csv | cut -d, -f1,2,3,4,5 | sed 's/,/_/g' | sed 's/^/>/g' | \
	sed 's,_\([GCAT]\),\n\1,g' > centromere_pattern_length_matches.fasta
cat cent_matches.csv | cut -d, -f1,2,3,4,6 | sed 's/,/_/g' | sed 's/^/>/g' | \
	sed 's,_\([GCAT]\),\n\1,g' > centromere_full_length_matches.fasta


