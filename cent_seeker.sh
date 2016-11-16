#!/bin/bash

#centromere_seeker

#this script looks for long repeats in pacbio data with trf
#and plots them visually to isolate potential centromeric
#repeats

#must have both tandem repeat finder (trf, https://tandem.bu.edu/trf/trf.html)
#and R (https://www.r-project.org/) installed

#INPUTS
#location of trf
#seq file (fasta or fastq)
TRF=$1
seqinput=$2

##add if loop here to turn fastq into fasta (grep out the ">" rows)
charcheck=`head -n3 $seqinput | tail -n1 | cut -c1`

if [ $charcheck = "+" ]; then 
	echo input is fastq
	echo converting to fasta
	cat $seqinput | grep "@" -A1 | sed 's,@,>,g' | grep [">"ACTG] > input.fasta
	seqinput="input.fasta"
else
	cat $seqinput > input.fasta
fi

head input.fasta | cut -c1-80

rm *dat

#STEP 1 RUN TRF
#
$TRF input.fasta 2 5 7 80 10 50 2000 -h
#
#trf409.linux64 $seqinput 2 5 7 80 10 50 2000 -l 6

echo "===========================TANDEM=REPEAT=FINDER=COMPLETE========================================="

trfoutput=`ls *dat`

grep -n Sequence $trfoutput | cut -d: -f1 > contig_line.count

rm trf_all_hits.csv

m=2

#CREATE contig_line.count = the line count of each "contig" in the trf output

totContigs=`wc -l contig_line.count | cut -d" " -f1`

for z in `seq 1 $totContigs`
do
	echo -n "_"
	sleep .01
done
echo

for start in `cat contig_line.count`
do
	endplus=`cat contig_line.count | awk 'NR=='$m''`
	end=$((endplus - 1))
	contig=`cat $trfoutput | awk 'NR=='$start'' | cut -d: -f2 | sed 's, ,,g'`
	#echo $contig,$start,$end
	echo -n "."
	cat $trfoutput | awk 'NR=='$start',NR=='$end'' | grep [ACTG] | grep -v Sequence | cut -d" " -f1,2,3,4,5,14,15 | sed 's/ /,/g' | head > trf_hits.tmp
	for h in `cat trf_hits.tmp`
	do 
		echo $contig,$h >> trf_all_hits.csv
	done
	m=$((m+1)) 
done
echo

rm input.fasta

#add some R stuff, output a graph, take input perhaps? a guess as to the length of repeat?

R --slave -f ../centromere_seeker/snitch.R

echo "===========================PLOTTING=======COMPLETE========================================="
echo Inspect the results, cent_seek_graph.pdf, "for" a pattern of repeats
echo the pattern will occur "in" multiples of decreasing height from left to right
echo
echo If you see a pattern, how long is the shortest "("x-axis")" and tallest "("y-axis")" peak"?"
echo This may be the length of your centromere, please enter it here to proceed:
read centLength

echo So your centromere is $centLength "?"
echo I will gather repeats of that length from the trf output"!"

nOneBEG=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((centLength - 1)), | head -n1 | cut -d: -f1`
#echo $nOneBEG
nOneEND=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -n ,$((centLength + 1)), | tail -n1 | cut -d: -f1`
#echo $nOneEND

if [ -z "$nOneBEG" ]; then 
	if [ -z "$nOneEND"]; then
		echo No Hits at 1x
	else
		cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nOneEND'' >> cent_matches.csv
	fi
elif [ -z "$nOneEND" ]; then
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nOneBEG'' >> cent_matches.csv
else
	cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | awk 'NR=='$nOneBEG',NR=='$nOneEND'' >> cent_matches.csv
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


cat cent_matches.csv | cut -d, -f1,2,3,4,5 | sed 's/,/_/g' | sed 's/^/>/g' | sed 's,_\([GCAT]\),\n\1,g' > centromere_pattern_length_matches.fasta
cat cent_matches.csv | cut -d, -f1,2,3,4,6 | sed 's/,/_/g' | sed 's/^/>/g' | sed 's,_\([GCAT]\),\n\1,g' > centromere_full_length_matches.fasta


