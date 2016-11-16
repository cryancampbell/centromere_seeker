#!/bin/bash

#centromere_seeker

#this script looks for long repeats in pacbio data with trf
#and plots them visually to isolate potential centromeric
#repeats

#must have both tandem repeat finder (trf, https://tandem.bu.edu/trf/trf.html)
#and R (https://www.r-project.org/) installed

#INPUTS
#seq file from PacBio
fastqinput=$1

#rm *dat

#STEP 1 RUN TRF
#trf409.linux64 $fastqinput 2 5 7 80 10 50 2000 -l 6

trfoutput=`ls *dat`

grep -n Contig $trfoutput | cut -d: -f1 > contig_line.count

rm trf_all_hits.csv

m=2

#CREATE contig_line.count = the line count of each "contig" in the trf output

for start in `cat contig_line.count`
do
	endplus=`cat contig_line.count | awk 'NR=='$m''`
	end=$((endplus - 1))
	contig=`cat $trfoutput | awk 'NR=='$start'' | cut -d: -f2 | sed 's, ,,g'`
	echo $contig,$start,$end
	cat $trfoutput | awk 'NR=='$start',NR=='$end'' | grep [ACTG] | grep -v Sequence | cut -d" " -f1,2,3,4,5,14,15 | sed 's/ /,/g' | head > trf_hits.tmp
	for h in `cat trf_hits.tmp`
	do 
		echo $contig,$h >> trf_all_hits.csv
	done
	m=$((m+1)) 
done

#add some R stuff, output a graph, take input perhaps? a guess as to the length of repeat?




#grep ,5[0123456],.*,5[0123456], trf_all_hits.csv | cut -d, -f1,4,5,6,7 | sed 's,Contig,>Contig,g' | sed 's/,5[0123456],\([ACTG]\)/\n\1/g' > as53_potential_patterns50-56.fasta