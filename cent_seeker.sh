#!/bin/bash

#centromere_seeker
#version 2.1 -- 161201

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

#cent_seeker_v2.1.sh -- 161201
#implemented v2.1
#	-verbose flag operational
#	-commenting complete

#cent_seeker_v2.0.sh -- 161121
#implemented v2.0
#	-the ability to pass flags into the program
#	-broke the program into parts that can be run separately



#usage function to output instruction/manual/help page
usage()
{
cat << EOF
usage: $0 options

cent_seeker, v2.0 161121 - crc
This script seeks centromeres.


OPTIONS:
   -h      Show this message
   -t      Executable location of TRF software
   -s      Input sequence file (fasta or fastq)
   -d      Input dat file (.dat, from TRF)
   -c      Input csv (from prior cent_seeker run)
   -m      Run mode, must be one of: 
   			'full' - run all
   			'trf' - run TRF and stop
   			'posttrf' - run everything after TRF, supply own .dat file
   			'convert' - convert .dat file to .csv file
   			'plot' - plot .csv results
   			'extract' - extract centromere hits into .fasta

   -l      Centromere length, required if running non-verbose (no -v)
   -v      Verbose
EOF
}

#define some variables for use later
TRF=
SEQUENCE=
DAT=
CSV=
MODE=
CENTLENGTH=
VERBOSE=0

#check flags for input, assign input to appropriate variable
while getopts “:ht:s:d:c:l:vm:” OPTION
do
     case $OPTION in
         h)
             usage
             exit 1
             ;;
         t)
             TRF=$OPTARG
             echo trf is at $TRF 
             ;;
         s)
             SEQUENCE=$OPTARG
             echo input seq is $SEQUENCE
             ;;
         d)
             DAT=$OPTARG
             echo input dat is $DAT
             ;;
         c)
             CSV=$OPTARG
             echo input csv is $CSV
             ;;
         l)
             CENTLENGTH=$OPTARG
             echo estimated length of repeat is $CENTLENGTH
             ;;
         v)
             VERBOSE=1
             ;;
         m)
             MODE=$OPTARG
             echo running mode is $MODE
             ;;
         ?)
             usage
             exit
             ;;
     esac
done

#confirm verbose mode
if [[ $VERBOSE = 0 ]]; then
	echo running in non-verbose mode... it said quietly
else
	echo running in VERBOSE mode!
fi

echo
echo
echo
echo
echo

#set of error messages for mis-used flags, no mode listed, incorrect options for a given mode
if [[ -z $MODE ]]; then
	echo "===ERROR=== - No mode listed"
	echo "Please use one of the prescribed modes listed in the help menu"
	echo
     usage
     exit 1
 elif [[ $MODE != "full" && $MODE != "trf" && $MODE != "posttrf" && $MODE != "convert" && $MODE != "plot" && $MODE != "extract"  ]]; then
     echo "===ERROR=== - Mode not found"
     echo "Please use one of the prescribed modes listed in the help menu below"
     echo
     usage
     exit 1
 elif [[ $MODE = "full" || $MODE = "trf" ]]; then
 	if  [[ -z $TRF ]]; then
 		echo "===ERROR=== need TRF path, -t flag"
 		echo
 		usage
 		exit 1
 	elif [[ -z $SEQUENCE ]]; then
 		echo "===ERROR=== need sequence input, -s flag"
 		echo
 		usage
 		exit 1
 	elif [[ $VERBOSE = 0 && -z $CENTLENGTH && $MODE = "full" ]]; then
 		echo "===ERROR=== to run in silent mode centromere length hypothesis is required, -l flag"
 		echo
 		usage
 		exit 1
 	fi
 elif [[ $MODE = "posttrf" || $MODE = "convert" ]]; then
 	if  [[ -z $DAT ]]; then
 		echo "===ERROR=== need .dat file, -d flag"
 		echo
 		usage
 		exit 1
 	elif [[ $VERBOSE = 0 && -z $CENTLENGTH && $MODE = "posttrf" ]]; then
 		echo "===ERROR=== to run in silent mode centromere length hypohthesis is required, -l flag"
 		echo
 		usage
 		exit 1
 	fi
 elif [[ $MODE = "extract" ]]; then
 	if  [[ -z $CSV ]]; then
 		echo "===ERROR=== need .csv file, -c flag"
 		echo
 		usage
 		exit 1
 	elif [[ -z $CENTLENGTH ]]; then
 		echo "===ERROR=== centromere length hypothesis is required, -l flag"
 		echo
 		usage
 		exit 1
 	fi
 elif [[ $MODE = "plot" ]]; then
 	if  [[ -z $CSV ]]; then
 		echo "===ERROR=== need .csv file, -c flag"
 		echo
 		usage
 		exit 1
 	fi
 fi
 
 #error checking has passed, and only satisfactory variable combinations should go past here
 if [[ $VERBOSE = 1 ]]; then
 echo "Proceeding to cent_seeker script"
 fi
###RUN TRF AND OUTPUT DAT FILE
 if [[ $MODE = "full" || $MODE = "trf" ]]; then
 	#check character to determine fasta/fastq
	charcheck=`head -n3 $SEQUENCE | tail -n1 | cut -c1`

	#if statement writes fastq to fasta for trf
	if [  $charcheck != ">"  ]; then 
		if [[ $char != A && $char != G && $char != C && $char != T ]]; then
			if [[ $VERBOSE = 1 ]]; then
				echo "input is fastq"
				echo "converting to fasta"
			fi
			cat $SEQUENCE | grep "@" -A1 | sed 's,@,>,g' | grep [">"ACTG] > input.fasta
			SEQUENCE="input.fasta"
			else
				if [[ $VERBOSE = 1 ]]; then
				echo "input is multi-line fasta"
				fi
		fi
	else
		cat $SEQUENCE > input.fasta
	fi

	#cat a few lines of the input file so the user can see
	if [[ $VERBOSE = 1 ]]; then
		head input.fasta | cut -c1-80
		echo "above is 80 characters from the first 10 lines of sequence from your input file"
		sleep 5 
	fi

	#remove the current/previous trf output
	rm input.fasta.2.5.7.80.10.50.2000.dat

	#run TRF verbose or not
	if [[ $VERBOSE = 1 ]]; then
		$TRF input.fasta 2 5 7 80 10 50 2000 -h
	else
		$TRF input.fasta 2 5 7 80 10 50 2000 -h >/dev/null 2>/dev/null
	fi

	#Run tandem repeat finder
	#$TRF input.fasta 2 5 7 80 10 50 2000 -h >/dev/null

	if [[ $VERBOSE = 1 ]]; then
		echo 
		echo "=========================TANDEM=REPEAT=FINDER=COMPLETE=================="
		echo
	fi

	DAT=`ls input.fasta.2.5.7.80.10.50.2000.dat`
	rm input.fasta
 fi

###CONVERT DAT FILE TO CSV 
 if [[ $MODE = "convert" || $MODE = "full" || $MODE = "posttrf" ]]; then

 	if [[ $VERBOSE = 1 ]]; then
 		echo "converting .dat file into searchable .csv"
 	fi

	#count the trf hits
	grep -n Sequence $DAT | cut -d: -f1 > contig_line.count

	rm trf_all_hits.csv

	m=2

	#number of trf hits to loop through
	totContigs=`wc -l contig_line.count | sed 's,contig_line.count,,g' | sed 's, ,,g'`
	fivePer=$(( totContigs / 20 ))
	count=1
	perCount=0

	#go through the trf hits and write them to csv, 
	#need to add a progress bar here, 5%, 10%... 100%
	
	for start in `cat contig_line.count`
	do
		endplus=`cat contig_line.count | awk 'NR=='$m''`
		end=$((endplus - 1))
		contig=`cat $DAT | awk 'NR=='$start'' | cut -d: -f2 | sed 's, ,,g'`
		
		cat $DAT | awk 'NR=='$start',NR=='$end'' | grep [ACTG] | 		\
			grep -v Sequence | cut -d" " -f1,2,3,4,5,14,15 | sed 's/ /,/g'  \
			| head > trf_hits.tmp
		
		for h in `cat trf_hits.tmp`
		do 
			echo $contig,$h >> trf_all_hits.csv
		done

		#counter for awking through trf output
		m=$((m+1)) 

		if [[ $VERBOSE = 1 ]]; then
			#this if loop prints a progress message at 5, 10, 15% etc
			if [[ $fivePer != 0 ]]; then
				if [[ `expr $count % $fivePer` = 0 ]]; then 
					echo $((5 * perCount))% done
					perCount=$((perCount+1))
				fi
			fi
		fi

		#counter for the overall progress loop
		count=$((count+1))

	done
	if [[ $VERBOSE = 1 ]]; then
		echo "100% done"
	fi
 fi

###PLOT A GRAPH IN R
 if [[ $MODE = "plot" || $MODE = "full" || $MODE = "posttrf"  ]]; then
 
 	#Use R to plot the results, and print to file
	R --slave -f ./snitch.R


	if [[ $VERBOSE = 1 ]]; then
	#ask user to assess the R plot and estimate centromere length
	echo
	echo "=======================PLOTTING=======COMPLETE=========================="
	echo
	fi
	
	if [[ $VERBOSE = 1 ]]; then

		#if no centromere length given in flag, ask in prompt
		if [[ -z $CENTLENGTH ]]; then

			echo Inspect the results, cent_seek_graph.pdf, "for" a pattern of repeats
			echo the occuring "in" multiples of decreasing height from left to right
			echo
			echo If you see a pattern, what is the x-axis value of the shortest and tallest
			echo peak"?"
			echo This may be the length of your centromere, enter it here to proceed:

			read CENTLENGTH
		fi

		echo So your centromere is $CENTLENGTH bases"?"
		echo
		echo I will gather repeats of that length from the trf output"!"

		sleep 3
	fi
 fi

###GET THE LENGTH MATCHED READS INTO A FASTA

 if [[ $MODE = "extract" || $MODE = "full" || $MODE = "posttrf"  ]]; then
 	seq $(( CENTLENGTH - 1 )) $(( CENTLENGTH + 1 )) > centLengths.list
	seq $(( 2 * CENTLENGTH - 2 )) $(( 2 * CENTLENGTH + 2 )) >> centLengths.list
	seq $(( 4 * CENTLENGTH - 4 )) $(( 4 * CENTLENGTH + 4 )) >> centLengths.list

	#for each length collect the hits into a single csv
	totalHits=0

	for l in `cat centLengths.list`
	do
		hits=`cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep -c ,$l,[ACTG][ACTG]`
		if [[ $VERBOSE = 1 ]]; then
			echo $l has $hits matches
		fi
		totalHits=$((totalHits + hits))
		cat trf_all_hits.csv | cut -d, -f1,2,3,6,7,8 | sort -k4 -t, -n | grep ,$l,[ACTG][ACTG] >> cent_matches.csv
	done

	#finally convert that csv to a fasta file
	cat cent_matches.csv | cut -d, -f1,2,3,4,5 | sed 's/,/_/g' | sed 's/^/>/g' | \
		sed 's,_\([GCAT]\),\n\1,g' > centromere_pattern_length_matches.fasta
	cat cent_matches.csv | cut -d, -f1,2,3,4,6 | sed 's/,/_/g' | sed 's/^/>/g' | \
		sed 's,_\([GCAT]\),\n\1,g' > centromere_full_length_matches.fasta

 fi

 echo cent_seeker is DONE, $totalHits potential centromeric sequences found!
