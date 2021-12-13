#!/bin/sh

input_file=$1 	# file with coordinateS
bam_ins=$2		# Bam file
sample=$3
#OUTPUT
#chr init end len diff win gen median_cx

while read chr init end len diff win gen; do 

	win_len=`echo "$len * 2" | bc`; 
	win_init=`echo "$init - $win_len" | bc`; 
	win_end=`echo "$end + $win_len" | bc`; 

	if [ $win_len -lt $diff ]; then
		
		cx=`./get_cx.sh $chr $win_init $init $bam_ins $end $win_end`;

        echo "$chr\t$init\t$end\t$len\t$diff\t$win\t$gen\t$cx" >> outputs/$sample"_"$chr"_readinfo.txt";

	elif [ $diff == 0 ]; then
		
		cx=`./get_cx.sh $chr $win_init $init $bam_ins $end $win_end`;
        echo "$chr\t$init\t$end\t$len\t$diff\t$win\t$gen\t$cx" > outputs/$sample"_"$chr"_readinfo.txt";

	else
		
	echo "$chr\t$init\t$end\t$len\t$diff\t$win\t$gen\t$cx\tCLOSE" >> outputs/$sample"_"$chr"_readinfo.txt";
	fi



done < $input_file 

