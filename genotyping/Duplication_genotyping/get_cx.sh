#!/bin/sh

chr=$1
win_init=$2
init=$3
bam_ins=$4
end=$5
win_end=$6


((./samtools-1.14/samtools depth -r $chr:$win_init-$init $bam_ins | cut -f3) && (./samtools-1.14/samtools depth -r $chr:$end-$win_end $bam_ins | cut -f3))| sort | awk '{arr[NR]=$1} END {x=int((NR)/2); print arr[x]}'

