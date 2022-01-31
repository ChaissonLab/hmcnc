#!/usr/bin/env bash


## requires bedtools 
## requires merged repeat annotations in HMM/annotations
## requires hg38 excluded regions bed in HMM/annotations


## usage:
## ./filter_finalize_call.sh <calls.bed>  HMM/annotation/region_EXCLUDE.no_alts.bed  HMM/annotation/repeatMask.merged.bed


call=$1
excl=$2
repeat=$3 

grep PASS $call |awk '$4>2' | cut -f 1-8 | intersectBed -v -a stdin -b $excl |\
intersectBed -a stdin -b $repeat -sorted  -loj| \
awk 'BEGIN{s=0;e=0;} {s=$10;e=$11;if($10==-1){print $0"\t"0;}  else {if(s<$2){s=$2;} if(e>$3){e=$3;} {print $0"\t"e-s;} }}' | \
groupBy -g 1-8 -c 12 -o sum |\
awk '{print $0"\t"$9/$6;}' > $call.filter.bed


grep PASS $call | awk '$4>2' |cut -f 1-8 |mergeBed| intersectBed -v -a stdin -b $excl | \
intersectBed -a stdin -b $repeat -sorted  -loj| \
awk 'BEGIN{s=0;e=0;} {s=$5;e=$6;if($5==-1){print $0"\t"0;}  else {if(s<$2){s=$2;} if(e>$3){e=$3;} {print $0"\t"e-s;} }}' |\
groupBy -g 1-3 -c 7 -o sum| awk '{print $0"\t"($3-$2);}'|awk '{print $0"\t"$4/$5;}' > $call.merge.filter.bed

intersectBed -loj -a  $call.filter.bed -b $call.merge.filter.bed| awk '$16<0.8' > dup.filtered_final.$call

grep PASS $call | awk '$4<2' |cut -f 1-8|intersectBed -v -a stdin -b $excl > del.filtered_final.$call

cat  del.filtered_final.$call dup.filtered_final.$call | sort -k1,1 -k2,2n |cut -f 1-8 > filtered_final.$call

