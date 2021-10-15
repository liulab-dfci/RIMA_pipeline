#!/bin/bash
##Script to downsample the reads for RSeQC

VALUE=$(echo $2 | sed -r 's/M//g')
FACTOR=$(samtools idxstats $1 | cut -f3 | awk -v COUNT=$((VALUE*1000000)) 'BEGIN {total=0} {total += $1} END {print COUNT/total}')

echo "$VALUE"
echo "$FACTOR"
if [[ $FACTOR > 1 ]]
  then
  echo '[ERROR]: Requested number of reads exceeds total read count in' $1 '-- exiting' && exit 1
fi

samtools view -s $FACTOR -b $1 > $1"_"$2"_downsampling.bam"
