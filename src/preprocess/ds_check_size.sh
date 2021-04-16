#!/bin/bash
size=$(cat $1 | awk 'FNR == 3 {{print}}' | grep -o '[[:digit:]]*')
size=$(($size/1000000))
if [ $size -lt 10 ]
then
   downsampling_size=$size'M'
else
   downsampling_size=$(($size*2/10))'M'
fi
echo $downsampling_size > $2

