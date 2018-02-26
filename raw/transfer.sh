#! /bin/bash

#This file should transfer everything new from the MIDAS computer and then unzip it.

rsync -ravpt york@134.158.196.38:/TapeData/201803/* .
#gzip -dk *.gz

for i in `find . -name 'R*.gz'`
do
#echo $i
j=${i%_*}
j=$j"_0"
#echo $j
if [ ! -e $j ]
then
	echo $i
	echo $j
	gzip -d < $i > $j
fi
done

