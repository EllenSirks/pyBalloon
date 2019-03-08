#!/bin/bash -l

file="/home/ellen/Desktop/SuperBIT/Flight_data/flights.txt"

while IFS= read -r line || [[ -n "$line" ]]
do
    echo $line
    python2 /home/ellen/Desktop/SuperBIT/pyBalloon-master/retrospective.py $line True &>/dev/null
done < $file

