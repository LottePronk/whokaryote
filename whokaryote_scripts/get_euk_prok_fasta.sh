#!/bin/bash

INFILE=$1
OUTFILE=$2
EUKHEAD=$3
PROKHEAD=$4
EUKFAST=$5
PROKFAST=$6

awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$INFILE" > "$OUTFILE"

for header in $(cat "$EUKHEAD");
do
	grep -A1 "$header" "$OUTFILE" >> "$EUKFAST"
done


for header in $(cat "$PROKHEAD");
do
	grep -A1 "$header" "$OUTFILE" >> "$PROKFAST"
done