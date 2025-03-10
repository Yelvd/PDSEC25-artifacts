#!/bin/bash
d=./results/das6-freq-analysis/1122467/

cd $d


for f in likwid_perfctr*; do
    filename=$(basename -- "$f")
    filename="${filename%.*}" 
    A=$(echo $filename | awk -F_ '{print $4}' )

    if [[ $(( $A % 24 )) == 0 ]]; then
        echo Keep $f
    else
        rm $f
    fi

done
