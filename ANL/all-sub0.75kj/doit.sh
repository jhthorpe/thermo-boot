#!/bin/bash
while read -r line; do
#    echo $line
    ./boot.exe $line | tail -1
#    echo "====================================="
done < ratio.csv
