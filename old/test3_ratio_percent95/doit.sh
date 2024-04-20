#!/bin/bash
while read -r line; do
    ./boot.exe $line | tail -1
done < ratio.csv
