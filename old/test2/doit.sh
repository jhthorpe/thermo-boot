#!/bin/bash
for i in 10 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 95 100 105 110 115 120 125 130 135 140; do 
    echo -n "$i : "
    ./boot.exe $i | tail -1 
done
