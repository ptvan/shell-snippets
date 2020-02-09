#!/bin/bash
# this script works on the CMU Faces image collection
# (https://archive.ics.uci.edu/ml/datasets/CMU+Face+Images)

if [ -d "faces" ]
then
    echo "Found 'faces' subdir, converting all PGM files found within to PNG..."
    for f in `find faces/ -name *.pgm`
    do 
        regex=${f%.*}
        newfile="$regex.jpg"
        echo "converting $f to $newfile"
        convert $f $newfile
    done

else
    echo "ERROR : 'faces' subdir not found, exiting...."
fi