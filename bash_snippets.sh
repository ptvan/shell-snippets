#!/bin/bash

# process a batch of files
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

# simply rename a lot of files
for f in *.png; do mv -n "$f" "${f/-0}"; done

# run a command on all results of a `find`, quick and dirty:
find ./ -type f -name "*.txt" -exec gedit "{}" \;

# nicer version: print0 gives a NULL after each entry to handle tricky filenames
# xargs batches args before passing down the pipe, a bit more efficient
find ./ -type f -name "*.txt" -print0 | xargs -0 gedit

