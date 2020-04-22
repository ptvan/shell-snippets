#!/bin/bash

# in-line variable substitution
echo ${str/foo/bar}

# in-line regex
if [[ $str =~ [0-9]+\.[0-9]+ ]]; then
    # do something
fi

# process a batch of files, in this case converting PGM > JPEG
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

# simply rename a lot of files using a regex
for f in *.png; do mv -n "$f" "${f/-0}"; done

# rename all files to lowercase
for i in *; do mv "$i" "${i,,}"; done

# run a command on all results of a `find`, quick and dirty:
find ./ -type f -name "*.txt" -exec gedit "{}" \;

# nicer version: print0 gives a NULL after each entry to handle tricky filenames
# xargs batches args before passing down the pipe, a bit more efficient
find ./ -type f -name "*.txt" -print0 | xargs -0 gedit

# by extension, move all files found in a find command to common new location
find ./ -name '*.txt' -exec mv {} /new/path/ \;

# trim a PDF to include only certain pages using qpdf
qpdf original.pdf --pages . 2-18 -- trimmed.pdf

# the last command run
!!

# the last command starting with `foo`
!foo

# the last command's argument, in this case the last command was `vi somefile.txt`
vi !$

# compare 2 directories
diff /dir1 /dir2