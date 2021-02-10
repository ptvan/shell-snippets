#!/bin/bash

##### FILE HANDLING

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

# compare 2 directories
diff /dir1 /dir2

#### REGEX
# in-line variable substitution
echo ${str/foo/bar}

# in-line regex
if [[ $str =~ [0-9]+\.[0-9]+ ]]; then
    # do something
fi

##### CONVERSIONS
# split paired CUE and FLAC file into individual FLAC files
shnsplit -f file.cue -t %n-%t -o flac file.flac

# convert FLAC to MP3 using parallelized ffmpeg
parallel ffmpeg -i {} -qscale:a 0 {.}.mp3 ::: ./*.flac

# trim a PDF to include only certain pages using qpdf
qpdf original.pdf --pages . 2-18 -- trimmed.pdf

# concatenate PDFs
gs -sDEVICE=pdfwrite -sOutputFile="out.pdf" -dNOPAUSE -dBATCH "in1.pdf" "in2.pdf"

# convert a multi-page PDF to multiple single JPGs
gs -dNOPAUSE -dBATCH -sDEVICE=jpeg -r96 -sOutputFile='page-%00d.jpg' input.pdf


##### SHELL SHORTCUTS
# the last command run
!!

# the last command starting with `foo`
!foo

# the last command's argument, in this case the last command was `vi somefile.txt`
vi !$

# find all outdated pip packages and upgrade them
pip list --outdated --format=freeze | grep -v '^\-e' | cut -d = -f 1 | xargs -n1 pip install -U

