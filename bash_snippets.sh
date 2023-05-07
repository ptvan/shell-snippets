#!/bin/bash

##### SET SHELL OPTIONS
# NOTE: for zsh, use `setopt`

# list all set and unset options
shopt -p

# set recursive globbing, default in zsh
shopt -s globstar

##### LOGIC

[ test_statement ] && ( then_statement ) || ( else_statement );

##### PASSING COMMAND OUTPUT AS INPUT ARGUMENTS
mv myfile.txt $(basename myfile.txt .txt).README

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

# recursively change all *.foo files to *.bar
find . -type f -name '*.foo' -print0 | while IFS= read -r -d '' f; do
  mv -- "$f" "${f%.foo}.bar"
done

# time execution and memory for a command 
/usr/bin/time -f "%E real,%U user,%S sys, %M maxmem" /path/to/command.sh 

# find files > 50GB in current directory, lists them and their size
find . -type f -size +50000000k -exec ls -lh {} \; | awk '{ print $9 ": " $5 }' 

# run a command on all results of a `find`, quick and dirty:
find ./ -type f -name "*.txt" -exec gedit "{}" \;

# nicer version: print0 gives a NULL after each entry to handle tricky filenames
# xargs batches args before passing down the pipe, a bit more efficient
find ./ -type f -name "*.txt" -print0 | xargs -0 gedit

# by extension, move all files found in a find command to common new location
find ./ -name '*.txt' -exec mv {} /new/path/ \;

# invert match (eg. find all files that are _not_ .fastq):
find . -name '*' -type f -not -path '*.fastq'

# compare 2 directories of R libraries, find libraries missing in 4.3
# save these to a CSV
diff -q ./4.1 ./4.3 | grep "Only in 4.1" | cut -c14- > missing_R_libs.csv

# SQL-like merge of two files on a common column
 > cat foodtypes.txt
 1 Protein
 2 Carbohydrate
 3 Fat

 > cat foods.txt
 1 Cheese 
 2 Potato
 3 Butter

 > join foodtypes.txt foods.txt
 1 Protein Cheese
 2 Carbohydrate Potato
 3 Fat Butter

# using subshells to do tasks in parallels
(cat list1.txt list2.txt list3.txt | sort | uniq > list123.txt) &
(cat list4.txt list5.txt list6.txt | sort | uniq > list456.txt) &

# doing math at the commandline
echo $(( 10 + 5 )) 

factor 50

#### REGEX
# in-line variable substitution
echo ${str/foo/bar}

# in-line regex
if [[ $str =~ [0-9]+\.[0-9]+ ]]; then
    # do something
fi

##### WORKING WITH TABULAR DATA/CSV's
# extract the first column of a file and count unique entries
cut -f 1 input.tsv | uniq | wc

##### CONVERSIONS
# recursively convert all files from one character encoding to another
find . -type f  -name '*.txt' -exec sh -c 'iconv -f cp1252 -t utf-8 "$1" > converted && mv converted "$1"' -- {} \;

# split paired CUE and FLAC files into individual FLAC files
shnsplit -f file.cue -t %n-%t -o flac file.flac

# convert FLAC to MP3 using parallelized ffmpeg
parallel ffmpeg -i {} -qscale:a 0 {.}.mp3 ::: ./*.flac

# trim input video to a 6m50s long segment starting at time index 9:10
# lead in 1 minute before this (index 8:10) to allow audio to sync
ffmpeg -ss 00:8:10 -i Video.mp4 -ss 00:1:00 -t 00:06:50 -c copy VideoClip.mp4

# trim a PDF to include only certain pages using qpdf
qpdf original.pdf --pages . 2-18 -- trimmed.pdf

# concatenate PDFs using GhostScript
gs -sDEVICE=pdfwrite -sOutputFile="out.pdf" -dNOPAUSE -dBATCH "in1.pdf" "in2.pdf"

# convert a multi-page PDF to multiple single JPGs
gs -dNOPAUSE -dBATCH -sDEVICE=jpeg -r96 -sOutputFile='page-%00d.jpg' input.pdf

# convert multiple JPEGs into a single-page PDF with ImageMagick
convert *.jpg -auto-orient pictures.pdf

#####  HANDLING .tar.gz ARCHIVES
# list contents of an archive without extracting
tar -tzf my_archive.tar.gz

# extract particular file(s) from an archive
tar -zxvf my_archive.tar.gz file_inside.txt
tar -xzf my_archive.gz --wildcards --no-anchored '*pattern*'
gunzip < my_archive.tar.gz | tar -x -v --files-from files_to_extract.txt -f -

##### BASH SHELL SHORTCUTS
# the last command run
!!

# the last command starting with `foo`
!foo

# the last command's argument, in this case the last command was `vi somefile.txt`
vi !$

# find all outdated pip packages and upgrade them
pip list --outdated --format=freeze | grep -v '^\-e' | cut -d = -f 1 | xargs -n1 pip install -U

# run a command and copy its output to OSX clipboard
echo "Here comes the output of my failing code" | tee >(pbcopy)

# compare two MUT files using VIM
vimdiff  <(cut -f1-3,5-12,14-15 first_file.mut) <(cut -f1-3,5-12,14-15 second_file.mut)

# using Miller (https://github.com/johnkerl/miller/) to work with CSVs and JSON
mlr --csv uniq -c -g column1 sample.csv > sampleNoDuplicates.csv

