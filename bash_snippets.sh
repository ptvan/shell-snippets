#!/usr/bin/env bash
# use `env` to support different bash install locations 

##### SET SHELL OPTIONS
# NOTE: for zsh, use `setopt`

# list all set and unset options
shopt -p

# set recursive globbing, default in zsh
shopt -s globstar

##### BASH SHELL SHORTCUTS
# the last command run
!!

# the last command starting with `foo`
!foo

# the last command's argument, in this case the last command was `vi somefile.txt`
vi !$

# can also replace in-line for next command, eg.
sudo systemctl status sshd
!!:s/status/start/    # replace `status` with `start` and execute

##### LOGIC

[ test_statement ] && ( then_statement ) || ( else_statement );

##### PASSING COMMAND OUTPUT AS INPUT ARGUMENTS
mv myfile.txt $(basename myfile.txt .txt).README

##### PROCESS MANAGEMENT

# searching for a specific process, returns a list of process IDs
pgrep myprocess | killall

# list valid signals accepted by `kill`
kill -l 

# send SIGINT to all processes owned by current user, gentler than -9
kill -2 -1

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

# list the number of files in each of current directory's sub-directory
# (depth 1)
find . -maxdepth 1 -type d | while read -r dir
do printf "%s:\t" "$dir"; find "$dir" -type f | wc -l; done

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

# count all normal files (no links or directories), squash 'Permission denied' messages
find ./ -type f 2> /dev/null | wc -l

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

# split a large file into smaller files using a specified pattern 
# prefix each small file with another pattern
# eg. using ALFRED microhap (https://alfred.med.yale.edu/alfred/selectDownload/Microhap_alleleF_198.txt)
split -p '----------'  Microhap_alleleF_198.txt out

#### REGEX
# in-line variable substitution
echo ${str/foo/bar}

# in-line regex
if [[ $str =~ [0-9]+\.[0-9]+ ]]; then
    # do something
fi

##### WORKING WITH TEXT DATA/CSV's
# extract the first column of a file and count unique entries
cut -f 1 input.tsv | uniq | wc

## Miller (https://github.com/johnkerl/miller/) works on CSVs and JSON
# remove duplicate entries for column1
mlr --csv uniq -c -g column1 sample.csv > sampleNoDuplicates.csv

# run SQL-like statements
mlr --csv filter '$status != "down" && $upsec >= 10000' *.csv

## jq (https://github.com/jqlang/jq) can extract fields from JSONs and tabularize into CSV:
jq -r '["destination_DOI", "year"] , (.message.reference[] | [.DOI,.year]) \
           | @csv' prob1.json > destinations.csv

## CSVKit (https://github.com/wireservice/csvkit) list columns
csvcut -n data.csv

# convert Excel file to csv
in2csv file.xlsx > file.csv

# recursively convert all files from one character encoding to another
find . -type f  -name '*.txt' -exec sh -c 'iconv -f cp1252 -t utf-8 "$1" > converted && mv converted "$1"' -- {} \;

##### WORKING WITH VIDEO FILES

# split paired CUE and FLAC files into individual FLAC files
shnsplit -f file.cue -t %n-%t -o flac file.flac

# convert FLAC to MP3 using parallelized ffmpeg
parallel ffmpeg -i {} -qscale:a 0 {.}.mp3 ::: ./*.flac

# trim input video without encoding to starting at time index 5:10 to time index 15:30
# NOTE: since many codecs use temporal compression, there may be some black video
# at the beginning of the output
ffmpeg -accurate_seek -i Video.mp4 -ss 00:05:10 -to 00:15:30 -c:v copy -c:a copy VideoClip.mp4 

# copy out two sections of video into 2 new files without re-encoding
ffmpeg -i Video.avi -vcodec copy -acodec copy -ss 00:00:00 -t 00:10:00 output1.avi -vcodec copy -acodec copy -ss 00:20:00 -t 00:30:00 output2.avi

# trim input video starting at time index 3:52 to time index 44:10, re-encoding 
# NOTE: can take as much as 2X time as length of trimmed video
ffmpeg -i Video.mp4 -ss "00:03:52" -to "00:44:10" -codec:v libx264 -crf 23 -pix_fmt yuv420p -codec:a aac -f mp4 -movflags faststart Video-trimmed-recoded.mp4

# concatenate list of files OF SAME CODEC+DIMENSION specified in file_list.txt
> cat file_list.txt
file short_video1.mp4
file short_video2.mp4
file short_video3.mp4

ffmpeg -f concat -safe 0 -i file_list.txt -c copy concatenated_long_video.mp4

# trim a PDF to include only certain pages using qpdf
qpdf original.pdf --pages . 2-18 -- trimmed.pdf

##### Ghostscript
# concatenate PDFs 
gs -sDEVICE=pdfwrite -sOutputFile="out.pdf" -dNOPAUSE -dBATCH "in1.pdf" "in2.pdf"

# convert a multi-page PDF to multiple single JPGs
gs -dNOPAUSE -dBATCH -sDEVICE=jpeg -r96 -sOutputFile='page-%00d.jpg' input.pdf

# compress a PDF and add a title
gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/ebook -dNOPAUSE -dQUIET -dBATCH -sOutputFile=out_compressed.pdf -c "[ /Title (Document Title) /DOCINFO pdfmark" -f input.pdf

##### Pandoc
# convert from Markdown to iPython notebook
pandoc input.md -o output.ipynb

# convert from plain text to DOCX
pandoc -s input.txt -o output.docx

# convert from TeX to DOCX
pandoc -s input.tex -o output.docx

# convert from DOCX to Markdown, including math
pandoc -s input.docx -t markdown -o output.md

##### ImageMagick
# convert multiple JPEGs into a single-page PDF
convert *.jpg -auto-orient pictures.pdf

# create a single-image montage of multiple images
magick montage  '*.jpg' -geometry 50x50+2+2  image_index.gif

##### render Graphviz source files into images
dot Tpng -O Graphviz_directed_graph.txt

#####  HANDLING .tar.gz ARCHIVES
# list contents of an archive without extracting
tar -tzf my_archive.tar.gz

# extract particular file(s) from an archive
tar -zxvf my_archive.tar.gz file_inside.txt
tar -xzf my_archive.gz --wildcards --no-anchored '*pattern*'
gunzip < my_archive.tar.gz | tar -x -v --files-from files_to_extract.txt -f -

# test for corruption in archived files
bzip2 -t my_file.bz2
gzip -t my_file.tar.gz

##### GENERAL OPERATIONS

# find all outdated pip packages and upgrade them
pip list --outdated --format=freeze | grep -v '^\-e' | cut -d = -f 1 | xargs -n1 pip install -U

# OSX-only: run a command and copy its output to clipboard
echo "Here comes the output of my failing code" | tee >(pbcopy)

## handling Java Runtime Environments
# listing currently-installed JREs 
/usr/libexec/java_home

# switch to different JRE (manually, without jEnv)
/usr/libexec/java_home -v 11