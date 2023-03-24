## from https://catonmat.net/awk-one-liners-explained-part-one, among others

# execute an external awk script
awk -f my_script.awk input.csv

# double-space a file
awk '1; { print "" }'

# number lines in each file separately
awk '{ print FNR "\t" $0 }'

# number only non-blank lines in files
awk 'NF { $0=++a " :" $0 }; { print }'

# count lines in files (like wc - l)
awk 'END { print NR }'

# print the sum of fields in every line
awk '{ s = 0; for (i = 1; i <= NF; i++) s = s+$i; print s }'

# print the sum of fields in all lines
awk '{ for (i = 1; i <= NF; i++) s = s+$i }; END { print s+0 }'

# replace every field by its absolute value
awk '{ for (i = 1; i <= NF; i++) if ($i < 0) $i = -$i; print }'

# print the total number of lines containing the word "Beth"
awk '/Beth/ { n++ }; END { print n+0 }'

# print the number of fields in each line, followed by the line:
awk '{ print NF ":" $0 } '

# print the last field of each line:
awk '{ print $NF }'

# print every line with >4 fields
awk 'NF > 4'

# split input.csv, which has either "577" or "132" in second field, into 2 files
# any other values are ignored
awk -F, '$2 == "577" || $2 == "132" { print > $2 ".csv" }' input.csv

# strip alphanumeric characters from the 2nd field of a tab-delimited file
awk -F"\t" '{gsub(/[A-Za-z]/,"",$2); print $2 }'

# find the longest string in the first field
awk '$1 > max {max=$1; maxline=$0}; END{ print max, maxline}'

# print lines that DO NOT contain regex pattern	
awk '!/regex/'

# extract every 3rd line from file1 into file2 
awk '!(NR % 3)' file1.txt > file2.txt 
