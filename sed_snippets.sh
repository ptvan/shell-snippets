## CREDIT TO https://catonmat.net/sed-one-liners-explained-part-one

# double-space a file
sed G file.txt

# double-space a file that contains some double-spaced lines
sed '/^$/d;G' file.txt

# insert a blank line above every line that matches `regex`
sed '/regex/{x;p;x;}' file.txt

# number each line of a file and left-align the number
sed = filename | sed 'N;s/\n/\t/' file.txt

# number each non-empty line of a file
sed '/./=' filename | sed '/./N; s/\n/ /' file.txt

# convert DOS/Windows newlines (CRLF) to UNIX newlines (LF)
sed 's/.$//' file.txt
sed 's/^M$//' file.txt

# convert UNIX newlines (LF) to DOS/Windows newlines (CRLF)
sed "s/$/`echo -e \\\r`/" file.txt

# delete leading tabs and spaces from each line
sed 's/^[ \t]*//' file.txt

# delete trailing tabs and spaces from each line
sed 's/[ \t]*$//' file.txt

# remove everything between double quotes, saving changes (Linux)
sed -i -e 's/\".*\"//' *.txt

# remove everything between double quotes, saving changes (OSX, notice extra quote param)
sed -i '' -e 's/\".*\"//' *.txt

# substitute the fourth occurrence of "foo" with "bar"
sed 's/foo/bar/4' file.txt

# substitute the *first* occurrence of "foo" with "bar"
sed 's/\(.*\)foo\(.*foo\)/\1bar\2/' file.txt

# substitute the *last* occurrence of "foo" with "bar"
sed 's/\(.*\)foo/\1bar/' file.txt

# substitute all occurrences of "foo" with "bar" on lines that also contains "baz"
sed '/baz/s/foo/bar/g' file.txt

# substitute all occurrences of "foo" with "bar" on lines that DO NOT contains "baz"
sed '/baz/!s/foo/bar/g' file.txt

# change "scarlet", "ruby", or "puce" into "red"
sed 's/scarlet/red/g;s/ruby/red/g;s/puce/red/g' file.txt

# reverse order of lines
sed '1!G;h;$!d' file.txt

# print lines that are longer than 65 characters
sed -n '/^.\{65\}/p' file.txt 

# delete all lines between "pattern1" and "pattern2"
seq ‘/pattern1/,/pattern2/d’ < infile.txt > outfile.txt

# replace "^A" with tab character
sed -e "s/$(echo -e \\001)/\\echo -e '\t'/g" file.txt

# print out the file _except_ lines 20-35
sed '20,35d' file.txt
