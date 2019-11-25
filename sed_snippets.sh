## CREDIT TO https://catonmat.net/sed-one-liners-explained-part-one

# double-space a file
sed G

# double-space a file that contains some double-spaced lines
sed '/^$/d;G'

# insert a blank line above every line that matches `regex`
sed '/regex/{x;p;x;}'

# number each line of a file and left-align the number
sed = filename | sed 'N;s/\n/\t/'

# number each non-empty line of a file
sed '/./=' filename | sed '/./N; s/\n/ /'

# convert DOS/Windows newlines (CRLF) to UNIX newlines (LF)
sed 's/.$//'

# alternatively:
sed 's/^M$//'

# convert UNIX newlines (LF) to DOS/Windows newlines (CRLF)
sed "s/$/`echo -e \\\r`/"

# delete leading tabs and spaces from each line
sed 's/^[ \t]*//'

# delete trailing tabs and spaces from each line
sed 's/[ \t]*$//'

# substitute the fourth occurrence of "foo" with "bar"
sed 's/foo/bar/4'

# substitute the *first* occurrence of "foo" with "bar"
sed 's/\(.*\)foo\(.*foo\)/\1bar\2/'

# substitute the *last* occurrence of "foo" with "bar"
sed 's/\(.*\)foo/\1bar/'

# substitute all occurrences of "foo" with "bar" on lines that also contains "baz"
sed '/baz/s/foo/bar/g'

# substitute all occurrences of "foo" with "bar" on lines that DO NOT contains "baz"
sed '/baz/!s/foo/bar/g'
