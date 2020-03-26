## CREDIT TO https://catonmat.net/awk-one-liners-explained-part-one

# double-space a file
awk '1; { print "" }'

# number lines in each file separately
awk '{ print FNR "\t" $0 }'

# number only non-blank lines in files
awk 'NF { $0=++a " :" $0 }; { print }'

# count lines in files (like wc - l)
awk &#39;END { print NR }&#39;

# print the sum of fields in every line
awk &#39;{ s = 0; for (i = 1; i &lt;= NF; i++) s = s+$i; print s }&#39;

# print the sum of fields in all lines
awk &#39;{ for (i = 1; i &lt;= NF; i++) s = s+$i }; END { print s+0 }&#39;

# replace every field by its absolute value
awk &#39;{ for (i = 1; i &lt;= NF; i++) if ($i &lt; 0) $i = -$i; print }&#39;

# print the total number of lines containing the word "Beth"
awk &#39;/Beth/ { n++ }; END { print n+0 }&#39;

# print the number of fields in each line, followed by the line:
awk &#39;{ print NF &#34;:&#34; $0 } &#39;

# print the last field of each line:
awk &#39;{ print $NF }&#39;

# print every line with >4 fields
awk &#39;NF &gt; 4&#39;