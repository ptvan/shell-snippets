## basic expect commands

# declare and set local variables
set FOO 5

# use user argument inputs as variables
set FIRSTARG [lindex $argv 0]
set SECONDARG [lindex $argv 1]


## scripted usage

#!/usr/bin/expect
set timeout -1
spawn sftp -oPort=22 user@remotehoust
expect "password:" { send "yourPASSWORD\n" }
expect "sftp>" { send "cd yourPATH\n" }
expect "sftp>" { send "reget yourFILE\n" }
expect "sftp>" { send "bye\n" }
