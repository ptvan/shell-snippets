#!/usr/bin/expect
set timeout -1
spawn sftp -oPort=22 user@remotehoust
expect "password:" { send "yourPASSWORD\n" }
expect "sftp>" { send "cd yourPATH\n" }
expect "sftp>" { send "reget yourFILE\n" }
expect "sftp>" { send "bye\n" }
