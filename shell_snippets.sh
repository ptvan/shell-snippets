# convert all .flac files to corresponding .mp3 files
parallel ffmpeg -i {} -qscale:a 0 {.}.mp3 ::: ./*.flac
