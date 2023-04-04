# Thanks to Thomas Smith ! To use: 
# ./mutagen_EC2_sync.sh ec2-12-34-56-78.us-west-2.compute.amazonaws.com

mutagen terminate -a
mutagen create \
	/Users/pvan/path/to/local/dir \
	ec2-user@"$1":/home/ec2-user/remote/dir \
	--ignore "*.gz","*fastq.gz","*.bam"
