# retrieve sequence index from FASTQ stored on S3
aws s3 cp s3://path/to/my/fastq.gz - | gzip -d - |  head -4
