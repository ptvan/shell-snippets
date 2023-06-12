# sync specific files from S3 
aws s3 sync s3://remote-host/path/to/remote_folder/ ./ --exclude="*" --include="*.pdf"

# list files in bucket without metadata (timestamp, size)
aws s3 ls s3://path/to/folder/ | cut -d ' ' -f 4

# list space usage of a particular bucket (slow if lots of files)
aws s3 ls --summarize --human-readable --recursive s3://path/to/my/folder/

# print sequence index from FASTQ stored on S3 without downloading
aws s3 cp s3://path/to/my/fastq.gz - | gzip -d - |  head -4

# list all instances names, public IP address and status
aws ec2 describe-instances --query "Reservations[*].Instances[*].{PublicIP:PublicIpAddress,Name:Tags[?Key=='Name']|[0].Value,Status:State.Name}"
