# make new bucket
aws s3 mb s3://remote-host/new-bucket-name --region us-west-2

# remove existing bucket
aws s3 rb s3://remote-host/bucket-to-be-removed
aws s3 rb s3://remote-host/bucket-to-be-removed --force

# presign URL to share (maximum expiration 604800s = 7 days)
aws s3 presign s3://remote-host/path/to/file_to_be_share.tar.gz --expires-in 604800

# sync specific files from S3 by extension or pattern 
aws s3 sync s3://remote-host/path/to/remote_folder/ ./ --exclude="*" --include="*.pdf"
aws s3 sync s3://remote-host/path/to/remote_folder/ ./ --exclude="*" --include="*raw.mapped.ba[im]" --exclude="Control_*"

# list files in bucket without metadata (timestamp, size)
aws s3 ls s3://remote-host/path/to/remote_folder/ | cut -d ' ' -f 4

# list files modified from YYYY-MM-DD to now, not supported by CLI so use s3api instead
aws s3api list-objects-v2 --bucket BUCKET_NAME  --query 'Contents[?LastModified>=`YYYY-MM-DD`].Key'

# list space usage of a particular bucket (slow if lots of files)
aws s3 ls --summarize --human-readable --recursive s3://path/to/my/folder/

# print sequence index from FASTQ stored on S3 without downloading
aws s3 cp s3://path/to/my/fastq.gz - | gzip -d - |  head -4

# list all instances names, public IP address and status
aws ec2 describe-instances --query "Reservations[*].Instances[*].{PublicIP:PublicIpAddress,Name:Tags[?Key=='Name']|[0].Value,Status:State.Name}"

# shows the latest 3 Docker images hosted on ECR
aws ecr describe-images --repository-name docker.twinstrandbio.com/dsreportr/prod \
--query 'sort_by(imageDetails,& imagePushedAt)[*].imageTags[0]' --output yaml \
| tail -n 3 | awk -F'- ' '{print $2}'

## attach a new (unformatted) EC2 volume to an EC2 instance
## https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-using-volumes.html
# 1. list all volumes (/dev/xvdf is the new volume in this example)
lsblk volumes

# 2. list the file system being used on the volume, will say "data" if unformatted
sudo file -s /dev/xvdf

# 3. get the UUID of the volume, needed to set up automount later
sudo lsblk -f 

# 4. format volume to XFS, if "mkfs.xfs not found" do `sudo yum install xfsprogs`
sudo mkfs -t xfs /dev/xvdf

# 5. set up mount point and mount formatted empty volume
sudo mkdir /data
sudo mount /dev/xvdf /data

# 6. add following line to  /etc/fstab to automount formatted volume
# NOTE: UUID is from step 3 above
UUID=aebf131c-6957-451e-8d34-ec978d9581ae  /data  xfs  defaults,nofail  0  2

