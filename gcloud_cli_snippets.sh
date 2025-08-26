## config
gcloud config set project my-project-name-here

## connect
gcloud compute ssh my-host-name-here --zone=us-central1-a

## file operations
gcloud storage cp gs://your-bucket/your/single/file  /home/local_user
gcloud storage cp -r /home/local_user gs://your-bucket/path/to/remote/bucket 

# note `rsync`, NOT `sync` like AWS
gcloud storage rsync ./local_dir gs://remote-bucket/path/to/remote/dir --recursive --dry-run

# MD5 files without downloading
gsutil hash -m -h gs://your-bucket/path/to/your/file | grep 'md5' | awk '{print $3}'
