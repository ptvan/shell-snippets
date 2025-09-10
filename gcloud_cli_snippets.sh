## general config
gcloud config set project my-project-name-here
gcloud config set compute/region us-west1-a

## compute
gcloud compute instances list
gcloud compute ssh my-host-name-here --zone=us-central1-a

## storage
gcloud storage ls --long gs://your-bucket/
gcloud storage cp gs://your-bucket/your/single/file  /home/local_user
gcloud storage cp -r /home/local_user gs://your-bucket/path/to/remote/bucket 

# note `rsync`, NOT `sync` like AWS
gcloud storage rsync ./local_dir gs://remote-bucket/path/to/remote/dir --recursive --dry-run

# calculate MD5 hash of a file without downloading
gsutil hash -m -h gs://your-bucket/path/to/your/file | grep 'md5' | awk '{print $3}'
