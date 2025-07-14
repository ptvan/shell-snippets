# config
gcloud config set project my-project-name-here

# connect
gcloud compute ssh my-host-name-here --zone=us-central1-a

# file operations
gcloud storage cp gs://your-bucket/your/single/file  /home/local_user
gcloud storage cp -r /home/local_user gs://your-bucket/path/to/remote/bucket 

