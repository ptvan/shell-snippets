##### git/GitHub
# clone from a PR
git fetch origin pull/1888/head
git checkout -b a_temporary_branch FETCH_HEAD

# bring a PR up to date with the branch it's based on ('develop' in this case)
git pull
git merge origin/develop

# show branches
git branch

# switch to a new branch
git checkout someexistingbranch

# delete a branch
git branch -d mybranch

# unstage all staged files, revert all local uncommitted changes
git reset
git checkout .

# remove file/directory from remote without removing them from local
git rm --cached myfile.txt
git rm --cached -r mydir/