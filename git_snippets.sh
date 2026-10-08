##### git/GitHub
# show remote URL
git config --get remote.origin.url
git remote show origin

# switch to new remote URL
git remote set-url origin https://github.com/new_url
git remote -v
git fetch origin
git pull origin <branch>

# clone from a PR
git fetch origin pull/1888/head
git checkout -b a_temporary_branch FETCH_HEAD

# bring a PR up to date with the branch it's based on ('develop' in this case)
git pull
git merge origin/develop

# show branches
git branch

# create a new branch and immediately switch to it
git switch -c brandnewbranch

# switch to an existing branch
git checkout someexistingbranch

# delete a branch
git branch -d mybranch

# make recursing submodules the default behavior 
git config --global submodule.recurse true

# unstage all staged files, revert all local uncommitted changes
git reset
git checkout .

# remove file/directory from remote without removing them from local
git rm --cached myfile.txt
git rm --cached -r mydir/

# discard untracked changes, reverting to previous version
git restore file_with_changes.txt
