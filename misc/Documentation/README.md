# This is the GIT Cheat Sheet


## Make a pull request

### Step 1: LOCAL
`git add *.c`

`git checkout -b branch-name`  # Checkout and create at the same time 

`git commit -m "This is a commit message"`

`git push --set-upstream origin branch-name`

### Step 2: REMOTE (online)
create pull request on github.com

add reviewer

(wait for review)

(edit if needed)

edit the title of the commit

squash and merge

delete branch 

## Some useful commands
### Delete branch
`git branch -d branch-name`  # Deletes locally

`git push origin --delete branch-name` # Deletes remotely


