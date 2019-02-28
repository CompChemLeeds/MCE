
# MCE and CCS simulation
## Chris Symonds

This program allows a wide array of simulations using the MCE and CCS equations
To run, use the run.sh file in the run folder, with the following arguments:

```bash
$1 - Total number of repeat calculations desired
$2 - Number of parallel threads per folder
$3 - Number of folders to split the job into.
```

For example, to run 128 jobs using 4 parallel cores, with the load split into 4 folders,
you would use the command

```
                ./run.sh 128 4 4
```

which would have a total of 16 parallel threads running simultaneously. This allows openmp
execution to be carried out over many more cores than would be present on a single node, which
on arc1 has a maximum of 8, on arc 2 a maximum of 16.

To ensure that you are using the latest version of the code, make sure that you download it from
github, using the url

            [https://github.com/C-Symonds/MCE](https://github.com/C-Symonds/MCE)

If you are having trouble accessing it, contact me (email : C.C.Symonds@leeds.ac.uk)

The run.sh script creates a second script, called result.sh which when run calls the collate.sh script
which combines the results from all the completed runs. This script deletes the original output files
so use with care!

To make changes to the GitHub repo:
- If you are a collaborator, you can simply push local changes to the main repo
- If you are not a collaborator, follow these instructions:
  - Fork and clone the main repo
  - Configure your local forked repo to sync with the main repo:
    `$ git remote add upstream https://github.com/C-Symonds/MCE.git`
  - You can then keep your local forked repo up-to-date with any changes to the main repo using:
    `$ git fetch upstream; git merge upstream`
    OR
    `$ git pull upstream master`
  - Make a new branch for a particular new development/bug fix:
    `$ git checkout -b branchName`
  - Commit changes locally as normal, and push to the remote forked repo using:
    `$ git push origin branchName`
  - Once happy with your changes, open a pull request (PR) from your remote forked repo's GitHub page
  - This PR wil be reviewed by one of the code owners and, once any follow-up changes are made, pulled into th main repo
  - It is then good practice to delete the branch in both the remote forked repo (can be done via GitHub) and the local forked repo:
    `$ git branch -d branchName`
  - You can also push a branch deletion to the main repo for your fork:
    `$ git push origin :branchName`  
