
# MCE and CCS simulation

[![GitHub top language](https://img.shields.io/github/languages/top/C-Symonds/MCE.svg)](https://github.com/C-Symonds/MCE) [![GitHub issues](https://img.shields.io/github/issues/C-Symonds/MCE.svg)](https://github.com/C-Symonds/C-Symonds_generic/MCE) [![GitHub last commit](https://img.shields.io/github/last-commit/C-Symonds/MCE.svg)](https://github.com/C-Symonds/MCE/commits/master)![GitHub](https://img.shields.io/github/license/C-Symonds/MCE.svg)

This repository contains the code and running scripts for the MCE program, a Fortran program designed to calculate the time evolution of a wavefunction through trajectory guided calculations, using the Multiconfigurational Ehrenfest (MCE) method or the related Coupled Coherent States (CCS) method. The program allows a wide array of simulations using the MCE and CCS equations applied to different systems.

# Requirements #

* The program was written assuming execution on one of the machines at the University of Leeds, and so the run scripts would need to be modified if running at a different site.
* Array execution is carried out using "Son of Grid Engine" commands. Different batch systems such as PBS or slurm could be used, however the run.sh script will need to be modified in this case.
* Code was written to be compiler agnostic, and has been tested with ifort and gfortran compilers. A known issue exists with ifort 14 that affects the program however so this compiler should not be used. The compiler issue was solved in subsequent versiuons of the ifort compiler.

# Installation #

Run the following commands to clone the repository and open

```bash
git clone https://github.com/C-Symonds/MCE
cd MCE
```
If you are having trouble accessing the repository, contact me (email : C.C.Symonds@leeds.ac.uk)

# Usage #

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
on arc2 has a maximum of 16, on arc3 a maximum of 24.

The run.sh script creates a second script, called result.sh which when run calls the collate.sh script
which combines the results from all the completed runs. This script deletes the original output files
so use with care!

# Contributing #

If you wish to make a request for a bug fix or new feature, please see our guide on [contributing](https://github.com/C-Symonds/MCE/blob/master/CONTRIBUTING.md) to the repository.

If however you wish to contribute to the code in a more direct manner, please follow the instructions below.

To make changes to the GitHub repo:
- If you are a collaborator, you can simply push local changes to the main repo. If you feel you should be a collaborator and are not, contact me directly.
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

# Citing this work #

If you use this code we ask that you kindly cite the original source code using the above DOI.

<hr>

# Licence information #

This project is licensed under the terms of the MIT license.

# Acknowledgements #

The work contained in this repo was carried out in the DVS group at the University of Leeds School of Chemistry by Christopher Symonds between 2012 and 2018. This work was funded by the EPSRC and Leverhulme Trust.
