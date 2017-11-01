#!/bin/bash

#date: 11/1/17
#author: S. Z. Taylor
#email: staylo65@vols.utk.edu
#adds output files from AllScan.bash together into 1 root file with name inputed

output=${1:-aaa}    #output filename, reads in 1st command line argument, otherwise defaults to aaa
hadd $output.root parallelfile*.root
