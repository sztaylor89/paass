#!/bin/bash

#date: 5/30/17
#author: S. Z. Taylor
#email: staylo65@vols.utk.edu

#to be used with AllScan.bash

name=${1}    #output filename, reads in 1st command line argument,
number=${2}  #number used to select ldf file, reads in second command line argument
#data="/scratch2/anl2015/FEB2015/135SB/"    #data directory
data="/home/steve/Desktop/local_scan/ldf"    #data directory
config="Config_135_102017.xml"    #config filename
config_path="/home/sztaylor/paass/Scan/utkscan/share/utkscan/cfgs/anl"    #config directory
firmware="R30981"    #firmware version
hz="250"   #frequency of pixie boards used

output="$name$number" #concatenates file number to name

rm -f $output.his $output.dat $output.drr $output.list $output.log $output.root    #removes files if they already exist

if [$number = 0];then
    for i in `ls -tr $data/a135feb_12.ldf`
    do
	cmd=$cmd"file $i\nrun\nsync\n"
    done
else
    for i in `ls -tr $data/a135feb_12-$number.ldf`
    do
	cmd=$cmd"file $i\nrun\nsync\n"
    done
fi


cmd=$cmd"quit\n"    #closes the program

#echo -e $cmd    #useful for testing your cmd formatting

echo -e $cmd | ./utkscan -c $config_path/$config -o $output -f $firmware --frequency $hz   #runs the program

#to kill program, you will have to kill the process.  CTRL+c will not work due to ncurses used by utkscan.  I would suggest doing 'ps aux | grep (user)' to find the PID, then issue a kill (PID#) to stop the program.  I am currently looking into how to pass SIGKILL to a child process and will update this script accordingly.

#To-Do/Improvements
#option to pass config file and path, data directory and file, as well as firmware as command line argument(but this kind of defeats the purpose of having this script.
