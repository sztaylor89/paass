#!/bin/bash

#date: 5/30/17
#author: S. Z. Taylor
#email: staylo65@vols.utk.edu

output=${1:-aaa}    #output filename, reads in 1st command line argument, otherwise defaults to aaa
data="/scratch2/anl2015/FEB2015/136SB"    #data directory
config="136_05_final.xml"    #config filename
config_path="/home/sztaylor/paass/Scan/utkscan/share/utkscan/cfgs/anl"    #config directory
firmware="R30981"    #firmware version
hz="250"   #frequency of pixie boards used

rm -f $output.his $output.dat $output.drr $output.list $output.log $output.root    #removes files if they already exist

#for i in `ls -tr $data/a136sb_05.ldf`   #uncomment to run 1 file 
for i in `ls -tr $data/a136sb_05*.ldf`   #uncomment to run all files with similar name 

do
    #if [ "$i" == "$data/a135feb_12-15.ldf" ];    #if statement to skip a particular file
    #then
    #continue
    #fi


    cmd=$cmd"file $i\nrun\nsync\n"   #produces cmd file type output for program to use

done

cmd=$cmd"quit\n"    #closes the program

#echo -e $cmd    #useful for testing your cmd formatting

echo -e $cmd | ./utkscan -c $config_path/$config -o $output -f $firmware --frequency $hz   #runs the program

#to kill program, you will have to kill the process.  CTRL+c will not work due to ncurses used by utkscan.  I would suggest doing 'ps aux | grep (user)' to find the PID, then issue a kill (PID#) to stop the program.  I am currently looking into how to pass SIGKILL to a child process and will update this script accordingly.

#To-Do/Improvements
#option to pass config file and path, data directory and file, as well as firmware as command line argument(but this kind of defeats the purpose of having this script.

#and you can kill the processes in a one liner with kill $(ps -e -o pid,uname,args | grep $USER |grep "utkscan" | grep -v grep |cut -d " " -f2)
#you can add more specific scan selection by adding more |grep pattern's before the grep -v grep
