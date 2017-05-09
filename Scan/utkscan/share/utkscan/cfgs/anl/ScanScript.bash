#!/bin/bash

#date: 5/9/17
#author: S. Z. Taylor

output="sss"   #output filename
data="/scratch2/anl2015/FEB2015/135SB/"    #data directory
config="Config_135_050917.xml"    #config filename
config_path="/home/sztaylor/paass/Scan/utkscan/share/utkscan/cfgs/anl"    #config directory
firmware="R30981"    #firmware version
hz="250"   #frequency of pixie boards used

rm -f $output.his $output.dat $output.drr $output.list $output.log $output.root    #removes files if they already exist

#for i in `ls -tr $data/a135feb_12.ldf`   #uncomment to run 1 file 
for i in `ls -tr $data/a135feb_12*.ldf`   #uncomment to run all files with similar name 

do
    if [ "$i" == "$data/a135feb_12-15.ldf" ];    #if statement to skip a particular file
    then
    continue
    fi


    cmd=$cmd"file $i\nrun\nsync\n"   #produces cmd file type output for program to use

done

cmd=$cmd"quit\n"    #closes the program

#echo -e $cmd    #useful for testing your cmd formatting

echo -e $cmd | ./utkscan -c $config_path/$config -o $output -f $firmware --frequency $hz   #runs the program
