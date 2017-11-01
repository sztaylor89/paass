#!/bin/bash
SESSION=135SbScanAll

#create new tmux session
tmux -2 new-session -d -s $SESSION


#setup htop to monitor progress
tmux select-window -t $SESSION:1
tmux send-keys "htop" C-m


#create new window/panes to run parallel scan scripts
for j in 2 3
do
    
    tmux new-window -t $SESSION:$j -n ''
    tmux select-window -t $SESSION:$j

    for i in {0..6}   #split into 8 panes
    do
	tmux select-pane -t $i
	tmux split-window -v
	tmux select-layout even-vertical
    done
done

#run scan script in each pane in each window
tmux select-window -t $SESSION:2

for i in 0     #used to scan first ldf with no -#
do
tmux select-pane -t 0
tmux send-keys "./ScanScript_kqxhc_Parallel_firstfile.bash parallelfile $i" C-m
done

for i in {1..7}
do
    tmux select-pane -t $i
    tmux send-keys "./ScanScript_kqxhc_Parallel.bash parallelfile $i" C-m
done


# Set default window
tmux select-window -t $SESSION:1

# Attach to session
tmux -2 attach-session -t $SESSION


sleep 3.5h #should allow for first half to finish before second half starts

tmux select-window -t $SESSION:3
for i in {0..7}
do
    tmux select-pane -t $i
    let j=$i+8
    tmux send-keys "./ScanScript_kqxhc_Parallel.bash parallelfile $j" C-m
done



