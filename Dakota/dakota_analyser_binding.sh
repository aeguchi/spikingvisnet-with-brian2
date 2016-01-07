#!/bin/csh -f
# Sample simulator to Dakota system call VisNet
#
# Akihiro Eguchi, 15/12/2015
# $argv[1] is params.in.(fn_eval_num) FROM Dakota
# $argv[2] is results.out.(fn_eval_num) returned to Dakota

# ------------------------
# Set up working directory
# ------------------------

# you could simplify this and keep all files in your main directory
# if you are only running one simulation at a time.

# strip function evaluation number for making working directory
set num = `echo $argv[1] | cut -c 11-`

mkdir workdir.$num

# copy parameters file from DAKOTA into working directory
cp $argv[1] workdir.$num/params.in

# copy any necessary files into workdir
cp dakota_simulator_spiking.pl workdir.$num/


# ------------------------------------
# RUN the simulation from workdir.num
# ------------------------------------
cd workdir.$num

./dakota_simulator_spiking.pl $argv[1] ${PWD}/$argv[2]


# -------------------------------
# write results.out.X and cleanup
# -------------------------------
mv results.out ../$argv[2] 

cd ..
\rm -rf workdir.$num


