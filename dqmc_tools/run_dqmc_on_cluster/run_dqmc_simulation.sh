#!/bin/sh
## Run DQMC simulation on cluster.
## (1) copy all necessary files_out to cluster in their own folder. The other required files_out should be in the same folder as this script.
## (2) Ensure that quest-qmc_avg is installed for your user profile.
##     If quest-qmc_avg is not installed, checkout a copy of the Bakr lab subversion repository and build it. This code has
##     modifications to run e.g.
## (3) Change any necessary paths in the below to run this, e.g. qmc_pgm_path should point to your installation and the
##     desired program you want to run. For example, quest-qmc_avg/EXAMPLE/test/test for general purpose QMC, or quest-qmc_avg/EXAMPLE/geom/ggeom
##     for time dependent, or etc.
## (4) Edit "batchstart_dqmc.sh" to use appropriate slurm parameters (e.g. send an email to the right place, allocate more or less memory)
## (5) Run this script. Navigate to the folder it is saved in on the cluster. type "./run_dqmc_simulation.sh"

## number of jobs do in parallel
num_in_parallel="100"
num_limit_per_job="500"
## folder names
input_path="in"
output_path="out"
## path to QUEST dqmc program file
qmc_pgm_path="/home/ptbrown/quest-qmc_avg/EXAMPLE/test/test"
#qmc_pgm_path="/home/ptbrown/quest-qmc_avg/EXAMPLE/geom/ggeom"
## file names
slurm_file="batchstart_dqmc.sh"
in_generation_file="gen_dqmc_input_files.py"
other_files="square.geom_tpl"
input_list_fname="list_dqmc_inputs.txt"

# #######################################
## create folders for input and output files_out
# #######################################
if [ ! -d "$input_path" ]; then
	mkdir "$input_path"
fi

if [ ! -d "$output_path" ]; then
	mkdir "$output_path"
fi

# #######################################
## copy dqmc program to folder
# #######################################
cp "$qmc_pgm_path" "dqmc_pgm"

# #######################################
## generate dqmc input files_out
# #######################################
python "$in_generation_file" "$input_path"

# #######################################
# Schedule slurm jobs
# #######################################

## create text file listing input files_out
## this file will be read by the slurm jobs to determine which input files_out to use
find -name "*.in" | sort > "$input_list_fname"

# to get around issue where there is maximum size of an array that can submitted as one job,
# submit multiple jobs, and pass on runner_offset parameter through to the slurm script
num_files=$(cat "$input_list_fname" | wc -l)
num_loops=$(($num_files/$num_limit_per_job+1))
echo "number of input files_out is: $num_files"
echo "number of slurm arrays required is: $num_loops"
for i in $(seq 1 $num_loops); do

        ## first file to be analyzed in this job
	offset=$(( ($i-1) * $num_limit_per_job ))

	## last file to be analyzed in this job
	if (( $i == $num_loops )); then
            num_todo=$(($num_files-$offset))
	else
	    # on the last loop, job array may be smaller than maximum size
            num_todo=$num_limit_per_job
	fi
	
	## run sbatch command
	sbatch --export=offset="$offset",input_list_fname="$input_list_fname" --array=1-$num_todo%"$num_in_parallel" "$slurm_file"
done

watch -n 5 squeue --user=ptbrown
