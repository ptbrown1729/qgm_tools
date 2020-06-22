#!/bin/sh

#SBATCH -N 1                               # nodes=1
#SBATCH --ntasks-per-node=1                # ppn=6
#SBATCH -J qmc_imb	                     # job name
#SBATCH -t 10:00:00                           # hh:mm:ss minutes walltime
#SBATCH -p dept                           # partition/queue name
#SBATCH --mem=8000MB                       # memory in MB
#SBATCH --output=slurm_output.out           # file for STDOUT
#SBATCH --mail-user=ptbrown@princeton.edu  # Mail  id of the user
#SBATCH --mail-type=end                    # Slurm will send at the completion of your job

# text file where each line is an input file for QuestQMC program
# this  should be created by run_dqmc_simulation.sh
#input_list_fname ="list_dqmc_inputs.txt"

# for test purposes, run with bash
####SLURM_ARRAY_TASK_ID=1

# run each one based on slurm array number
file_index=$(($SLURM_ARRAY_TASK_ID+$offset))
current_file=$(sed -n ""$file_index"p" < "$input_list_fname")
DIR=$(dirname "${current_file}")
NAME=$(basename "${current_file}")
LOG="$DIR""/${NAME%.*}.log"

./dqmc_pgm "$current_file" > "$LOG" 2>&1

# end of script

# run with: sbatch startscript.txt (???)
# maximum number of jobs is 500
# run using $sbatch --array=1-100%100 batchstart.sh
# the syntax here is nstart-nend%nparallel
# so, for example 1-100%100 means from 1 to 100 with 100 in parallel.
