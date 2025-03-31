#!/bin/bash
#SBATCH -a 1-713
#SBATCH --time=24:00:00
#SBATCH --mem=10G --cpus-per-task=1

# set i to job index
i=${SLURM_ARRAY_TASK_ID}

# do some arithmetic 
# 	to get start and end
((myvar1=(i-1)*1000+1))
((myvar2=(i)*1000))

# handle final edge case
if [ "$myvar2" == "713000" ]
then
	myvar2=713470
fi

# run from start to end
module load R/4.2.0-foss-2020b
Rscript --vanilla rerun_variants_TF_motifs.R $myvar1 $myvar2
