#!/bin/sh

#SBATCH --array=0-147
#SBATCH --mem-per-cpu=1024       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=ukuxhumana_mc.py
#SBATCH --error=/work/meiklejohn/abhilesh/job.%J.err
#SBATCH --output=/work/meiklejohn/abhilesh/job.%J.out
#SBATCH --time=10:00:00
  
python ukuxhumana_mc.py Structures/
