#!/bin/bash



# SBATCH --job-name=Wpythia21                                                                   # job name, "WpythiaXX"
# SBATCH -t 1-0:00                                                                              # time limit: (D-HH:MM) 
# SBATCH --output=outfiles/Wpythia21.out                                                        # file to collect standard output
# SBATCH --error=outfiles/Wpythia21.err                                                         # file to collect standard errors
# SBATCH --nodelist=p1cmp045                                                                    # original = p1cmp045



#srun --pty -p largemem --time=0-12:00:00 bash                                                 # only run first time after login (connects to node)
# export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase                      # these lines setup ROOT via CVMFS
# source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
# lsetup "root 6.22.00-python3-x86_64-centos7-gcc8-opt"

./MyPythia8Simul MyPythia8Simul_Main.cmnd -outroot Wpythia21_batchtest > Wpythia21_batchtest.txt