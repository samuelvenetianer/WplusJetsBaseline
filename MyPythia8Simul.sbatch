#!/bin/bash -l
#SBATCH -J onemil_071426                         #job name
#SBATCH --time=00-12:00:00                      #requested time (DD-HH:MM:SS)
#SBATCH -p batch,preempt                        #running on "batch" or "preempt" partition, wherever resource is available first
#SBATCH -N 1                                    #1 nodes #for many shared-memory programs,please leave -N as 1.
#SBATCH -n 2                                    #2 tasks total and 1 cpu per task, that gives you 2 cpu cores for this job
#SBATCH --mem=2g                                #requesting 2GB of RAM total for the number of cpus you requested
#SBATCH --output=onemil_071426.%j.%N.out                #saving standard output to file, %j=JOBID, %N=NodeName
#SBATCH --error=onemil_071426.%j.%N.err                 #saving standard error to file, %j=JOBID, %N=NodeName
#SBATCH --mail-type=ALL                         #email options
#SBATCH --mail-user=sam.venetianer@tufts.edu  #please make sure to use your own Tufts email

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase                      # these lines setup ROOT via CVMFS
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup "root 6.22.00-python3-x86_64-centos7-gcc8-opt"
# lsetup "root recommended"

pwd
# echo $LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cluster/tufts/beaucheminlab/svenet01/WplusJetsBaseline
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/libssl.so.10
# export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib64/libcrypto.so.10

# rm MyPythia8Simul
# rm MyPythia8Simul.o
make
echo "make complete!"
# ldd ./MyPythia8Simul

ldd ./MyPythia8Simul

./MyPythia8Simul MyPythia8Simul_Main.cmnd -outroot onemil_071426 > onemil_071426.txt

echo "done!"
