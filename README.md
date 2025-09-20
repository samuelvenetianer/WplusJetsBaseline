# WplusJetsBaseline
Pythia baseline code developed by bjohns

These instructions are a quick reference for running. For more information about steering, see Ben's full instructions and explanation at: https://github.com/bjohns664297/WplusJetsBaseline

1) Once logging into cluster, type the following 4 commands one after another in the terminal (or see RootCopyPasta.txt in WplusJetsAnalysis): 

srun --pty -p largemem --time=0-12:00:00 bash

OR

srun --pty -c 4 --time=0-12:00:00 bash

THEN 

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh

End with:

lsetup "root 6.22.00-python3-x86_64-centos7-gcc8-opt"
(for Pythia8.302) -> CURRENTLY WORKING

2) Submit a job on the cluster or via terminal:

Via cluster:
Set output file names in MyPythia8Simul.sbatch
In terminal, type: sbatch MyPythia8Simul.sbatch

Via terminal: 
In terminal, type: ./MyPythia8Simul MyPythia8Simul_Main.cmnd -outroot name_of_root_output > name_of_txt_output.txt 
***Note: you don't need to type ".root" for name_of_root_output

3) Once this runs, the txt output will be saved under "outfiles" and the root output will be saved under /WplusJetsAnalysis/pythia-outputs

This isall for WplusJetsBaseline. To continue processing, move to WplusJetsAnalysis.


