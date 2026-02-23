# WplusJetsBaseline
Pythia baseline code developed by bjohns

Generates events using Pythia simulator and implements the detector effects.
Outputs a file with truth-level and reco-level data.

These instructions are a quick reference for running. For more information about steering, see Ben's full instructions and explanation at: https://github.com/bjohns664297/WplusJetsBaseline

### Setup Instructions<br />

1) Once logging into cluster, type the following 4 commands one after another in the terminal (or see RootCopyPasta.txt in WplusJetsAnalysis): 

srun --pty -p largemem --time=0-12:00:00 bash

OR

#THIS CURRENTLY DOESN'T WORK
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

If submit via cluster, can check job status by typing in terminal: squeue --me

3) Once this runs, the txt output will be saved under "outfiles" and the root output will be saved under /WplusJetsAnalysis/pythia-outputs

This is all for WplusJetsBaseline. To continue processing, move to WplusJetsAnalysis.

### Directions for adding new functions & variables<br />
Examples are included.<br /> 

`ANA_utils.cc`: <br />
1. FUNCTION: Add new function and short description into revelant comment block (i.e. Truth Level Particle Reconstruction Functions) <br />
2. FUNCTION: Initialize new function and revelant inputs in body of code: <br />
    void ANA_utils::Get_Tau_Info(Pythia8::Event event, std::vector<int> vecboson_index, std::vector<TruthPart>* p_Tau_Coll, std::vector<TruthPart>* p_AntiTau_Coll, std::vector<TruthPart>* p_TauDecay_Coll, std::vector<TruthPart>* p_AntiTauDecay_Coll){} <br />

`ANA_utils.h`:<br />
1. FUNCTION: Add function to ANA_utils class:<br />
    void Get_Tau_Info(Pythia8::Event event, std::vector<int> vecboson_index, std::vector<TruthPart>* p_Tau_Coll, std::vector<TruthPart>* p_AntiTau_Coll, std::vector<TruthPart>* p_TauDecay_Coll, std::vector<TruthPart>* p_AntiTauDecay_Coll)<br />


`MyPythia8Simul.cc`:<br />
1. VARIABLES: Initialize branches and pointers for new kinematic variables & multiplicities:<br />
    tree->Branch("tau_born_pt",&tau_born_pt)<br />
    tree->Branch("nTauBorn", &nTauBorn)<br />
2. FUNCTION: Set pointers for "Colls" that were defined<br />
    p_Tau_Coll = &Tau_Coll;<br />
3. FUNCTION: Add function call where relevant:<br />
     if (vecbosindex.size() > 0) <br />
      {<br />
        myUtils.Get_Tau_Info(event, vecbosindex, p_Tau_Coll, p_AntiTau_Coll, p_TauDecay_Coll, p_AntiTauDecay_Coll);<br />
      }<br />
4. FUNCTON & VARIABLES: Set multiplicity variable to zero and loop through Coll containers, saving variables info to ntuple:<br />
    nTauBorn = 0;<br />
    if (p_Tau_Coll->size() != 0)<br />
      {<br />
        for (size_t i = 0; i < p_Tau_Coll->size(); i++)<br />
          {<br />
            tau_born_pt.push_back((Tau_Coll[i]).Pt())<br />
5. VARIABLES: Clear event-based vectors:<br />
    tau_born_pt.clear()<br />

`MyPythiaSimul.h`:<br />
1. FUNCTION: Declare containers & pointers:<br />
    std::vector<TruthPart> Tau_Coll;<br /> 
    std::vector<TruthPart>* p_Tau_Coll;<br />
2. VARIABLES: Declare multiplicities and kinematic veriables:<br />
    int nTauBorn<br />
    std::vector<float> tau_born_pt<br />

`MyPythia8Simul_RunParameters.cmnd`:<br />
1. Specify number of events<br />

**Don't forget to run "make" in terminal after making changes!<br />**
