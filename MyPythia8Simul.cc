//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: MyPythia8Simul_main
//
// Purpose:  This is the main code for generating events with pythia 8, as well as
//           storing relevant output variables in Root Trees and histograms. Jets are
//           reconstructed with FastJet. The goal is to produce datasets with all the
//           relevant variables for doing theoretical predictions and sensitivity
//           studies.
//
//    Note: The main() method is used to generate events. It should not be modified
//          by users. General flags setting, histogram definitions and object
//          declarations should be made in the MyAnalysis::init(). The analysis to be
//          run in the loop must be written in void MyAnalysis::analyze(Event& event).
//          Note that the analysis here is only to calculate the physics quantities
//          to store in the Ntuple. Finally, the Pythia settings needed for a specific
//          production are set in MyPythiaSimul.cmnd.
//
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


//*************************************************************************************
// Include headers
//*************************************************************************************

#include "MyPythia8Simul.h"


using namespace Pythia8;


//*************************************************************************************






//************************************************************************************
// void Usage(const char* exename)
//
// Function to be called giving the list of all flags and key-value pairs defined. 
// Options available are:
//
//      <root_output>    = [your choice of output name] (no need to add .root, this is
//                                                  done in the code)
//
//      <hepmc_output>   = [your choice of hepmc name] (no need to add .dat, this is                                         
//                                                  done in the code)
//  
//************************************************************************************
void Usage(const char* exename) 
{
  cout << "Usage : MyPythia8Simul <inputfile.cmnd> -outroot <root_output> -outhepmc <hepmc_output>" << endl;
  cout << "Look at MyPythia8Simul.cc for valid <key> values." << endl;
}

//*************************************************************************************






//*************************************************************************************
// Initialization Code
//*************************************************************************************
void MyAnalysis::init()
{


  // Initialize counters and other global variables
  // ----------------------------------------------

  // Note: These must be define in MyPythiaSimul.h to be in-scope for all functions


  // Debug flag
  // ..........

  debug = false;


  // Number of events
  // ................

  nEvt = 0;

  // Book Event-based Histograms
  // ---------------------------


  // Book Ntuples
  // ------------

  tree = new TTree("ParticleTree","Particle Data");

  tree->Branch("NumHardJets",&NumHardJets);     // BSJ

  tree->Branch("event_weights",&event_weights);

  tree->Branch("top_pt",&top_pt);
  tree->Branch("top_eta",&top_eta);
  tree->Branch("top_phi",&top_phi);
  tree->Branch("top_E",&top_E);

  tree->Branch("neutrino_pt",&neutrino_pt);
  tree->Branch("neutrino_eta",&neutrino_eta);
  tree->Branch("neutrino_phi",&neutrino_phi);
  tree->Branch("neutrino_E",&neutrino_E);
  tree->Branch("neutrino_PdgId",&neutrino_PdgId);

  tree->Branch("muon_bare_pt",&muon_bare_pt);
  tree->Branch("muon_bare_eta",&muon_bare_eta);
  tree->Branch("muon_bare_phi",&muon_bare_phi);
  tree->Branch("muon_bare_E",&muon_bare_E);
  tree->Branch("muon_bare_charge",&muon_bare_charge);

  tree->Branch("muon_dress_pt",&muon_dress_pt);
  tree->Branch("muon_dress_eta",&muon_dress_eta);
  tree->Branch("muon_dress_phi",&muon_dress_phi);
  tree->Branch("muon_dress_E",&muon_dress_E);
  tree->Branch("muon_dress_charge",&muon_dress_charge);

  tree->Branch("muon_born_pt",&muon_born_pt);     // BSJ
  tree->Branch("muon_born_eta",&muon_born_eta);
  tree->Branch("muon_born_phi",&muon_born_phi);
  tree->Branch("muon_born_E",&muon_born_E);
  tree->Branch("muon_born_charge",&muon_born_charge);

  tree->Branch("electron_bare_pt",&electron_bare_pt);
  tree->Branch("electron_bare_eta",&electron_bare_eta);
  tree->Branch("electron_bare_phi",&electron_bare_phi);
  tree->Branch("electron_bare_E",&electron_bare_E);
  tree->Branch("electron_bare_charge",&electron_bare_charge);

  tree->Branch("electron_dress_pt",&electron_dress_pt);
  tree->Branch("electron_dress_eta",&electron_dress_eta);
  tree->Branch("electron_dress_phi",&electron_dress_phi);
  tree->Branch("electron_dress_E",&electron_dress_E);
  tree->Branch("electron_dress_charge",&electron_dress_charge);

  tree->Branch("electron_born_pt",&electron_born_pt);   // BSJ
  tree->Branch("electron_born_eta",&electron_born_eta);
  tree->Branch("electron_born_phi",&electron_born_phi);
  tree->Branch("electron_born_E",&electron_born_E);
  tree->Branch("electron_born_charge",&electron_born_charge);

  tree->Branch("boson_pt",&boson_pt);
  tree->Branch("boson_eta",&boson_eta);
  tree->Branch("boson_phi",&boson_phi);
  tree->Branch("boson_E",&boson_E);
  tree->Branch("boson_ID",&boson_ID);

  tree->Branch("promptphoton_pt",&promptphoton_pt);       // BSJ
  tree->Branch("promptphoton_eta",&promptphoton_eta);
  tree->Branch("promptphoton_phi",&promptphoton_phi);
  tree->Branch("promptphoton_E",&promptphoton_E);

  tree->Branch("lightjet_bare_pt",&lightjet_bare_pt);
  tree->Branch("lightjet_bare_eta",&lightjet_bare_eta);
  tree->Branch("lightjet_bare_phi",&lightjet_bare_phi);
  tree->Branch("lightjet_bare_E",&lightjet_bare_E);
  tree->Branch("lightjet_bare_nPart",&lightjet_bare_nPart);

  tree->Branch("jet_dress_pt",&jet_dress_pt);
  tree->Branch("jet_dress_eta",&jet_dress_eta);
  tree->Branch("jet_dress_phi",&jet_dress_phi);
  tree->Branch("jet_dress_E",&jet_dress_E);
  tree->Branch("jet_dress_nPart",&jet_dress_nPart);

  tree->Branch("bjet_bare_pt",&bjet_bare_pt);
  tree->Branch("bjet_bare_eta",&bjet_bare_eta);
  tree->Branch("bjet_bare_phi",&bjet_bare_phi);
  tree->Branch("bjet_bare_E",&bjet_bare_E);
  tree->Branch("bjet_bare_nPart",&bjet_bare_nPart);

  tree->Branch("jet_born_pt",&jet_born_pt);
  tree->Branch("jet_born_eta",&jet_born_eta);
  tree->Branch("jet_born_phi",&jet_born_phi);
  tree->Branch("jet_born_E",&jet_born_E);
  tree->Branch("jet_born_nPart",&jet_born_nPart);

  tree->Branch("largeRjet_bare_pt", &largeRjet_bare_pt);
  tree->Branch("largeRjet_bare_eta",&largeRjet_bare_eta);
  tree->Branch("largeRjet_bare_phi",&largeRjet_bare_phi);
  tree->Branch("largeRjet_bare_E",  &largeRjet_bare_E);

  tree->Branch("largeRjet_dress_pt", &largeRjet_dress_pt);
  tree->Branch("largeRjet_dress_eta",&largeRjet_dress_eta);
  tree->Branch("largeRjet_dress_phi",&largeRjet_dress_phi);
  tree->Branch("largeRjet_dress_E",  &largeRjet_dress_E);

  tree->Branch("largeRjet_born_pt", &largeRjet_born_pt);
  tree->Branch("largeRjet_born_eta",&largeRjet_born_eta);
  tree->Branch("largeRjet_born_phi",&largeRjet_born_phi);
  tree->Branch("largeRjet_born_E",  &largeRjet_born_E);

  tree->Branch("lightpartonjet_pt",&lightpartonjet_pt);
  tree->Branch("lightpartonjet_eta",&lightpartonjet_eta);
  tree->Branch("lightpartonjet_phi",&lightpartonjet_phi);
  tree->Branch("lightpartonjet_E",&lightpartonjet_E);

  tree->Branch("bpartonjet_pt",&bpartonjet_pt);
  tree->Branch("bpartonjet_eta",&bpartonjet_eta);
  tree->Branch("bpartonjet_phi",&bpartonjet_phi);
  tree->Branch("bpartonjet_E",&bpartonjet_E);

  if(BareKtSplittingScales)
    {
      tree->Branch("BarekTSplittingScale1_R04",&BarekTSplittingScale1_R04);
      tree->Branch("BarekTSplittingScale2_R04",&BarekTSplittingScale2_R04);
      tree->Branch("BarekTSplittingScale3_R04",&BarekTSplittingScale3_R04);
      tree->Branch("BarekTSplittingScale1_R10",&BarekTSplittingScale1_R10);
      tree->Branch("BarekTSplittingScale2_R10",&BarekTSplittingScale2_R10);
      tree->Branch("BarekTSplittingScale3_R10",&BarekTSplittingScale3_R10);
    }

  if(DressKtSplittingScales)
    {
      tree->Branch("DresskTSplittingScale1_R04",&DresskTSplittingScale1_R04);
      tree->Branch("DresskTSplittingScale2_R04",&DresskTSplittingScale2_R04);
      tree->Branch("DresskTSplittingScale3_R04",&DresskTSplittingScale3_R04);
      tree->Branch("DresskTSplittingScale1_R10",&DresskTSplittingScale1_R10);
      tree->Branch("DresskTSplittingScale2_R10",&DresskTSplittingScale2_R10);
      tree->Branch("DresskTSplittingScale3_R10",&DresskTSplittingScale3_R10);
    }

  if(BornKtSplittingScales)
    {
      tree->Branch("BornkTSplittingScale1_R04",&BornkTSplittingScale1_R04);
      tree->Branch("BornkTSplittingScale2_R04",&BornkTSplittingScale2_R04);
      tree->Branch("BornkTSplittingScale3_R04",&BornkTSplittingScale3_R04);
      tree->Branch("BornkTSplittingScale1_R10",&BornkTSplittingScale1_R10);
      tree->Branch("BornkTSplittingScale2_R10",&BornkTSplittingScale2_R10);
      tree->Branch("BornkTSplittingScale3_R10",&BornkTSplittingScale3_R10);
    }

  tree->Branch("nTop",&nTop);
  tree->Branch("nNeutrino",&nNeutrino);
  tree->Branch("nMuonBare",&nMuonBare);
  tree->Branch("nElectronBare",&nElectronBare);
  tree->Branch("nMuonDress",&nMuonDress);
  tree->Branch("nElectronDress",&nElectronDress);
  tree->Branch("nMuonBorn",&nMuonBorn);
  tree->Branch("nElectronBorn",&nElectronBorn);
  tree->Branch("nLightjetBare",&nLightjetBare);
  tree->Branch("nJetDress",&nJetDress);
  tree->Branch("nBjetBare",&nBjetBare);
  tree->Branch("nJetBorn",&nJetBorn);
  tree->Branch("nLargeRjetBare",&nLargeRjetBare);
  tree->Branch("nLargeRjetdress",&nLargeRjetDress);
  tree->Branch("nLargeRjetBorn",&nLargeRjetBorn);
  tree->Branch("nLightpartonjet",&nLightpartonjet);
  tree->Branch("nBpartonjet",&nBpartonjet);
  tree->Branch("nBoson",&nBoson);
  tree->Branch("nPromptPhoton",&nPromptPhotons);
  tree->Branch("Met",&Met);
  tree->Branch("Met_phi",&Met_phi);

  tree->Branch("glob_TransvSphericity",&glob_TransSpher);
  tree->Branch("glob_TransvThrustMajor",&glob_TransThrustMaj);
  tree->Branch("glob_TransvThrustMinor",&glob_TransThrustMin);
  tree->Branch("glob_TransvThrustMajorWithResidue",&glob_TransThrustMajRes);
  tree->Branch("glob_TransvThrustMinorWithResidue",&glob_TransThrustMinRes);
  tree->Branch("glob_TransvThrustMajorWithSuppress",&glob_TransThrustMajSup);
  tree->Branch("glob_TransvThrustMinorWithSuppress",&glob_TransThrustMinSup);
  tree->Branch("glob_SumMassWithResidue",&glob_SumMassRes);
  tree->Branch("glob_HeavyMassWithResidue",&glob_HeavyMassRes);
  tree->Branch("glob_SumMassWithSuppress",&glob_SumMassSup);
  tree->Branch("glob_HeavyMassWithSuppress",&glob_HeavyMassSup);
  tree->Branch("glob_TotalBroadeningsWithResidue",&glob_TotBroadRes);
  tree->Branch("glob_WideBroadeningsWithResidue",&glob_WideBroadRes);
  tree->Branch("glob_TotalBroadeningsWithSuppress",&glob_TotBroadSup);
  tree->Branch("glob_WideBroadeningsWithSuppress",&glob_WideBroadSup);
  tree->Branch("glob_SuperSpherocity",&glob_SuperSphero);


}

//*************************************************************************************






//*************************************************************************************
// Analysis Code
//*************************************************************************************
void MyAnalysis::analyze(Event& event, Event& partonevent, std::vector<double> EventWeights, int ihardjet = -999)
  {

    // Declare an Analysis Utilities class object
    // ------------------------------------------

    // Note: To be able to access the functions define there

    ANA_utils myUtils;


    // Fill truth particle and jets information
    // ----------------------------------------

    // Set pointers
    // ............

    // Note: In each case, we first need to set a pointer to the vector of TruthPart containing the relevant particles,
    //       and then we call the function using the pointer.

    p_Top_Coll = &Top_Coll;
    p_Vecboson_Coll = &Vecboson_Coll;
    p_Promptphoton_Coll = &Promptphoton_Coll;   // BSJ
    p_PromptLeptonBare_Coll = &PromptLeptonBare_Coll;
    p_LeptonConversion_Coll = &LeptonConversion_Coll;
    p_LeptonDress_Coll = &LeptonDress_Coll;
    p_LeptonBorn_Coll = &LeptonBorn_Coll;
    p_FSRPhoton_Coll = &FSRPhoton_Coll;
    p_DressPhoton_Coll = &DressPhoton_Coll;
    p_Neutrino_Coll = &Neutrino_Coll;
    p_TruthBareSmallRJets_Coll = &TruthBareSmallRJets_Coll;
    p_TruthDressSmallRJets_Coll = &TruthDressSmallRJets_Coll;
    p_TruthBornSmallRJets_Coll = &TruthBornSmallRJets_Coll;
    p_TruthBareLargeRJets_Coll = &TruthBareLargeRJets_Coll;
    p_TruthDressLargeRJets_Coll = &TruthDressLargeRJets_Coll;
    p_TruthBornLargeRJets_Coll = &TruthBornLargeRJets_Coll;
    p_PartonJets_Coll = &PartonJets_Coll;           // Small-R parton jets
    p_BareKtSplittingScale1_R04 = &BareKtSplittingScale1_R04;
    p_BareKtSplittingScale2_R04 = &BareKtSplittingScale2_R04;
    p_BareKtSplittingScale3_R04 = &BareKtSplittingScale3_R04;
    p_BareKtSplittingScale1_R10 = &BareKtSplittingScale1_R10;
    p_BareKtSplittingScale2_R10 = &BareKtSplittingScale2_R10;
    p_BareKtSplittingScale3_R10 = &BareKtSplittingScale3_R10;
    p_DressKtSplittingScale1_R04 = &DressKtSplittingScale1_R04;
    p_DressKtSplittingScale2_R04 = &DressKtSplittingScale2_R04;
    p_DressKtSplittingScale3_R04 = &DressKtSplittingScale3_R04;
    p_DressKtSplittingScale1_R10 = &DressKtSplittingScale1_R10;
    p_DressKtSplittingScale2_R10 = &DressKtSplittingScale2_R10;
    p_DressKtSplittingScale3_R10 = &DressKtSplittingScale3_R10;
    p_BornKtSplittingScale1_R04 = &BornKtSplittingScale1_R04;
    p_BornKtSplittingScale2_R04 = &BornKtSplittingScale2_R04;
    p_BornKtSplittingScale3_R04 = &BornKtSplittingScale3_R04;
    p_BornKtSplittingScale1_R10 = &BornKtSplittingScale1_R10;
    p_BornKtSplittingScale2_R10 = &BornKtSplittingScale2_R10;
    p_BornKtSplittingScale3_R10 = &BornKtSplittingScale3_R10;

  
    // Get Tops
    // ........

    myUtils.Get_Tops(event, p_Top_Coll);

  
    // Get Vector bosons
    // .................

    myUtils.Get_VectorBosons(event, p_Vecboson_Coll);



    // Get Prompt photons   (BSJ)
    // ..................

    myUtils.Get_PromptPhotons(event, p_Promptphoton_Coll);


    // Bare, Dressed, and Born leptons
    // ...............................

    // Note: Need to find the indices of the vector bosons found above. Don't call the function if none are found.

    std::vector<int> vecbosindex;
    for (int vb_i = 0; vb_i < Vecboson_Coll.size(); vb_i++)
      {
        vecbosindex.push_back( (Vecboson_Coll[vb_i]).Index() );
      }

    if (vecbosindex.size() > 0) 
      {
        myUtils.Get_BarePromptLepton(event, vecbosindex, p_PromptLeptonBare_Coll, p_LeptonConversion_Coll, p_FSRPhoton_Coll, p_Neutrino_Coll);
        myUtils.Get_BornPromptLepton(event, vecbosindex, p_LeptonBorn_Coll);
      }

    if (PromptLeptonBare_Coll.size()!=0) myUtils.Get_DressPromptLepton(PromptLeptonBare_Coll,FSRPhoton_Coll, p_LeptonDress_Coll, p_DressPhoton_Coll);

				   

    // True Jets
    // .........

    // Note 1: A list of stable particles not to be clustered in jets must first be defined for truth jets.

    // Note 2: A different function is called for final state particle jets and for pre-hadronization parton jets.

    // Note 3: The Kt Splitting scale flags are defined in the header file.

    // Note 4: For bare jets, we only exclude the bare prompt leptons. For dressed jet, we also exlcude dressing FSR photons. For Born jets,
    //         we exclude prompt photon, all fsr photons, and conversion leptons. In principle conversion leptons should be removed from
    //         dressed jets, but they are not in ATLAS, so they are left here too. 


    // Remove Prompt leptons from any jet definition
    std::vector<int> skippart_born;   // born = fsr
    std::vector<int> skippart_bare;   // bare = blank (no suffix)
    std::vector<int> skippart_dress;
    for (int i_part = 0; i_part < PromptLeptonBare_Coll.size(); i_part++) 
      {
        skippart_born.push_back((PromptLeptonBare_Coll[i_part]).Index());
        skippart_bare.push_back((PromptLeptonBare_Coll[i_part]).Index());
        skippart_dress.push_back((PromptLeptonBare_Coll[i_part]).Index());
      }

    // BSJ vv
    for (int i_phot = 0; i_phot < Promptphoton_Coll.size(); i_phot++) 
      {
        skippart_born.push_back((Promptphoton_Coll[i_phot]).Index());
        skippart_bare.push_back((Promptphoton_Coll[i_phot]).Index());
        skippart_dress.push_back((Promptphoton_Coll[i_phot]).Index());
      
      }

  
    // Particle bare small-R jets
    float jetcut = 30.; // jet pT cut
    myUtils.TrueJetsReco(event, skippart_bare, p_TruthBareSmallRJets_Coll, jetcut, BareKtSplittingScales, p_BareKtSplittingScale1_R04, p_BareKtSplittingScale2_R04, p_BareKtSplittingScale3_R04, 0.4);

  
    // Also remove dressing FSR photons from dressed jet definition
    for (int i_dress = 0; i_dress < DressPhoton_Coll.size(); i_dress++)
    {
      skippart_dress.push_back((DressPhoton_Coll[i_dress]).Index());
    }

  
    // Particle dress small-R jets
    myUtils.TrueJetsReco(event, skippart_dress, p_TruthDressSmallRJets_Coll, jetcut, DressKtSplittingScales, p_DressKtSplittingScale1_R04, p_DressKtSplittingScale2_R04, p_DressKtSplittingScale3_R04, 0.4);


  
    // Remove any photons and conversion leptons for Born jets
    for (int i_fsr = 0; i_fsr < FSRPhoton_Coll.size(); i_fsr++)
      {
        skippart_born.push_back((FSRPhoton_Coll[i_fsr]).Index());
      }

    for (int i_part2 = 0; i_part2 < LeptonConversion_Coll.size(); i_part2++) 
      {
        skippart_born.push_back((LeptonConversion_Coll[i_part2]).Index());
      }

    // Particle born small-R jets
    myUtils.TrueJetsReco(event, skippart_born, p_TruthBornSmallRJets_Coll, jetcut, BornKtSplittingScales, p_BornKtSplittingScale1_R04, p_BornKtSplittingScale2_R04, p_BornKtSplittingScale3_R04, 0.4);

    // Particle bare large-R jets
    myUtils.TrueJetsReco(event, skippart_bare, p_TruthBareLargeRJets_Coll, jetcut, BareKtSplittingScales, p_BareKtSplittingScale1_R10, p_BareKtSplittingScale2_R10, p_BareKtSplittingScale3_R10, 1.0);
  
    // Particle dress large-R jets
    myUtils.TrueJetsReco(event, skippart_dress, p_TruthDressLargeRJets_Coll, jetcut, DressKtSplittingScales, p_DressKtSplittingScale1_R10, p_DressKtSplittingScale2_R10, p_DressKtSplittingScale3_R10, 1.0);
  
    // Particle born large-R jets
    myUtils.TrueJetsReco(event, skippart_born, p_TruthBornLargeRJets_Coll, jetcut, BornKtSplittingScales, p_BornKtSplittingScale1_R10, p_BornKtSplittingScale2_R10, p_BornKtSplittingScale3_R10, 1.0);

    // Parton small-R jets
    myUtils.PartonJetsReco(event, partonevent, p_PartonJets_Coll, jetcut);

  
    // Fill ntuples and histograms
    // ---------------------------

    // Event weights
    // .............

    event_weights = EventWeights;


    // vv BSJ -> not necessary? -> Pretty sure this is only used for CKKW merge 
    // Hard Jet multiplicity
    // .....................

    // Note: Is not -999 only if obtained from lhe files merged

    NumHardJets = ihardjet;




    // W and Z
    // .......

    nBoson = p_Vecboson_Coll->size();
    if (nBoson != 0)
      {
        for (size_t i = 0; i < p_Vecboson_Coll->size(); i++)
          {
            boson_pt.push_back((Vecboson_Coll[i]).Pt());
            boson_eta.push_back((Vecboson_Coll[i]).Eta());
            double i_phi = (Vecboson_Coll[i]).Phi();
            if (i_phi < 0.) i_phi = (Vecboson_Coll[i]).Phi() + 6.283185307;
            boson_phi.push_back(i_phi);
            boson_E.push_back((Vecboson_Coll[i]).E());
            boson_ID.push_back((Vecboson_Coll[i]).Pdgid());
          }
      }


    // Prompt Photons
    // ..............
  
    nPromptPhotons = p_Promptphoton_Coll->size();
    if (nPromptPhotons != 0)
      {
        for (size_t i = 0; i < p_Promptphoton_Coll->size(); i++)
          {
            promptphoton_pt.push_back((Promptphoton_Coll[i]).Pt());
            promptphoton_eta.push_back((Promptphoton_Coll[i]).Eta());
            double i_phi = (Promptphoton_Coll[i]).Phi();
            if (i_phi < 0.) i_phi = (Promptphoton_Coll[i]).Phi() + 6.283185307;
            promptphoton_phi.push_back(i_phi);
            promptphoton_E.push_back((Promptphoton_Coll[i]).E());
          }
      }


    // Top quarks
    // ..........

    nTop = p_Top_Coll->size();
    if (nTop != 0)
      {
        for (size_t i = 0; i < p_Top_Coll->size(); i++)
          {
            top_pt.push_back((Top_Coll[i]).Pt());
            top_eta.push_back((Top_Coll[i]).Eta());
            double i_phi = (Top_Coll[i]).Phi();
            if (i_phi < 0.) i_phi = (Top_Coll[i]).Phi() + 6.283185307;	      
            top_phi.push_back(i_phi);
            top_E.push_back((Top_Coll[i]).E());
          }
      }


    // Neutrinos and Met
    // .................

    nNeutrino = p_Neutrino_Coll->size();
    if (nNeutrino != 0)
      {
        for (size_t i = 0; i < p_Neutrino_Coll->size(); i++)
          {
            neutrino_pt.push_back((Neutrino_Coll[i]).Pt());
            neutrino_eta.push_back((Neutrino_Coll[i]).Eta());
            double i_phi = (Neutrino_Coll[i]).Phi();
            if (i_phi < 0.) i_phi = (Neutrino_Coll[i]).Phi() + 6.283185307;	      
            neutrino_phi.push_back(i_phi);
            neutrino_E.push_back((Neutrino_Coll[i]).E());
            neutrino_PdgId.push_back((Neutrino_Coll[i]).Pdgid());
          }
      }

    float sqsum;
    float sum_y=0;
    float sum_x=0;

    for (auto i : Neutrino_Coll)
    //for (std::vector<double>::iterator i : Neutrino_Coll)
      {
        sum_y += i.Py();
        sum_x += i.Px();
      }

    sqsum = pow(sum_x, 2.0) + pow(sum_y, 2.0);
    Met = sqrt(sqsum);
    Met_phi = atan2(sum_y,sum_x);
    if (Met_phi < 0.) Met_phi = Met_phi + 6.283185307;

  




    // light and b-jets for particle and parton jets
    // .............................................

    // Note: For b-jets, we use the b-quark tag for parton jets, and the b-hadron tag for particle jets.

    // Note: While we keep light and b-jets separate for bare jets, for Born jets we keep them in the same collection

  
    // Bare Small-R light and b-jets for particle jets

    nLightjetBare  = 0;
    nBjetBare      = 0;
    if (p_TruthBareSmallRJets_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_TruthBareSmallRJets_Coll->size(); i++)
          {
            if ((TruthBareSmallRJets_Coll[i]).BHTag())
              {
                bjet_bare_pt.push_back((TruthBareSmallRJets_Coll[i]).Pt());
                bjet_bare_eta.push_back((TruthBareSmallRJets_Coll[i]).Eta());
                bjet_bare_phi.push_back((TruthBareSmallRJets_Coll[i]).Phi());
                bjet_bare_E.push_back((TruthBareSmallRJets_Coll[i]).E());
                bjet_bare_nPart.push_back((TruthBareSmallRJets_Coll[i]).Npart());
                nBjetBare += 1;
              }
            else
              {
                lightjet_bare_pt.push_back((TruthBareSmallRJets_Coll[i]).Pt());
                lightjet_bare_eta.push_back((TruthBareSmallRJets_Coll[i]).Eta());
                lightjet_bare_phi.push_back((TruthBareSmallRJets_Coll[i]).Phi());
                lightjet_bare_E.push_back((TruthBareSmallRJets_Coll[i]).E());
                lightjet_bare_nPart.push_back((TruthBareSmallRJets_Coll[i]).Npart());
                nLightjetBare += 1;
              }
          }
      }

    // Bare large-R jets
    nLargeRjetBare = 0;
    if (p_TruthBareLargeRJets_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_TruthBareLargeRJets_Coll->size(); i++)
          {
            largeRjet_bare_pt.push_back((TruthBareLargeRJets_Coll[i]).Pt());
            largeRjet_bare_eta.push_back((TruthBareLargeRJets_Coll[i]).Eta());
            largeRjet_bare_phi.push_back((TruthBareLargeRJets_Coll[i]).Phi());
            largeRjet_bare_E.push_back((TruthBareLargeRJets_Coll[i]).E());
            nLargeRjetBare += 1;
          }
      }

  
    // Dress small-R jets (include both light and B-jets)
    nJetDress = 0;
    if (p_TruthDressSmallRJets_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_TruthDressSmallRJets_Coll->size(); i++)
          {
          jet_dress_pt.push_back((TruthDressSmallRJets_Coll[i]).Pt());
          jet_dress_eta.push_back((TruthDressSmallRJets_Coll[i]).Eta());
          jet_dress_phi.push_back((TruthDressSmallRJets_Coll[i]).Phi());
          jet_dress_E.push_back((TruthDressSmallRJets_Coll[i]).E());
          jet_dress_nPart.push_back((TruthDressSmallRJets_Coll[i]).Npart());
          nJetDress += 1;
          }
      }

  
    // Dress large-R jets
    nLargeRjetDress = 0;
    if (p_TruthDressLargeRJets_Coll->size() != 0)
    {
      for (size_t i = 0; i < p_TruthDressLargeRJets_Coll->size(); i++)
        {
        largeRjet_dress_pt.push_back((TruthDressLargeRJets_Coll[i]).Pt());
        largeRjet_dress_eta.push_back((TruthDressLargeRJets_Coll[i]).Eta());
        largeRjet_dress_phi.push_back((TruthDressLargeRJets_Coll[i]).Phi());
        largeRjet_dress_E.push_back((TruthDressLargeRJets_Coll[i]).E());
        nLargeRjetDress += 1;
        }
    }

  

    // Born small-R jets (include both light and B-jets)
    nJetBorn = 0;
    if (p_TruthBornSmallRJets_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_TruthBornSmallRJets_Coll->size(); i++)
          {
	          jet_born_pt.push_back((TruthBornSmallRJets_Coll[i]).Pt());
	          jet_born_eta.push_back((TruthBornSmallRJets_Coll[i]).Eta());
	          jet_born_phi.push_back((TruthBornSmallRJets_Coll[i]).Phi());
	          jet_born_E.push_back((TruthBornSmallRJets_Coll[i]).E());
	          jet_born_nPart.push_back((TruthBornSmallRJets_Coll[i]).Npart());
	          nJetBorn += 1;
          }
      }

  
    // Born large-R jets
    nLargeRjetBorn = 0;
    if (p_TruthBornLargeRJets_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_TruthBornLargeRJets_Coll->size(); i++)
          {
            largeRjet_born_pt.push_back((TruthBornLargeRJets_Coll[i]).Pt());
            largeRjet_born_eta.push_back((TruthBornLargeRJets_Coll[i]).Eta());
            largeRjet_born_phi.push_back((TruthBornLargeRJets_Coll[i]).Phi());
            largeRjet_born_E.push_back((TruthBornLargeRJets_Coll[i]).E());
            nLargeRjetBorn += 1;
          }
      }

  

    // Small-R parton light and b-jets  

    nLightpartonjet = 0;
    nBpartonjet = 0;
    if (p_PartonJets_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_PartonJets_Coll->size(); i++)
          {
            if ((PartonJets_Coll[i]).BQTag())
              {
                bpartonjet_pt.push_back((PartonJets_Coll[i]).Pt());
                bpartonjet_eta.push_back((PartonJets_Coll[i]).Eta());
                bpartonjet_phi.push_back((PartonJets_Coll[i]).Phi());
                bpartonjet_E.push_back((PartonJets_Coll[i]).E());
                nBpartonjet += 1;
              }
            else
              {
                lightpartonjet_pt.push_back((PartonJets_Coll[i]).Pt());
                lightpartonjet_eta.push_back((PartonJets_Coll[i]).Eta());
                lightpartonjet_phi.push_back((PartonJets_Coll[i]).Phi());
                lightpartonjet_E.push_back((PartonJets_Coll[i]).E());
                nLightpartonjet += 1;
              }
          }
      }


  
    // kT Splitting scales

    // Bare
    int nSmallRjetBare = nLightjetBare + nBjetBare;
    if(p_BareKtSplittingScale1_R04->size() != 0)
      {
        for (size_t i = 0; i < p_BareKtSplittingScale1_R04->size(); i++)
          {
            BarekTSplittingScale1_R04.push_back(BareKtSplittingScale1_R04[i]);
          }
      }
    if(p_BareKtSplittingScale2_R04->size() != 0)
      {
        for (size_t i = 0; i < p_BareKtSplittingScale2_R04->size(); i++)
          {
            BarekTSplittingScale2_R04.push_back(BareKtSplittingScale2_R04[i]);
          }
      }
    if(p_BareKtSplittingScale3_R04->size() != 0)
      {
        for (size_t i = 0; i < p_BareKtSplittingScale3_R04->size(); i++)
          {
            BarekTSplittingScale3_R04.push_back(BareKtSplittingScale3_R04[i]);
          }
      }
  
    if(p_BareKtSplittingScale1_R10->size() != 0)
      {
        for (size_t i = 0; i < p_BareKtSplittingScale1_R10->size(); i++)
          {
            BarekTSplittingScale1_R10.push_back(BareKtSplittingScale1_R10[i]);
          }
      }
    if(p_BareKtSplittingScale2_R10->size() != 0)
      {
        for (size_t i = 0; i < p_BareKtSplittingScale2_R10->size(); i++)
          {
            BarekTSplittingScale2_R10.push_back(BareKtSplittingScale2_R10[i]);
          }
      }

    if(p_BareKtSplittingScale3_R10->size() != 0)
      {
        for (size_t i = 0; i < p_BareKtSplittingScale3_R10->size(); i++)
          {
            BarekTSplittingScale3_R10.push_back(BareKtSplittingScale3_R10[i]);
          }
      }

  
    // Dress
    int nSmallRjetDress = nJetDress;
    if(p_DressKtSplittingScale1_R04->size() != 0)
      {
        for (size_t i = 0; i < p_DressKtSplittingScale1_R04->size(); i++)
          {
            DresskTSplittingScale1_R04.push_back(DressKtSplittingScale1_R04[i]);
          }
      }
    if(p_DressKtSplittingScale2_R04->size() != 0)
      {
        for (size_t i = 0; i < p_DressKtSplittingScale2_R04->size(); i++)
          {
            DresskTSplittingScale2_R04.push_back(DressKtSplittingScale2_R04[i]);
          }
      }
    if(p_DressKtSplittingScale3_R04->size() != 0)
      {
        for (size_t i = 0; i < p_DressKtSplittingScale3_R04->size(); i++)
          {
            DresskTSplittingScale3_R04.push_back(DressKtSplittingScale3_R04[i]);
          }
      }
  
    if(p_DressKtSplittingScale1_R10->size() != 0)
      {
        for (size_t i = 0; i < p_DressKtSplittingScale1_R10->size(); i++)
          {
            DresskTSplittingScale1_R10.push_back(DressKtSplittingScale1_R10[i]);
          }
      }
    if(p_DressKtSplittingScale2_R10->size() != 0)
      {
        for (size_t i = 0; i < p_DressKtSplittingScale2_R10->size(); i++)
          {
            DresskTSplittingScale2_R10.push_back(DressKtSplittingScale2_R10[i]);
          }
      }

    if(p_DressKtSplittingScale3_R10->size() != 0)
      {
        for (size_t i = 0; i < p_DressKtSplittingScale3_R10->size(); i++)
          {
            DresskTSplittingScale3_R10.push_back(DressKtSplittingScale3_R10[i]);
          }
      }

  
    // Born
    if(p_BornKtSplittingScale1_R04->size() != 0)
      {
        for (size_t i = 0; i < p_BornKtSplittingScale1_R04->size(); i++)
          {
            BornkTSplittingScale1_R04.push_back(BornKtSplittingScale1_R04[i]);
          }
      }
  
    if(p_BornKtSplittingScale2_R04->size() != 0)
      {
        for (size_t i = 0; i < p_BornKtSplittingScale2_R04->size(); i++)
          {
            BornkTSplittingScale2_R04.push_back(BornKtSplittingScale2_R04[i]);
          }
      }
  
    if(p_BornKtSplittingScale3_R04->size() != 0)
      {
        for (size_t i = 0; i < p_BornKtSplittingScale3_R04->size(); i++)
          {
            BornkTSplittingScale3_R04.push_back(BornKtSplittingScale3_R04[i]);
          }
      }

    if(p_BornKtSplittingScale1_R10->size() != 0)
      {
        for (size_t i = 0; i < p_BornKtSplittingScale1_R10->size(); i++)
          {
            BornkTSplittingScale1_R10.push_back(BornKtSplittingScale1_R10[i]);
          }
      }
  
    if(p_BornKtSplittingScale2_R10->size() != 0)
      {
        for (size_t i = 0; i < p_BornKtSplittingScale2_R10->size(); i++)
          {
            BornkTSplittingScale2_R10.push_back(BornKtSplittingScale2_R10[i]);
          }
      }
  
    if(p_BornKtSplittingScale3_R10->size() != 0)
      {
        for (size_t i = 0; i < p_BornKtSplittingScale3_R10->size(); i++)
          {
            BornkTSplittingScale3_R10.push_back(BornKtSplittingScale3_R10[i]);
          }
      }

    
    // Electrons and muons
    // ...................

    // Bare
    nElectronBare = 0;
    nMuonBare     = 0;
    if (p_PromptLeptonBare_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_PromptLeptonBare_Coll->size(); i++)
          {
            if ((PromptLeptonBare_Coll[i]).Pdgid() == 11 || (PromptLeptonBare_Coll[i]).Pdgid() == -11)
              {
                electron_bare_pt.push_back((PromptLeptonBare_Coll[i]).Pt());
                electron_bare_eta.push_back((PromptLeptonBare_Coll[i]).Eta());
                double i_phi = (PromptLeptonBare_Coll[i]).Phi();
                if (i_phi < 0.) i_phi = (PromptLeptonBare_Coll[i]).Phi() + 6.283185307;	      
                electron_bare_phi.push_back(i_phi);
                electron_bare_E.push_back((PromptLeptonBare_Coll[i]).E());
                electron_bare_charge.push_back((PromptLeptonBare_Coll[i]).Charge());
                nElectronBare += 1;
              }
            if ((PromptLeptonBare_Coll[i]).Pdgid() == 13 || (PromptLeptonBare_Coll[i]).Pdgid() == -13)
              {
                muon_bare_pt.push_back((PromptLeptonBare_Coll[i]).Pt());
                muon_bare_eta.push_back((PromptLeptonBare_Coll[i]).Eta());
                double i_phi = (PromptLeptonBare_Coll[i]).Phi();
                if (i_phi < 0.) i_phi = (PromptLeptonBare_Coll[i]).Phi() + 6.283185307;	      
                muon_bare_phi.push_back(i_phi);	      
                muon_bare_E.push_back((PromptLeptonBare_Coll[i]).E());
                muon_bare_charge.push_back((PromptLeptonBare_Coll[i]).Charge());
                nMuonBare += 1;
              }
          }
      }

    // Dress
    nElectronDress = 0;
    nMuonDress     = 0;
    if (p_LeptonDress_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_LeptonDress_Coll->size(); i++)
          {
            if ((LeptonDress_Coll[i]).Pdgid() == 11 || (LeptonDress_Coll[i]).Pdgid() == -11)
              {
                electron_dress_pt.push_back((LeptonDress_Coll[i]).Pt());
                electron_dress_eta.push_back((LeptonDress_Coll[i]).Eta());
                double i_phi = (LeptonDress_Coll[i]).Phi();
                if (i_phi < 0.) i_phi = (LeptonDress_Coll[i]).Phi() + 6.283185307;	      
                electron_dress_phi.push_back(i_phi);
                electron_dress_E.push_back((LeptonDress_Coll[i]).E());
                electron_dress_charge.push_back((LeptonDress_Coll[i]).Charge());
                nElectronDress += 1;
              }
            if ((LeptonDress_Coll[i]).Pdgid() == 13 || (LeptonDress_Coll[i]).Pdgid() == -13)
              {
                muon_dress_pt.push_back((LeptonDress_Coll[i]).Pt());
                muon_dress_eta.push_back((LeptonDress_Coll[i]).Eta());
                double i_phi = (LeptonDress_Coll[i]).Phi();
                if (i_phi < 0.) i_phi = (LeptonDress_Coll[i]).Phi() + 6.283185307;	      
                muon_dress_phi.push_back(i_phi);	      
                muon_dress_E.push_back((LeptonDress_Coll[i]).E());
                muon_dress_charge.push_back((LeptonDress_Coll[i]).Charge());
                nMuonDress += 1;
              }
          }
      }

    // Born
    nElectronBorn = 0;
    nMuonBorn = 0;
    if (p_LeptonBorn_Coll->size() != 0)
      {
        for (size_t i = 0; i < p_LeptonBorn_Coll->size(); i++)
          {
            if ((LeptonBorn_Coll[i]).Pdgid() == 11 || (LeptonBorn_Coll[i]).Pdgid() == -11)
              {
                electron_born_pt.push_back((LeptonBorn_Coll[i]).Pt());
                electron_born_eta.push_back((LeptonBorn_Coll[i]).Eta());
                double i_phi = (LeptonBorn_Coll[i]).Phi();
                if (i_phi < 0.) i_phi = (LeptonBorn_Coll[i]).Phi() + 6.283185307;	      
                electron_born_phi.push_back(i_phi);
                electron_born_E.push_back((LeptonBorn_Coll[i]).E());
                electron_born_charge.push_back((LeptonBorn_Coll[i]).Charge());
                nElectronBorn += 1;
              }
            if ((LeptonBorn_Coll[i]).Pdgid() == 13 || (LeptonBorn_Coll[i]).Pdgid() == -13)
              {
                muon_born_pt.push_back((LeptonBorn_Coll[i]).Pt());
                muon_born_eta.push_back((LeptonBorn_Coll[i]).Eta());
                double i_phi = (LeptonBorn_Coll[i]).Phi();
                if (i_phi < 0.) i_phi = (LeptonBorn_Coll[i]).Phi() + 6.283185307;	      
                muon_born_phi.push_back(i_phi);	      
                muon_born_E.push_back((LeptonBorn_Coll[i]).E());
                muon_born_charge.push_back((LeptonBorn_Coll[i]).Charge());
                nMuonBorn += 1;
              }
          }
      }
  


    // Global variables
    // ................

    std::vector<TLorentzVector> InputPart;
    std::vector<TLorentzVector> InputPart_Cent;
    std::vector<TLorentzVector> InputPart_Fwd;
    std::vector<float> MetVec;
    map<string, double> sphe_var;
    map<string, double> thru_var;
    map<string, double> Cent_Thru_var;
    map<string, double> Cent_UpDownHemisphere_var;
    map<string, double> Cent_corrections;
    double forward_suppression;


    float eta_min = 0;
    float eta_max = 4.9;
    bool skiplepton = true;
    bool include_chglep = true;
    bool include_etmiss = true;


    MetVec.push_back(sum_x);
    MetVec.push_back(sum_y);


    // Find all detectable hadrons (BSJ -> should be bare???)
    InputPart = myUtils.FindParticles(event, eta_min, eta_max, skiplepton, skippart_born);

    // Find all central visible particles
    eta_max = 2.5;
    skiplepton = false;
    InputPart_Cent = myUtils.FindParticles(event, eta_min, eta_max, skiplepton, skippart_bare);

    // Find all visible forward particles
    eta_min = 2.5;
    eta_max = 100.;
    InputPart_Fwd = myUtils.FindParticles(event, eta_min, eta_max, skiplepton, skippart_bare);

    // Get pseudo-global variables (up to eta = 4.9)
    sphe_var = myUtils.TransverseSphericity(InputPart, include_chglep, PromptLeptonBare_Coll, include_etmiss, MetVec);
    thru_var = myUtils.TransverseThrust(InputPart, include_chglep, PromptLeptonBare_Coll, include_etmiss, MetVec);


    // Get central global variables
      // Note: The inputs already include the charged lepton so they don't need to be added by hand here => include_chglep = false
    include_chglep = false;
    Cent_Thru_var = myUtils.TransverseThrust(InputPart_Cent, include_chglep, PromptLeptonBare_Coll, include_etmiss, MetVec);


    // Get corrections to central global variables
    Cent_corrections= myUtils.CentralGlobalCorr(InputPart_Cent,  include_chglep, PromptLeptonBare_Coll, include_etmiss, MetVec);
    forward_suppression = myUtils.ForwardSuppressionCorr(InputPart_Fwd, Cent_corrections["QTC"], Cent_corrections["qrap_cent"]);


    // Up and Down hemisphere quantities
    TVector3 CentralThrustAxis;
    CentralThrustAxis.SetXYZ(Cent_Thru_var["TransvThrustAxisX"],Cent_Thru_var["TransvThrustAxisY"],Cent_Thru_var["TransvThrustAxisZ"]);

    Cent_UpDownHemisphere_var = myUtils.CentralUpAndDownHemisphereQuantities(InputPart_Cent, Cent_corrections["QTC"], CentralThrustAxis);


    // Fill ntuple variables
    glob_TransSpher =  sphe_var["TransvSphericity"];
    glob_TransThrustMaj = 1. - thru_var["TransvThrustMaj"];
    glob_TransThrustMin = thru_var["TransvThrustMin"];

    glob_TransThrustMajRes = 1. - Cent_Thru_var["TransvThrustMaj"] + Cent_corrections["Residue"];
    glob_TransThrustMinRes = Cent_Thru_var["TransvThrustMin"] + Cent_corrections["Residue"];
    glob_TransThrustMajSup = 1. - Cent_Thru_var["TransvThrustMaj"] + forward_suppression;
    glob_TransThrustMinSup = Cent_Thru_var["TransvThrustMin"] + forward_suppression;

    glob_SumMassRes = Cent_UpDownHemisphere_var["CentralMass"] + Cent_corrections["Residue"];
    glob_HeavyMassRes = Cent_UpDownHemisphere_var["HeavyMass"] + Cent_corrections["Residue"];
    glob_SumMassSup = Cent_UpDownHemisphere_var["CentralMass"] + forward_suppression;
    glob_HeavyMassSup = Cent_UpDownHemisphere_var["HeavyMass"] + forward_suppression;

    glob_TotBroadRes = Cent_UpDownHemisphere_var["TotBroadenings"] + Cent_corrections["Residue"];
    glob_WideBroadRes = Cent_UpDownHemisphere_var["WideBroadenings"] + Cent_corrections["Residue"];
    glob_TotBroadSup = Cent_UpDownHemisphere_var["TotBroadenings"] + forward_suppression;
    glob_WideBroadSup = Cent_UpDownHemisphere_var["WideBroadenings"] + forward_suppression;

    glob_SuperSphero = Cent_UpDownHemisphere_var["SuperSpherocity"];



    // Fill the analysis tree
    // ......................


    tree->Fill();

 
    // Clear event-based vectors
    // -------------------------

    top_pt.clear();
    top_eta.clear();
    top_phi.clear();
    top_E.clear();

    neutrino_pt.clear();
    neutrino_eta.clear();
    neutrino_phi.clear();
    neutrino_E.clear();
    neutrino_PdgId.clear();

    muon_bare_pt.clear();
    muon_bare_eta.clear();
    muon_bare_phi.clear();
    muon_bare_E.clear();
    muon_bare_charge.clear();

    muon_dress_pt.clear();
    muon_dress_eta.clear();
    muon_dress_phi.clear();
    muon_dress_E.clear();
    muon_dress_charge.clear();

    muon_born_pt.clear();
    muon_born_eta.clear();
    muon_born_phi.clear();
    muon_born_E.clear();
    muon_born_charge.clear();

    electron_bare_pt.clear();
    electron_bare_eta.clear();
    electron_bare_phi.clear();
    electron_bare_E.clear();
    electron_bare_charge.clear();

    electron_dress_pt.clear();
    electron_dress_eta.clear();
    electron_dress_phi.clear();
    electron_dress_E.clear();
    electron_dress_charge.clear();

    electron_born_pt.clear();
    electron_born_eta.clear();
    electron_born_phi.clear();
    electron_born_E.clear();
    electron_born_charge.clear();

    bjet_bare_pt.clear();
    bjet_bare_eta.clear();
    bjet_bare_phi.clear();
    bjet_bare_E.clear();
    bjet_bare_nPart.clear();

    lightjet_bare_pt.clear();
    lightjet_bare_eta.clear();
    lightjet_bare_phi.clear();
    lightjet_bare_E.clear();
    lightjet_bare_nPart.clear();

    jet_dress_pt.clear();
    jet_dress_eta.clear();
    jet_dress_phi.clear();
    jet_dress_E.clear();
    jet_dress_nPart.clear();

    jet_born_pt.clear();
    jet_born_eta.clear();
    jet_born_phi.clear();
    jet_born_E.clear();
    jet_born_nPart.clear();

    largeRjet_bare_pt.clear();
    largeRjet_bare_eta.clear();
    largeRjet_bare_phi.clear();
    largeRjet_bare_E.clear();

    largeRjet_born_pt.clear();
    largeRjet_born_eta.clear();
    largeRjet_born_phi.clear();
    largeRjet_born_E.clear();

    largeRjet_dress_pt.clear();
    largeRjet_dress_eta.clear();
    largeRjet_dress_phi.clear();
    largeRjet_dress_E.clear();

    bpartonjet_pt.clear();
    bpartonjet_eta.clear();
    bpartonjet_phi.clear();
    bpartonjet_E.clear();

    lightpartonjet_pt.clear();
    lightpartonjet_eta.clear();
    lightpartonjet_phi.clear();
    lightpartonjet_E.clear();

    boson_pt.clear();
    boson_eta.clear();
    boson_phi.clear();
    boson_E.clear();
    boson_ID.clear();

    promptphoton_pt.clear();    // BSJ 
    promptphoton_eta.clear();
    promptphoton_phi.clear();
    promptphoton_E.clear();

    BarekTSplittingScale1_R04.clear();
    BarekTSplittingScale2_R04.clear();
    BarekTSplittingScale3_R04.clear();
    BarekTSplittingScale1_R10.clear();
    BarekTSplittingScale2_R10.clear();
    BarekTSplittingScale3_R10.clear();

    DresskTSplittingScale1_R04.clear();
    DresskTSplittingScale2_R04.clear();
    DresskTSplittingScale3_R04.clear();
    DresskTSplittingScale1_R10.clear();
    DresskTSplittingScale2_R10.clear();
    DresskTSplittingScale3_R10.clear();

    BornkTSplittingScale1_R04.clear();
    BornkTSplittingScale2_R04.clear();
    BornkTSplittingScale3_R04.clear();
    BornkTSplittingScale1_R10.clear();
    BornkTSplittingScale2_R10.clear();
    BornkTSplittingScale3_R10.clear();

    event_weights.clear();

    Top_Coll.clear();
    Vecboson_Coll.clear();
    Promptphoton_Coll.clear();    // BSJ
    PromptLeptonBare_Coll.clear();
    LeptonDress_Coll.clear();
    LeptonBorn_Coll.clear();
    LeptonConversion_Coll.clear();
    FSRPhoton_Coll.clear();
    DressPhoton_Coll.clear();
    Neutrino_Coll.clear();


    TruthBareSmallRJets_Coll.clear();
    TruthBornSmallRJets_Coll.clear();
    TruthDressSmallRJets_Coll.clear();
    TruthBareLargeRJets_Coll.clear();
    TruthDressLargeRJets_Coll.clear();
    TruthBornLargeRJets_Coll.clear();

    BareKtSplittingScale1_R04.clear();
    BareKtSplittingScale2_R04.clear();
    BareKtSplittingScale3_R04.clear();
    BareKtSplittingScale1_R10.clear();
    BareKtSplittingScale2_R10.clear();
    BareKtSplittingScale3_R10.clear();
  
    DressKtSplittingScale1_R04.clear();
    DressKtSplittingScale2_R04.clear();
    DressKtSplittingScale3_R04.clear();
    DressKtSplittingScale1_R10.clear();
    DressKtSplittingScale2_R10.clear();
    DressKtSplittingScale3_R10.clear();
  
    BornKtSplittingScale1_R04.clear();
    BornKtSplittingScale2_R04.clear();
    BornKtSplittingScale3_R04.clear();
    BornKtSplittingScale1_R10.clear();
    BornKtSplittingScale2_R10.clear();
    BornKtSplittingScale3_R10.clear();
    PartonJets_Coll.clear();

    vecbosindex.clear();
    skippart_born.clear();
    skippart_bare.clear();
    skippart_dress.clear();
    InputPart.clear();
    InputPart_Cent.clear();
    InputPart_Fwd.clear();
    MetVec.clear();

  } 

//*************************************************************************************






//*************************************************************************************
// Finishing Code
//*************************************************************************************
void MyAnalysis::finish()
  {

    // Normalize histograms
    // --------------------

    //Double_t norm = lepton_px->GetEntries();
    //lepton_px->Scale(1/norm);


    // Print histograms
    // ----------------


    // Write trees

    tree->Write();

    // Final Print Statements if needed
    // --------------------------------


  }

  //*************************************************************************************






//*************************************************************************************
// Main Code
//*************************************************************************************
int main(int argc, char* argv[])
  {


    //===========================================================================
    // Initialization Phase
    //===========================================================================


    // Safety checks to make sure the code is run properly
    // ---------------------------------------------------

    // Check that correct number of command-line arguments
    // ...................................................

    // if (argc < 2)
    // {
    //   cerr << " Unexpected number of command-line arguments. \n"
    //   << " You are expected to provide at least an input file name. \n"
    //   << " Program stopped! " << endl;
    //   return 1;
    // }


    // Check that the provided file name corresponds to an existing file
    // .................................................................

    // ifstream is(argv[1]);
    // if (!is) {
    //   cerr << " Command-line file " << argv[1] << " was not found. \n"
    //   << " Program stopped! " << endl;
    //   return 1;
    // }


    // Set other running option keys and flags
    // ---------------------------------------


     // Declare StdArg objects with all flags and keys used in command line
     // ...................................................................

  
    StdArg arg(argc,argv);


     // Print the keys and flags used in command line
     // .............................................
 
    cout << "=================== Begin Arguments ===================================================" << endl;
    for (int i=0;i<argc;i++)
      cout << argv[i] << " ";
    cout << endl;
    cout << "=================== End Arguments ===================================================" << endl;
   

     // Enter all possible keys
     // .......................

    arg.keys  << "-outroot" ;
    arg.keys  << "-outhepmc" ;


     // Instantiate variables that will correspond to key values or correspond to flag booleans
     // .......................................................................................

    string root_output    = "OutputHistos";
    string hepmc_output   = "OutputHepMC";


     // Use Process function of StdArg to parse command line and assign input value or bools to variables
     // .................................................................................................
                                                                              
    try 
      { 

        // Process the argument inputed on command line
        // ............................................
    
        arg.Process();


        // Assign value to keys
        // ....................

        // Note: Check if the key has been parsed in command line and if it is, than get the value 
        //       parsed and assign it to the proper variable

        if ( arg.Key("-outroot")   ) root_output      = arg.Get<string>("-outroot");
        if ( arg.Key("-outhepmc")   ) hepmc_output      = arg.Get<string>("-outhepmc");

      } // end try


     // Use the function BadInput of StdArg to report any error on the flags and keys used in command line
     // ..................................................................................................

    catch (StdArg::BadInput) 
      {
        if (argc > 2) cout<< "Input error" <<endl;
    
        // Call usage function in case of error or no parameters
        // .....................................................
    
        Usage(argv[0]);
        return 1;
    
      } // end catch
  

    // Complain if input file not found                                                                                        
    // ................................


    ifstream is(argv[1]);
    if (!is)
      {
        cerr << " Command-line input file " << argv[1] << " was not found. \n"
            << " Program stopped! " << endl;
        return 1;
      }


    // Prepare output file          
    // ...................

    std::string myOutName;
    myOutName = "../WplusJetsAnalysis/pythia-outputs/2025/"+root_output+".root";

    TFile *myfile = TFile::Open(myOutName.c_str(),"recreate");


    std::string myHepMCName;
    myHepMCName = hepmc_output+".dat";




    // Pythia 8 initialization
    // -----------------------


    // Interface for conversion from Pythia8::Event to HepMC event
    // ...........................................................

    HepMC::Pythia8ToHepMC ToHepMC;
    ToHepMC.set_print_inconsistency(false);

    // Specify file where HepMC events will be stored
    // ..............................................

    HepMC::IO_GenEvent ascii_io(myHepMCName, std::ios::out);


    // Declare a pythia object
    // .......................

    Pythia pythia;


    // Specification of Pythia settings
    // ................................

    // Note: These settings are specified in the input file


    pythia.readFile("MyPythia8Simul_RunParameters.cmnd");
  
    pythia.readFile("MyPythia8Simul_PhysicsParameters.cmnd");




    pythia.readFile(argv[1]);


    // Initialize Pythia
    // .................

    pythia.init();


    // Analysis initialization
    // -----------------------

    // Declare user analysis class
    // ...........................

    MyAnalysis myAnalysis;
    ANA_utils myUtilAna;


    // Read in number of event
    // .......................

    int nEvent = pythia.mode("Main:numberOfEvents");

    // Some global variables
    // ---------------------

    int nListEvts = 2;


    // Declare Event Variables
    // -----------------------

    // An event record for parton level particles
    // ..........................................

    // Note: Partons will be taken at the end of the parton evolution
    //       i.e. just before the hadronization.

    Event partonLevelEvent;
    partonLevelEvent.init("Parton Level event record", &pythia.particleData);

  
     // Uncertainty band settings
     // .........................

        // Note: The number of weights is set by the dimension of "UncertaintyBandList" in input
        //       file + 1 for the nominal weight. The string stored are the names of the variations
        //       stored in this List.

    int numOfWeights = pythia.info.nWeights();
    std::cout << "N weights = " << numOfWeights << std::endl;
    vector<double> sumOfWeights;
    vector<string> names;
    vector<string> weightStrings = pythia.settings.wvec("UncertaintyBands:List");

    for (int iWeight=0; iWeight < numOfWeights; ++iWeight)
      {
        names.push_back( (iWeight==0)
                      ? "baseline" : myUtilAna.weightLabel(weightStrings[iWeight-1]));
        sumOfWeights.push_back(0.);
      }


     // Histograms for event weight sum
     // ...............................

        // Note: This histograms will keep track of the Sum of weights for all events for all systematic
        //       variations in order to properly normalize the varied samples. These event weight sums can
        //       also be obtained by summing over all events the weights stored in the ntuple.

    TH1D *sumweights = new TH1D("sumw","Event Sum Weights", numOfWeights, 0, numOfWeights);


    // Initialize the analysis
    // -----------------------

     // Note: This where ntuple variables and histograms are defined

    myAnalysis.init();


    //===========================================================================
    // Loop to Generate Events
    //===========================================================================


    // Begin event loop. Generate event. Skip if error
    // -----------------------------------------------

    for (int iEvent = 0; iEvent < nEvent; ++iEvent)
      {
        if (!pythia.next()) 
          {
            if ( pythia.info.atEndOfFile() ) break;
            else
              {
                continue;
              } 
          }   


        //std::cout << "Event " << iEvent << std::endl;


        // Get the parton-level event
        // --------------------------

        // Declare an Analysis Utilities class object
        // ------------------------------------------

        // Note: To be able to access the functions define there
    
        ANA_utils myUtilsMain;



        myUtilsMain.getPartonLevelEvent(pythia.event, partonLevelEvent);


        // Display some event info
        // -----------------------

        // List first few events
        // .....................

        if (iEvent < nListEvts)
          {
            //pythia.event.list();
            //partonLevelEvent.list();
          }


        // User Analysis of current event
        // ------------------------------

        std::vector<double> EvtWeights;

        for (int iWeight = 0; iWeight < numOfWeights; ++iWeight)
          {

            // Get weight
            // ..........

            double w = pythia.info.weight(iWeight);


            // Eliminate events with a negative weight
            // .......................................

            if (w<0. || w>10.)
              {
                w = 0.;
                std::cout << "Negative or too large weight for variation: " << names[iWeight] << std::endl;
              }


            // Add the weight of the current event to the wsum of weights
            // ..........................................................

            sumOfWeights[iWeight]  += w;

            // Store event weights in a vector
            // ...............................

            EvtWeights.push_back(w);

          }


        // Run the analysis
        // ................

        myAnalysis.analyze(pythia.event, partonLevelEvent, EvtWeights);

        EvtWeights.clear();



        // Produce HepMC output
        // --------------------

        // Construct new empty HepMC event and fill it
        // ...........................................

         // Note: Units will be as chosen for HepMC build, but can be changed by
         //       arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)

        HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
        ToHepMC.fill_next_event( pythia, hepmcevt );


        // Write the HepMC event to file. Done with it
        // ...........................................

        ascii_io << hepmcevt;
        delete hepmcevt;


      }   // End of event loop.


    //===========================================================================
    // Control Output and Run Information
    //===========================================================================


    // Pythia Statistics display
    // -------------------------

    pythia.stat();


    // User finishing
    // --------------

    myAnalysis.finish();

     // Fill sum of weight histogram
     // ............................

    for (int iWeight=0; iWeight < numOfWeights; ++iWeight)
      {
        sumweights->SetBinContent(iWeight,sumOfWeights[iWeight]);
      }

    // Write info to file
    // ------------------

    myfile->Write();
    myfile->Close();

    return 0.;
  }

//*************************************************************************************
