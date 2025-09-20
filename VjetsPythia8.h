//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: VjetsPythia8.h
// 
// Purpose:  Header file for the Pythia8 production and analysis
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// Dependencies (#includes)
// ========================

//// Pythia

#include "Pythia8/Pythia.h"


//// HepMC

#include "Pythia8Plugins/HepMC2.h"


//// FastJet

//   Note: The FastJet3.h header enables automatic initialisation of fastjet::PseudoJet 
//         objects from Pythia8 Particle and Vec4 objects, as well as advanced features
//         such as access to (a copy of) the original Pythia 8 Particle directly from
//         the PseudoJet, and fastjet selectors that make use of the Particle properties.
//         See the extensive comments in the header file for further details and examples.

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "Pythia8Plugins/FastJet3.h"

//// Root

#include "TH1.h"
#include "TVirtualPad.h"     // Interactive graphics.
#include "TApplication.h"
#include "TFile.h"           // Saving file.
#include "TTree.h"
#include "TNtuple.h"
#include "TLorentzVector.h" 
#include "TMath.h"
#include "TRandom3.h"
#include "TRandom.h" 

//// My includes

#include "ANA_utils.h"
#include "TruthPart.h"
#include "TruthJets.h"

//// C/C++ includes

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include "StdArg.hpp"


using namespace Pythia8;

// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
// class MyAnalysis
//
//   Note: Main analysis class definition
//  
// -----------------------------------------------------------------------------------
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

class MyAnalysis {

 public:

  // Constructor and destructor
  // --------------------------

  MyAnalysis() {}

  // Declare functions
  // -----------------
  // Note: Utility functions are defined in ANA_utils.cc, while the main generator
  //       functions and the analysis are defined in VjetsPythia.cc.

  // Initialization actions
  // ......................
  
  void init();

  // Analysis of each new event
  // ..........................
  
  void analyze(Event& event, Event& partonevent, std::vector<double> EventWeights, int i_hardjet);

  // Show final results
  // ..................
  
  void finish();

  // Compute and save kT splitting scales
  bool BareKtSplittingScales = true;
  bool BornKtSplittingScales = true;
  bool DressKtSplittingScales = true;         // BSJ
  
 private:

  // Declare variables and objects that span init - analyze - finish
  // ---------------------------------------------------------------
  
  // Global variables
  // ................
  
  int nEvt;
  int nEventAccept;
  int vetoCount[4];

  bool firstEvent;
  bool debug;

  // Lepton and jet level
  // ....................

  string level;
  
  // Declare Objects and Vectors of Objects
  // --------------------------------------

  std::vector<TruthPart> Top_Coll;
  std::vector<TruthPart>* p_Top_Coll;

  std::vector<TruthPart> Vecboson_Coll;
  std::vector<TruthPart>* p_Vecboson_Coll;

  std::vector<TruthPart> Promptphoton_Coll;
  std::vector<TruthPart>* p_Promptphoton_Coll;

  std::vector<TruthPart> PromptLeptonBare_Coll;
  std::vector<TruthPart>* p_PromptLeptonBare_Coll;

  std::vector<TruthPart> LeptonDress_Coll;          // BSJ
  std::vector<TruthPart>* p_LeptonDress_Coll;

  std::vector<TruthPart> LeptonBorn_Coll;
  std::vector<TruthPart>* p_LeptonBorn_Coll;

  std::vector<TruthPart> LeptonConversion_Coll;       // BSJ
  std::vector<TruthPart>* p_LeptonConversion_Coll;

  std::vector<TruthPart> FSRPhoton_Coll;
  std::vector<TruthPart>* p_FSRPhoton_Coll;

  std::vector<TruthPart> DressPhoton_Coll;        // BSJ
  std::vector<TruthPart>* p_DressPhoton_Coll;

  std::vector<TruthPart> Neutrino_Coll;
  std::vector<TruthPart>* p_Neutrino_Coll;

  std::vector<TruthJets> TruthBareSmallRJets_Coll;
  std::vector<TruthJets>* p_TruthBareSmallRJets_Coll;

  std::vector<TruthJets> TruthDressSmallRJets_Coll;     // BSJ
  std::vector<TruthJets>* p_TruthDressSmallRJets_Coll;

  std::vector<TruthJets> TruthBornSmallRJets_Coll;
  std::vector<TruthJets>* p_TruthBornSmallRJets_Coll;

  std::vector<TruthJets> TruthBareLargeRJets_Coll;
  std::vector<TruthJets>* p_TruthBareLargeRJets_Coll;

  std::vector<TruthJets> TruthDressLargeRJets_Coll;       // BSJ
  std::vector<TruthJets>* p_TruthDressLargeRJets_Coll;

  std::vector<TruthJets> TruthBornLargeRJets_Coll;
  std::vector<TruthJets>* p_TruthBornLargeRJets_Coll;

  std::vector<TruthJets> PartonJets_Coll;    // Small-R only
  std::vector<TruthJets>* p_PartonJets_Coll; // Small-R only

  // kT splitting scales using as inputs the constituents of small-R (light and b-) jets
  std::vector<double> BareKtSplittingScale1_R04;
  std::vector<double> BareKtSplittingScale2_R04;
  std::vector<double> BareKtSplittingScale3_R04;
  std::vector<double>* p_BareKtSplittingScale1_R04;
  std::vector<double>* p_BareKtSplittingScale2_R04;
  std::vector<double>* p_BareKtSplittingScale3_R04;
  std::vector<double> DressKtSplittingScale1_R04;           // BSJ
  std::vector<double> DressKtSplittingScale2_R04;           
  std::vector<double> DressKtSplittingScale3_R04;
  std::vector<double>* p_DressKtSplittingScale1_R04;
  std::vector<double>* p_DressKtSplittingScale2_R04;
  std::vector<double>* p_DressKtSplittingScale3_R04;
  std::vector<double> BornKtSplittingScale1_R04;
  std::vector<double> BornKtSplittingScale2_R04;
  std::vector<double> BornKtSplittingScale3_R04;
  std::vector<double>* p_BornKtSplittingScale1_R04;
  std::vector<double>* p_BornKtSplittingScale2_R04;
  std::vector<double>* p_BornKtSplittingScale3_R04;

  // kT splitting scales using as inputs the constituents of large-R jets
  std::vector<double> BareKtSplittingScale1_R10;
  std::vector<double> BareKtSplittingScale2_R10;
  std::vector<double> BareKtSplittingScale3_R10;
  std::vector<double>* p_BareKtSplittingScale1_R10;
  std::vector<double>* p_BareKtSplittingScale2_R10;
  std::vector<double>* p_BareKtSplittingScale3_R10;
  std::vector<double> DressKtSplittingScale1_R10;       // BSJ
  std::vector<double> DressKtSplittingScale2_R10;
  std::vector<double> DressKtSplittingScale3_R10;
  std::vector<double>* p_DressKtSplittingScale1_R10;
  std::vector<double>* p_DressKtSplittingScale2_R10;
  std::vector<double>* p_DressKtSplittingScale3_R10;
  std::vector<double> BornKtSplittingScale1_R10;
  std::vector<double> BornKtSplittingScale2_R10;
  std::vector<double> BornKtSplittingScale3_R10;
  std::vector<double>* p_BornKtSplittingScale1_R10;
  std::vector<double>* p_BornKtSplittingScale2_R10;
  std::vector<double>* p_BornKtSplittingScale3_R10;
  
  // Histograms
  // ..........

  // TTree variables
  // ...............
  // Note: For cases with more than one particle of the same kind, produce a vector for each branch
  
  TTree *tree;

  int NumHardJets, nTop, nNeutrino, nMuonBare, nElectronBare, nMuonDress, nElectronDress, nMuonBorn, nElectronBorn, nLightjetBare, nBjetBare, nJetDress, nJetBorn, nLargeRjetDress, nLargeRjetBare, nLargeRjetBorn, nBoson, nPromptPhotons, nLightpartonjet, nBpartonjet;  // BSJ
  int nMuonReco, nElectronReco, nJetReco;       // BSJ

  double Met, Met_phi;
  double RecoMet, RecoMet_phi;        // BSJ

  double glob_TransSpher, glob_TransThrustMaj, glob_TransThrustMin;
  double glob_TransThrustMajRes, glob_TransThrustMinRes, glob_TransThrustMajSup, glob_TransThrustMinSup;
  double glob_SumMassRes, glob_HeavyMassRes, glob_SumMassSup, glob_HeavyMassSup, glob_TotBroadRes;
  double glob_WideBroadRes, glob_TotBroadSup, glob_WideBroadSup, glob_SuperSphero;

  std::vector<int> neutrino_PdgId;
  std::vector<int> lightjet_bare_nPart;
  std::vector<int> jet_dress_nPart;    //BSJ
  std::vector<int> jet_born_nPart;
  std::vector<int> bjet_bare_nPart;
  std::vector<int> boson_ID;

  std::vector<float> top_pt, top_eta, top_phi, top_E;
  std::vector<float> neutrino_pt, neutrino_eta, neutrino_phi, neutrino_E;
  std::vector<float> muon_bare_pt, muon_bare_eta, muon_bare_phi, muon_bare_E, muon_bare_charge;
  std::vector<float> muon_dress_pt, muon_dress_eta, muon_dress_phi, muon_dress_E, muon_dress_charge;       // BSJ
  std::vector<float> muon_born_pt, muon_born_eta, muon_born_phi, muon_born_E, muon_born_charge;
  std::vector<float> electron_bare_pt, electron_bare_eta, electron_bare_phi, electron_bare_E, electron_bare_charge;
  std::vector<float> electron_dress_pt, electron_dress_eta, electron_dress_phi, electron_dress_E, electron_dress_charge;     // BSJ
  std::vector<float> electron_born_pt, electron_born_eta, electron_born_phi, electron_born_E, electron_born_charge;
  std::vector<float> lightjet_bare_pt, lightjet_bare_eta, lightjet_bare_phi, lightjet_bare_E;
  std::vector<float> jet_dress_pt, jet_dress_eta, jet_dress_phi, jet_dress_E;         // BSJ
  std::vector<float> bjet_bare_pt, bjet_bare_eta, bjet_bare_phi, bjet_bare_E;
  std::vector<float> jet_born_pt, jet_born_eta, jet_born_phi, jet_born_E;
  std::vector<float> largeRjet_bare_pt, largeRjet_bare_eta, largeRjet_bare_phi, largeRjet_bare_E;
  std::vector<float> largeRjet_dress_pt, largeRjet_dress_eta, largeRjet_dress_phi, largeRjet_dress_E;   // BSJ
  std::vector<float> largeRjet_born_pt, largeRjet_born_eta, largeRjet_born_phi, largeRjet_born_E;
  std::vector<double> BarekTSplittingScale1_R04, BarekTSplittingScale2_R04, BarekTSplittingScale3_R04;
  std::vector<double> BarekTSplittingScale1_R10, BarekTSplittingScale2_R10, BarekTSplittingScale3_R10;
  std::vector<double> DresskTSplittingScale1_R04, DresskTSplittingScale2_R04, DresskTSplittingScale3_R04;    // BSJ
  std::vector<double> DresskTSplittingScale1_R10, DresskTSplittingScale2_R10, DresskTSplittingScale3_R10;
  std::vector<double> BornkTSplittingScale1_R04, BornkTSplittingScale2_R04, BornkTSplittingScale3_R04;
  std::vector<double> BornkTSplittingScale1_R10, BornkTSplittingScale2_R10, BornkTSplittingScale3_R10;
  std::vector<float> lightpartonjet_pt, lightpartonjet_eta, lightpartonjet_phi, lightpartonjet_E;
  std::vector<float> bpartonjet_pt, bpartonjet_eta, bpartonjet_phi, bpartonjet_E;
  std::vector<float> boson_pt, boson_eta, boson_phi, boson_E;
  std::vector<float> promptphoton_pt, promptphoton_eta, promptphoton_phi, promptphoton_E;
  
  std::vector<double> event_weights;

  std::vector<float> muon_reco_pt, muon_reco_eta, muon_reco_phi, muon_reco_E, muon_reco_charge;            // BSJ
  std::vector<float> electron_reco_pt, electron_reco_eta, electron_reco_phi, electron_reco_E, electron_reco_charge;
  std::vector<float> jet_reco_pt, jet_reco_eta, jet_reco_phi, jet_reco_E;
  
};
