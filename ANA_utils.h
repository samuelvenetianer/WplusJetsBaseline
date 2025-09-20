//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: ANA_utils.h
// 
// Purpose:  Header file for analysis
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// Dependencies (#includes)
// ------------------------

//// C++

#include <cmath>
#include <iostream>
#include <sstream>
#include <set>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>

//// Root

#include "TLorentzVector.h" 
#include "TMatrixD.h"
#include "TDecompSVD.h"

//// FastJet

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "Pythia8Plugins/FastJet3.h"

//// My includes

#include "TruthPart.h"
#include "TruthJets.h"

//// Pythia
#include "Pythia8/Pythia.h"

using namespace std;
using namespace Pythia8;


// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
// -----------------------------------------------------------------------------------
// class ANA_utils:
//
//       Functions to be used in the analysis part of the generation process in order
//       to store the relevant information in the ntuples, ttrees, and histograms.
//  
// -----------------------------------------------------------------------------------
// $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

class ANA_utils {

 public :

  // Constructor and destructor
  // --------------------------

  ANA_utils(){}
  
  ~ANA_utils(){}

  // Declare functions
  // -----------------

  //// Utility functions:

  double  delta_phi(double phi1, double phi2);

  
  //// Event handling functions:

  double compute_weight(const int NB, const double x_sect, const double inst_lum);

  void getPartonLevelEvent( Pythia8::Event& event, Event& partonLevelEvent);

  std::vector<TLorentzVector> FindParticles(Pythia8::Event event, float etamin, float etamax, bool skip_lep, std::vector<int> partskipped);

  string weightLabel(string weightString);
  

  //// Kinematic functions:

  map<string, double> TransverseSphericity(std::vector<TLorentzVector> input, bool addChglep, std::vector<TruthPart> LeptonBare, bool addNeut, std::vector<float> ETmiss);

  map<string, double> TransverseThrust(std::vector<TLorentzVector> input, bool addChglep, std::vector<TruthPart> LeptonBare, bool addNeut, std::vector<float> ETmiss);

  map<string, double> CentralGlobalCorr(std::vector<TLorentzVector> input, bool addChglep, std::vector<TruthPart> LeptonBare, bool addNeut, std::vector<float> ETmiss);

  double ForwardSuppressionCorr(std::vector<TLorentzVector> input, double QTC, double eta_c);

  map<string, double> CentralUpAndDownHemisphereQuantities(std::vector<TLorentzVector> input, double QTC, TVector3 ThrustAxis_c);
  
    
  //// Truth Level Analysis Functions

  void Fill_TruthPart(Pythia8::Event event, int index, TruthPart* p_TruthPart);

  void Get_Tops(Pythia8::Event event, std::vector<TruthPart>* p_Top_Coll);

  void Get_VectorBosons(Pythia8::Event event, std::vector<TruthPart>* p_VecBoson_Coll);

  void Get_PromptPhotons(Pythia8::Event event, std::vector<TruthPart>* p_PromptPhotons_Coll);

  void Get_BarePromptLepton(Pythia8::Event event, std::vector<int> vecboson_index, std::vector<TruthPart>* p_PromptBareLept_Coll, std::vector<TruthPart>* p_ConversionLept_Coll, std::vector<TruthPart>* p_PhotonFSR_Coll, std::vector<TruthPart>* p_Neutrino_Coll);  // BSJ -> should have 5 or 6? (no conversion originally)

  void Get_DressPromptLepton(std::vector<TruthPart> BareLept_Coll, std::vector<TruthPart> PhotonFSR_Coll, std::vector<TruthPart>* p_DressLept_Coll, std::vector<TruthPart>* p_PhotonDress_Coll);    // BSJ

  void Get_BornPromptLepton(Pythia8::Event event, std::vector<int> vecbosindex, std::vector<TruthPart>*  p_LeptonBorn_Coll);
    
  void Bare_Welectron4Vec(Pythia8::Event event, std::vector<TLorentzVector>* p_Born_Coll);

  void Bare_WelectronTruePart(Pythia8::Event event, std::vector<TruthPart>* p_BornElec_Coll);

  void TrueJetsReco(Pythia8::Event event, std::vector<int> partskipped, std::vector<TruthJets>* p_TruthJets_Coll, float ptcut, bool doKtSplitting, std::vector<double>* p_KtSplittingScale1 = NULL, std::vector<double>* p_KtSplittingScale2 = NULL, std::vector<double>* p_KtSplittingScale3 = NULL, double R = 0.4);

  void PartonJetsReco(Pythia8::Event event, Pythia8::Event partonevent, std::vector<TruthJets>* p_PartonJets_Coll, float ptcut);

  void Get_BottomQuarks(Pythia8::Event event, std::vector<TruthPart>* p_BQuarks_Coll);

  void Get_BottomHadrons(Pythia8::Event event, std::vector<TruthPart> BQuark_Coll, std::vector<TruthPart>* p_BHadrons_Coll);

};
