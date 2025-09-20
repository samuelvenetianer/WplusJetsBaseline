//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// File: ANA_utils.cc
// 
// Purpose: Provide functions to fill physics objects with proper corrections and
//          systematics. Utility functions for various kinematic calculations too.
//
//
//&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

// Dependencies (#includes)
// ------------------------

#include "ANA_utils.h"

//// C++

#include <cmath>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cstring>
#include <stdio.h>
#include <ios>
#include <fstream>
#include <cstdio>

using namespace std;


//*************************************************************************************************************
//*************************************************************************************************************
//
//    Utility Functions:
//    =================
//
//         -Calculate Delta Phi between 0 and Pi
//
//
//*************************************************************************************************************
//*************************************************************************************************************

// =========================================================================================================================
double  ANA_utils::delta_phi(double phi1, double phi2) {
//
//    Note: These are functions used in kinematic calculations
//
// =========================================================================================================================

  const double PI=2.0*acos(0.);
  const double TWOPI=2.0*PI;

  if (phi1<0) phi1= phi1+TWOPI;
  if (phi2<0) phi2= phi2+TWOPI;

  double PHI=fabs(phi1-phi2);

  return (PHI<=PI)? PHI : TWOPI-PHI;
}

// =========================================================================================================================





//*************************************************************************************************************
//*************************************************************************************************************
//
//    Event Handling Functions:
//    ========================
//
//         -compute_weight: Weight event for lumi
//
//         -getPartonLevelEvent: fill a pythia Event object with particles at parton level
//
//         -FindParticles: find final state particles and store them in a 4-vector container
//
//         -weightLabel: Small helper function to get uncertainty variation names.
//
//*************************************************************************************************************
//*************************************************************************************************************


// ============================================================================================================
double ANA_utils::compute_weight(const int NB, const double x_sect, const double instLum) {
//
// To compute the event weigth allowing to combine various samples
//
// ============================================================================================================
  double weight;

  weight = x_sect/NB*instLum; 
  
  return weight;
  
} 

//------------------------------------------------------------------------------------------------------------



// ============================================================================================================
void  ANA_utils::getPartonLevelEvent( Event& event, Event& partonLevelEvent) {
//
// A generic routine to extract the particles that existed right before the hadronization machinery is invoked
// by Pythia. This is useful to get parton-level quantities, but after the entire shower history, and not
// randomly between the hard process and the end-point of the evolution.
//
// ============================================================================================================

  partonLevelEvent.reset();

  // Loop over the entire event to select parton level particles just before hadronization
  // -------------------------------------------------------------------------------------
 
  for (int i = 0; i < event.size(); ++i)
    {
      bool accept = false;

      // Only partons after full evolution and before hadronization are kept
      // ...................................................................
      
      if (event[i].isFinalPartonLevel()) accept = true;

      // Don't keep neutrinos and charged leptons because they come from hard process
      // ............................................................................

      int idAbs = event[i].idAbs();
      if (idAbs >10 && idAbs < 17) accept = false;

      // Reject particles outside calorimeter acceptance
      // ...............................................

      if (event[i].eta() > 4.9) accept = false;

      if (accept == true)
	{
	  int iNew = partonLevelEvent.append( event[i] );

          // Set copied properties more appropriately
	  // Note: Set a positive status, original location as "mother", and with no daughters.

	  partonLevelEvent[iNew].statusPos();
	  partonLevelEvent[iNew].mothers( i, i);
	  partonLevelEvent[iNew].daughters( 0, 0);
	}
    }
}

//--------------------------------------------------------------------------



// ============================================================================================================
std::vector<TLorentzVector> ANA_utils::FindParticles(Pythia8::Event event, float etamin, float etamax, bool skip_lep, std::vector<int> partskipped) {
//
// This function simply loop over all the particles in the Pythia event record and store the final state particles
// in a vector of 4-vectors (TLorentzVector), one for each relevant particle. Some particles can be skipped,
// depending on their pseudo-rapidity reach (etamax) or if they are charged lepton decay product of vector bosons.
// This is controled by a boolean flag. Neutrino are not included here. If they are needed in a global event
// observable, the Etmiss of neutrino transverse momentum vector will be added directly to the relevant functions.
//
//    Note: These serve as input to global event variables, and other quantities that are obtained from the full
//          final state event records (modulo possible exceptions mentioned above). For more details, see the
//          relevant global event variable function.
//
// ============================================================================================================

  std::vector<TLorentzVector> InputPartColl;

  for (int i = 0; i < event.size(); ++i)
    {

      // Final state only
      if (!event[i].isFinal())        continue;

      // No neutrinos
      if (event[i].idAbs() == 12 || event[i].idAbs() == 14 || event[i].idAbs() == 16)     continue;

      // Only |eta| < eta_max
      if (fabs(event[i].eta()) > etamax || fabs(event[i].eta()) < etamin) continue;

      // Do not include the stable particles listed in the partskipped vector if skipe_lep is set to true
      // Note: This is used to not include the decay product of the W, Z, H and top for example

      bool skip_i = false;

      for (int j = 0; j < partskipped.size(); j++){
        if (partskipped[j] == i) skip_i = true;
      }

      if (skip_lep == true && skip_i == true) continue;

      TLorentzVector v;
      double vPt  = event[i].pT();
      double vEta = event[i].eta();
      double vPhi = event[i].phi();
      double vM = event[i].m();

      v.SetPtEtaPhiM(vPt,vEta,vPhi,vM);

      InputPartColl.push_back(v);
    }

  return InputPartColl;

}

//--------------------------------------------------------------------------



// ============================================================================================================
string ANA_utils::weightLabel(string weightString) {
//
// Small helper function to get uncertainty variation names.
//
// ============================================================================================================

  // Strip leading whitespace
  // ------------------------

  weightString.erase(0, weightString.find_first_not_of(" \t\n\r\f\v"));

  // Find first blank and use this to isolate weight label
  // -----------------------------------------------------

  int iBlank = weightString.find(" ", 0);
  return weightString.substr(0, iBlank);

}

//--------------------------------------------------------------------------





//*************************************************************************************************************
//*************************************************************************************************************
//
//    Global Quantity Functions:
//    =========================
//
//         -TransverseSphericity: Calculate the transverse Sphericity on final state particles.
//
//         -TransverseThrust: Calculate the transverse Thrust on final state particles.
//
//         -CentralGlobalCorr: Calculate corrections to apply to central quantites to make them global
//
//         -ForwardSuppressionCorr: Calculate exponential suppression factor for forward particles
//
//         -CentralUpAndDownHemisphereQuantities: Calculate central mass, jet-broadenings, and superspherocity
//
//
//
//*************************************************************************************************************
//*************************************************************************************************************

// ============================================================================================================
map<string, double> ANA_utils::TransverseSphericity(std::vector<TLorentzVector> input, bool addChglep, std::vector<TruthPart> LeptonBare, bool addNeut, std::vector<float> ETmiss) {
//
// Reconstruct the transverse sphericity of each event using final state particles. Different versions of these
// quantities could be calculated. Here only the transverse sphericity is calculated because the total sphericity
// makes sense only in the cms of the collision and we cannot reconstruct this system with a neutrino in the final
// state. In addition, underlying event and non-4pi coverage rends cms boost calculation hard to exactly obtain
// even if there is no real ETmiss in the event.  Highly boosted collision will look not spherical even if they
// are, in the cms the event. The full event sphericity in the lab frame is therefore biased, and we can't boost
// to the cms beyond crude approximations that we are not doing here. Choice of including charged lepton and ETmiss
// can be made at the analysis level. This function can therefore accomodate:
//
//       -Transverse calculation using all visible and Met particles
//       -Transverse calculation using only visible particles
//       -Transverse calculation using only hadronic particles
//
// In addition, by definition, even at LEP, sphericity is NOT infrared safe. This makes the observable more
// sensitive to parton shower and hadronization. It is also not resummable at NLL in phase space regions where
// the sphericity is small and the logs are large. It is impossible to make reliable pQCD predictions for this
// observable beyond the Leading Order. Anyway, pythia is not ideal for describing global event observables
// when the value of the variable is large and is better described by fixed-order ME calculations because
// Pythia is only LO in 2->2 process and therefore does not include a ME description of higher jet multiplicity
// final state for which sphericity takes larger values.
//
// As a consequence, the discriminative power between mercedes-like patterns, typically modeled by ME, and
// collinear splitting topologies, of the realm of QCD radiation effects (modeled by PS) is largely reduced.
// It is therefore not the ideal observable for measuring the modeling of parton shower and QCD emission
// effects. We would get better sensitivity with other global event observables, but it is still interesting
// to implement sphericity to see what it gives. Note however that none of the other global event variable
// proposes something that captures the sphericity event topology as was intended by the sphericity variable.
// Banfi, Salam, and Zanderighi therefore proposed an alternative variable, the SuperSphericity which should
// capture the real distribution of object in the 3D volume, while keeping resummed accuracy incorporated in
// analytic predictions or in parton shower. This has yet to be studied. It requires the thrust axis for being
// calculated, as well as the projection on UP or DOWN hemisphere. This is therefore computed in a different
// function.
//
//    Note 1: For reference: Banfi, Salam and Zanderighi, "Phenomenology of event shapes at hadron colliders",
//                         JHEP 1006 (2010) 038, arXiv:1001.4082 [hep-ph];
//                         Banfi, Salam and Zanderighi, "Resummed event shapes at hadron-hadron colliders",
//                         JHEP 0408 (2004) 062, arXiv:hep-ph/0407287v3.
//
//    Note 2: The transverse sphericity can be calculated with jets or partons as input, provided jets are
//            4-vectors used to build a vector of these objects to serve as input. The function is enough
//            general to accomodate for such possibility.
//
//
// ============================================================================================================
  
  // Declare variables internal to this function
  // -------------------------------------------

  map<string, double> Variables;

  Variables["TransvSphericity"] = -999.*1000.;

  TMatrixD TransvMomentumTensor(2,2);
  double   PT2Sum = 0;

  double TransvSphericity = -1;

  // Loop over input objects
  // -----------------------

  for (int i = 0; i < input.size(); ++i) {

    // Calculate tensors
    // -----------------

    // Form the transverse tensor
    // ..........................

    TransvMomentumTensor(0,0) += input[i].Px()*input[i].Px();
    TransvMomentumTensor(0,1) += input[i].Px()*input[i].Py();
    TransvMomentumTensor(1,0) += input[i].Py()*input[i].Px();
    TransvMomentumTensor(1,1) += input[i].Py()*input[i].Py();


    // Include charged leptons from Vector Bosons
    // ..........................................

    if (addChglep == true) {
      for (int ilep = 0; ilep < LeptonBare.size(); ++ilep) {
        TransvMomentumTensor(0,0) += LeptonBare[ilep].Px()*LeptonBare[ilep].Px();
        TransvMomentumTensor(0,1) += LeptonBare[ilep].Px()*LeptonBare[ilep].Py();
        TransvMomentumTensor(1,0) += LeptonBare[ilep].Py()*LeptonBare[ilep].Px();
        TransvMomentumTensor(1,1) += LeptonBare[ilep].Py()*LeptonBare[ilep].Py();
      }
    }

    // Include Etmiss from neutrinos
    // .............................

    if (addNeut == true){
      TransvMomentumTensor(0,0) += ETmiss[0]*ETmiss[0];
      TransvMomentumTensor(0,1) += ETmiss[0]*ETmiss[1];
      TransvMomentumTensor(1,0) += ETmiss[1]*ETmiss[0];
      TransvMomentumTensor(1,1) += ETmiss[1]*ETmiss[1];
    }

  } // end loop on input objects

  // Calculate the denominators
  // ..........................

  PT2Sum = TransvMomentumTensor(0,0)+ TransvMomentumTensor(1,1);

  // Calculate Transverse Sphericity and Acoplanarity
  // ------------------------------------------------

  // Normalize the tensor
  // ....................

  if(PT2Sum > 0) {
    
    const double inv_PT2Sum = 1. / PT2Sum;

    for(int i=0; i<2; i++) {
      for(int j=0; j<2; j++) {
        TransvMomentumTensor(i,j) *= inv_PT2Sum;
      }
    } // End on loops on tensor components

    // Calculate eigenvalues
    // .....................

    TDecompSVD * aTSVD = new TDecompSVD(TransvMomentumTensor);
    TVectorD TLambda = aTSVD->GetSig();

    // Get Sphericity and Acoplanarity
    // ...............................

    TransvSphericity = 2.0*TLambda[1] / (TLambda[0]+TLambda[1]);

    delete aTSVD;

    // Fill the map
    // ............

    Variables["TransvSphericity"] = TransvSphericity;

  }

  // Return variables
  // ----------------

  return Variables;

}
//-------------------------------------------------------------------------------------------------------------



// ============================================================================================================
map<string, double> ANA_utils::TransverseThrust(std::vector<TLorentzVector> input, bool addChglep, std::vector<TruthPart> LeptonBare, bool addNeut, std::vector<float> ETmiss) {
//
// Similarly as for sphericity, only transverse thrust is not biased by the boost of the center-of-mass system
// and thus only transverse thrust makes sense at hadron colliders. It can be computed with or without charged
// leptons and/or neutrinos (ETmiss). Flags are used to decide. The outputs of the function are:
//
//       - Transverse Thrust Major value
//       - Transverse Thrust Major axis vector
//       - Transverse Thrust Minor value
//       - Transverse Thrust Minor axis vector
//
// Corrections to make this function infrared safe and global will follow from other functions.
//
//           Note: The algorithm is adapted from:
//                   ATLAS: https://acode-browser1.usatlas.bnl.gov/lxr/source/athena/Reconstruction/Jet/JetSubStructureUtils/Root/Thrust.cxx
//                   PYTHIA: HEP 05 (2006) 026, hep-ph/0603175, pages 523-525.
//
// ============================================================================================================

  // Declare variables internal to this function
  // -------------------------------------------

  // Variables to be returned by the function
  // ........................................

  map<string, double> Variables;
  Variables["TransvThrustMin"] = -999. * 1000.;
  Variables["TransvThrustMaj"] = -999. * 1000.;
  Variables["TransvThrustAxisX"] = -999. * 1000.;
  Variables["TransvThrustAxisY"] = -999. * 1000.;

  // Define a vector for the thrust principal axis
  // .............................................

  TVector3 transvthrust(0,0,0);

  // Other internal variables
  // ........................

  double transvthrust_major = -1;
  double transvthrust_minor = -1;

  int agree = 0;
  int disagree = 0;

  // Define vectors and arrays for testing the convergence of the minimization process
  // .................................................................................
  // Note: The idea here is to find the vector maximizing the thrust by starting from the most energetic
  //       input object as a seed, and comparing the max obtained with the ones that follow from different
  //       seeds, namely from using the 2nd, 3rd and 4th most energetic input objects. If the max reached
  //       is not the same, the procedure is repeated with 16 combinations of these 4 input objects as
  //       seeds to find out which is the real max.

  // Principal thrust axis for each seed
  TVector3 nT_0[20];

  // Seed variations
  short add0[20] = { 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1,-1,-1,-1,-1,-1,-1,-1,-1 };
  short add1[20] = { 0, 1, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1, 1, 1, 1, 1,-1,-1,-1,-1 };
  short add2[20] = { 0, 0, 1, 0, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1, 1, 1,-1,-1 };
  short add3[20] = { 0, 0, 0, 1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1, 1,-1 };

  // Modify the input if chaged leptons from Vector Boson and/or neturinos are to be included
  // ----------------------------------------------------------------------------------------

  // Get all particles obtained from ANA_utils::FindParticles
  // ........................................................

  std::vector<TLorentzVector> input_fin;

  for (int k=0; k<input.size(); k++) {
    input_fin.push_back(input[k]);
  }

  // Include charged leptons from Vector Bosons
  // ..........................................
  // Note: They are stored in a vector of TruthPart, we need to copy them first in a TLorentzVector

  if (addChglep == true) {
    for (int ilep = 0; ilep < LeptonBare.size(); ++ilep) {
      TLorentzVector Lepton_i;
      Lepton_i.SetPtEtaPhiE(LeptonBare[ilep].Pt(), LeptonBare[ilep].Eta(), LeptonBare[ilep].Phi(), LeptonBare[ilep].E());
      input_fin.push_back( Lepton_i);
    }
  }

  // Include neutrino using ETmiss
  // .............................
  // Note: We assume that pz and mass are null. This makes the energy wrong for the 4-vector. Below
  //       only the 3-momentum is used, so it is not a problem.

  if (addNeut == true) {
    TLorentzVector ETmissPart;
    ETmissPart.SetPxPyPzE(ETmiss[0], ETmiss[1], 0., sqrt(ETmiss[0]*ETmiss[0]+ETmiss[1]*ETmiss[1]));
    input_fin.push_back(ETmissPart);
  }

  // Set the seed, allowing for different starting points
  // ----------------------------------------------------

  // The 4 first initial seeds are the 4 most energetic input objects
  // ................................................................

  std::vector<TLorentzVector> v_copy(4);
  for(int i=0; i<4; i++) {
    v_copy[i] = input_fin[i];
  }

  // Find the first and last input objects in the list of input
  // ..........................................................

  std::vector<TLorentzVector>::const_iterator iBeg = input_fin.begin();
  std::vector<TLorentzVector>::const_iterator iEnd = input_fin.end();

  // Number of tests for convergence
  // ...............................
  // Note 1: If the 4 first seeds converge to the same unit vector, that will be declared as the Thrust axis
  //         otherwise, more seeds must be tested to find convergence. At most 20 different seeds will be tested,
  //         with less than this if the number of input objects in the thrust calculation is 4 < N < 19.
  // Note 2: The std member "distance" calculates the number of elements between the two iterators iBeg and iEnd

  int n_tests = 0;
  int max_tests = min<int>(20, distance(iBeg, iEnd));

  // Find the Thrust Principal Axis using PYTHIA's method
  // ----------------------------------------------------
  // Note: The Thrust axis is the unit vector vec{n} for which the max of Sum | vec{n}.vec{p_i} | / Sum |p_i| is obtained.
  //       This max is the Major Thrust value. There are multiple ways to implement such optimization process, but typically
  //       they require from 4n^2 to 2^(n-1) possibilities to be tested for exact maximization. Such approaches are
  //       prohibitive if the number of input objects is large (e.g. ~10). PYTHIA authors used an alternative method which
  //       is much faster, but which does not guarantee that a global max is found. Trying with different seeds significantly
  //       increases the chance that the max is not just a local max, and this is the method adopted here. According to this
  //       approach, the thrust axis is given by the vector vec{n} for which the following iterative process converge to
  //       a given value, i.e. that the thrust axis is given by vec{n}_(j) when vec{n}_(j+1) = vec{n}_(j), where:
  //
  //                   vec{n}_(j+1) = Sum_i [e(vec{n}_j . p_i) p_i ] / || Sum_i [e(vec{n}_j . p_i) p_i ] ||
  //
  //       where e(vec{n}_j . p_i) = 1 when vec{n}_j . p_i > 0, and -1 otherwise. It can be shown that T_(j+1) >= T_j
  //       and therefore the vec{n} that max T is found after only a few iteration (normally 2 to 4). It is to avoid
  //       convergence toward a local max that different seeds are tried. This is the approach used here.

  // Loop over the seeds to be tested
  // ................................

  do {
    // Re-initialized the unit vector to be used as the thrust axis to (0,0,0) at each pass (each test)

    nT_0[n_tests]=TVector3(0,0,0);

    // Set the seed
    // ............
    // Note: For the first test, the seed is the first input object; from the 5th seed, it is a combination of the
    //       4 first input object starting from all possible combination of axis with only 4 input object.

    nT_0[n_tests] +=
      add0[n_tests] * TVector3(v_copy.at(0).Px(), v_copy.at(0).Py(), v_copy.at(0).Pz()) +
      add1[n_tests] * TVector3(v_copy.at(1).Px(), v_copy.at(1).Py(), v_copy.at(1).Pz()) +
      add2[n_tests] * TVector3(v_copy.at(2).Px(), v_copy.at(2).Py(), v_copy.at(2).Pz()) +
      add3[n_tests] * TVector3(v_copy.at(3).Px(), v_copy.at(3).Py(), v_copy.at(3).Pz());

    // Transverse thrust: set the z-component of the seed to 0.

    nT_0[n_tests].SetZ(0);

    // Normalize the seed vector
    // Note: Divide each component by the magnitude of the vector

    if (nT_0[n_tests].Mag() > 0) nT_0[n_tests] *= 1/nT_0[n_tests].Mag();

    // Find the thrust principal axis for iteration j
    // ..............................................

    int loop = 0;
    bool run = false;
    bool runT = false;

    // Calculate: Sum_i [e(vec{n}_j . p_i) p_i ]
    do {
      TVector3 nT_1(0,0,0);

      for (std::vector<TLorentzVector>::const_iterator i = iBeg; i != iEnd; ++i) {
        const TLorentzVector &itr = *i;

        if ((itr).Px() * nT_0[n_tests].x() +
            (itr).Py() * nT_0[n_tests].y() +
            (itr).Pz() * nT_0[n_tests].z() > 0)
          nT_1 += TVector3((itr).Px(), (itr).Py(), (itr).Pz());
        else
          nT_1 -= TVector3((itr).Px(), (itr).Py(), (itr).Pz());
      }

      // Set the z-component to 0 for transverse thrust
      nT_1.SetZ(0);

      // Normalize to get vec{n}_(j+1)
      if (nT_1.Mag() > 0) nT_1 *= 1/nT_1.Mag();

      // Compare iteration j+1 with iteration j: we found the thrust axis if no change is found
      // ........................................................................................
      // Note 1: If vec{n}^(j+1) != vec{n}^j, try at most ten times ( to make sure to converge for
      //         each seed because the thrust axis has two fold ambiguity and 2 to 4 iterations are
      //         required in general).
      // Note 2: If we have convergence, the vector nT_1 is the major axis of the transverse thrust,
      //         otherwise keep iterating after having changed the seed to nT_1 of this iteration.

      run = (nT_0[n_tests] != nT_1) && (-nT_0[n_tests] != nT_1) && loop++ < 10;

      // Set the thrust axis to its jth iteration
      if (run) nT_0[n_tests] = nT_1;

   } while (run);

   // Check convergence of the found thrust axis
   // ..........................................
   // Note 1: Compare the results of the new seeds with the first seed. How many seeds agree (and disagree) with first result ?
   //         Stop if four first tries give same result (no test for first try!). If not, try at most max_tests combinations.
   //         The full set of seeds will be calculated if any of the thrust or transverse thrust seed disagree with the first
   //         seed considered.
   // Note 2: The thrust axis, amoung all the axis found with the different seeds, that will constitute THE thrust axis, is the
   //         one that actually maximize the Major Thrust value. This final choice is therefore postpone later. However, if the
   //         first 4 seeds converge to the same value, there is no need to calculate the axis for other seeds.
   // Note 3: Make sure to account for thrust sign ambiguity

   if (n_tests > 0 && (nT_0[0] == nT_0[n_tests] || nT_0[0] == -nT_0[n_tests])) agree++;
   if (n_tests > 0 &&  nT_0[0] != nT_0[n_tests] && nT_0[0] != -nT_0[n_tests])  disagree++;

  } while ( (disagree > 0 || agree < 4) && ++n_tests < max_tests);

  // Calculate the Major and Minor Thrust values
  // -------------------------------------------
  // Note: Now that we have the thrust axis, we determine the thrust value if the various calculations of the thrust axes disagree,
  //       try all and take the maximum, calculate minor and mayor axis. nT_0[n_test] keeps the thrust axis found for all iterations
  //       so here we can calculate the thrust for each iteration, if the axis disagreed between the 4 first, and the thrust will be
  //       the one that is the largest.

  n_tests = 0;

  do {
    double denominator = 0;
    double numerator_t = 0;
    double numerator_m = 0;

    // Calculate the numerator and denominator major and minor Thrusts
    // ...............................................................
    // Note: Make sure the denominator does not diverge

    // Get the vector for each particles to be used in the thrust calculation
    for(std::vector<TLorentzVector>::const_iterator i = iBeg; i != iEnd; ++i) {
      
      const TLorentzVector & itr = *i;
      TLorentzVector v((itr).Px(), (itr).Py(), (itr).Pz(), (itr).E());
      TVector3 c(v.Vect());

      // Only consider transverse directions
      c.SetZ(0);

      // Calculate numerator and denominator of major and minor thrusts
      numerator_t += fabs(c.Dot(nT_0[n_tests]));
      numerator_m += (c.Cross(nT_0[n_tests])).Mag();
      denominator += c.Mag();

    }

    if(denominator < 1e-20) {
      //FPE
      return Variables;
    }

    // Calculate the major and minor thrust
    // ....................................
    // Note: Verify that the major trust is the max, otherwise use the thrust axis found from the next seed,
    //       until max is found unless the first 4 seeds converged to the same axis.

    if (numerator_t / denominator > transvthrust_major) {
      transvthrust_major = numerator_t / denominator;
      transvthrust_minor = numerator_m / denominator;
      transvthrust=nT_0[n_tests];
    }

  } while (disagree > 0 && ++n_tests < max_tests);

  // Fill the map and return variables
  // ---------------------------------

  Variables["TransvThrustMin"] = transvthrust_minor;
  Variables["TransvThrustMaj"] = 1. - transvthrust_major;
  Variables["TransvThrustAxisX"] = transvthrust(0);
  Variables["TransvThrustAxisY"] = transvthrust(1);
  Variables["TransvThrustAxisZ"] = transvthrust(2);

  return Variables;

}
//-------------------------------------------------------------------------------------------------------------



// ============================================================================================================
map<string, double> ANA_utils::CentralGlobalCorr(std::vector<TLorentzVector> input, bool addChglep, std::vector<TruthPart> LeptonBare, bool addNeut, std::vector<float> ETmiss) {
//
// This function is used to calculate corrections that, added to the central version of global observables,
// will make them actually global, despite the limited eta range. The corrected central global quantities must
// be compared to their pseudo-global version. Small differences would mean that the global quantities limited
// to the detector coverage are sufficiently large for not running into resummation issues. On the contrary,
// large differences, especially at low values, would mean that the pseudo-global version is not adequate for
// NNL resummed predictions. The pieces computed here are:
//
//       - QTc: Central qT scalar sum
//       - Central transverse residue
//       - qT-averaged central rapidity
//
// These pieces are of no interest in themselves; only when used to correct centrally restricted global variables
// will these terms make sense. For more details, look at:
//
//           Banfi, Sala,, Zanderighi, "Resummed event shapes at hadron-hadron colliders", JHEP 0408 (2004) 062,
//           arXiv:hep-ph/0407287v3
//
//      Note 1: Only central particles are needed as input.
//
// ============================================================================================================

  // Declare variables internal to this function
  // -------------------------------------------

  // Variables to be returned by the function
  // ........................................

  map<string, double> Variables;
  Variables["QTC"] = -999. * 1000.;
  Variables["Residue"] = -999. * 1000.;
  Variables["qrap_cent"] = -999. * 1000.;

  // Internal variables
  // ..................

  TLorentzVector sumvectqt_c;

  double sumscalqt_c = 0.;
  double sumweightrap_c = 0.;

  double res_c = -1.;
  double eta_c = -1.;

  // Modify the input if chaged leptons from Vector Boson and/or neturinos are to be included
  // ----------------------------------------------------------------------------------------

  // Get all particles obtained from ANA_utils::FindParticles
  // ........................................................

  std::vector<TLorentzVector> input_fin;

  for (int k=0; k<input.size(); k++) {
    input_fin.push_back(input[k]);
  }

  // Include charged leptons from Vector Bosons
  // ..........................................
  // Note: They are stored in a vector of TruthPart, we need to copy them first in a TLorentzVector

  if (addChglep == true) {
    for (int ilep = 0; ilep < LeptonBare.size(); ++ilep) {
      TLorentzVector Lepton_i;
      Lepton_i.SetPtEtaPhiE(LeptonBare[ilep].Pt(), LeptonBare[ilep].Eta(), LeptonBare[ilep].Phi(), LeptonBare[ilep].E());
      input_fin.push_back( Lepton_i);
    }
  }

  // Include neutrino using ETmiss
  // .............................
  // Note: We assume that pz and mass are null. This makes the energy wrong for the 4-vector. Below
  //       only the 3-momentum is used, so it is not a problem.

  if (addNeut == true) {
    TLorentzVector ETmissPart;
    ETmissPart.SetPxPyPzE(ETmiss[0], ETmiss[1], 0., sqrt(ETmiss[0]*ETmiss[0]+ETmiss[1]*ETmiss[1]));
    input_fin.push_back(ETmissPart);
  }

  // Perform the calculations
  // ------------------------

  // Loop over the input particles
  // .............................

  std::vector<TLorentzVector>::const_iterator iBeg = input_fin.begin();
  std::vector<TLorentzVector>::const_iterator iEnd = input_fin.end();
  for (std::vector<TLorentzVector>::const_iterator i = iBeg; i != iEnd; ++i) {

    const TLorentzVector &itr = *i;

    // Calculate Sum qT, Sum vec(qT) and Sum qT*eta_C
    // ..............................................
    // Note: Adding or not the ETmiss to the weighted rapidity is irrelevant because pz is assumed to be 0 for ETmiss, and so eta = 0 and does
    //       not contributed to the weighted sum anyway.

    TLorentzVector v((itr).Px(), (itr).Py(), (itr).Pz(), (itr).E());

    sumvectqt_c += v;
    sumscalqt_c += v.Pt();
    sumweightrap_c += v.Pt()*v.Eta();

  }

  // Calculate the residue and eta_c
  // -------------------------------

  if (sumscalqt_c > 0.) {
    res_c = sumvectqt_c.Pt() / sumscalqt_c;
    eta_c = sumweightrap_c / sumscalqt_c;
  }

  // Fill the map and return variables
  // ---------------------------------

  Variables["QTC"] = sumscalqt_c;
  Variables["Residue"] = res_c;
  Variables["qrap_cent"] = eta_c;

  return Variables;

}
//-------------------------------------------------------------------------------------------------------------



// ============================================================================================================
double ANA_utils::ForwardSuppressionCorr(std::vector<TLorentzVector> input, double QTC, double eta_c) {
//
// This function is used to calculate corrections that, added to the central version of global observables,
// will make them actually global, despite the limited eta range. The approach here is different than the residue
// approach: it consists on suppressing the contribution from forward particles using an exponential growing
// decreasing with with the rapidity separation increase between the forward particles and the qt weighted
// rapidity average for central particles. That requires that all particles are used. In practice this is not
// possible, but we can compare the impact of suppressing from central to 4.9 with suppressing from central to
// 100. Chances are that the impact will be very small. This can be easily tested by simply inputing particles
// limited to 4.9 in pseudo-rapidity rather than taking them all. The function works in all cases. The output
// of the function is:
//
//       - The forward suppression correction factor
//
// Once again, this piece has no interest in itself; only when used to correct centrally restricted global variables
// will this term makes sense. For more details, look at:
//
//           Banfi, Sala,, Zanderighi, "Resummed event shapes at hadron-hadron colliders", JHEP 0408 (2004) 062,
//           arXiv:hep-ph/0407287v3
//
//      Note 1: All final state particles are needed as input except neutrinos.
//
// ============================================================================================================

  // Declare variables internal to this function
  // -------------------------------------------

  // Variables to be returned by the function
  // ........................................

  double Variables;

  // Internal variables
  // ..................

  double num = 0.;

  // Perform the calculations
  // ------------------------

  // Loop over the input particles
  // .............................

  std::vector<TLorentzVector>::const_iterator iBeg = input.begin();
  std::vector<TLorentzVector>::const_iterator iEnd = input.end();
  for (std::vector<TLorentzVector>::const_iterator i = iBeg; i != iEnd; ++i) {

    const TLorentzVector &itr = *i;

    // Calculate Sum qT, Sum vec(qT) and Sum qT*eta_C
    // ..............................................

    TLorentzVector v((itr).Px(), (itr).Py(), (itr).Pz(), (itr).E());

    double exp_i = -999.;
    double eta_i = v.Eta();
    double q_i = v.Pt();
    exp_i = exp(-fabs(eta_i - eta_c));

    num += q_i * exp_i;

  }

  // Calculate the variable and return
  // ---------------------------------

  if (QTC > 0.) {
    Variables = num/QTC;
  }

  return Variables;

}
//-------------------------------------------------------------------------------------------------------------



// ============================================================================================================
map<string, double> ANA_utils::CentralUpAndDownHemisphereQuantities(std::vector<TLorentzVector> input, double QTC, TVector3 ThrustAxis_c) {
//
// Having determined a central transverse thrust axis nT_C, one can separate the central region C into an up part
// C_U consisting of all particles in C with pT.nT_C > 0 and a down part C_D, particles in C with pT.nT_C < 0. We
// can use this to define a normalized squared invariant mass in each hemisphere as:
//
//          rho_C,X = [1 / (QT_C)^2] * [ Sum_i = q_i ]^2, where i is in C_X, and X = either U or D
//
// From these, we can obtain the non-global central sum of masses and heavy-mass as:
//
//            rho_C,S = rho_C,U + rho_C,D     and     rho_C,H = max{rho_C,U; rho_C,D}
//
// These quantities can be globalised by adding the suppression correction or the residue corrections calculated
// in functions above. We can also define jet broadening quantities defined in each hemisphere as:
//
//          B_C,X = [1 / (2*QT_C)] * [Sum_i = sqrt ( (eta_i - eta_C,X)^2 + (phi_i - phi_C,X)^2 ), where i is in C,X.
//
// From these, we can obtain the central total broadenings and the wide broadenings as:
//
//            B_Tot,C =  B_C,U + B_C,D       and       B_Wide,C = max{B_C,U; B_C,D}
//
// Finally, we can obtain the superspherocity as:
//
//                          SuperSphero = min{lambda_2,U; lambda_2,D}
//
// where lambda_2,X is the lowest eigenvalue of the transverse sphericity matrix obtained from the hemishpere X.
//
// Once again, these pieces have no interest in themselves; only when used to correct centrally restricted global
// variables will these term makes sense (except maybe the super-sphericity). For more details, look at:
//
//           Banfi, Sala,, Zanderighi, "Resummed event shapes at hadron-hadron colliders", JHEP 0408 (2004) 062,
//           arXiv:hep-ph/0407287v3
//
//      Note 1: All final state particles are needed as input except neutrinos.
//
//      Note 2: It is assumed in the paper that phi range between 0 and 2*Pi, so we must make sure that when
//              the sums of phi, this is satisfied beforehand.  
//
// ============================================================================================================

  // Declare variables internal to this function
  // -------------------------------------------

  // Variables to be returned by the function
  // ........................................

  map<string, double> Variables;
  Variables["CentralMass"] = -999. * 1000.;
  Variables["HeavyMass"] = -999. * 1000.;
  Variables["TotBroadenings"] = -999. * 1000.;
  Variables["WideBroadenings"] = -999. * 1000.;
  Variables["SuperSpherocity"] = -999. * 1000.;

  // Internal variables
  // ..................

  TLorentzVector Num4vecCU;
  TLorentzVector Num4vecCD;

  double NumEtaCU = 0.;
  double NumEtaCD = 0.;
  double NumRapCU = 0.;
  double NumRapCD = 0.;
  double NumPhiCU = 0.;
  double NumPhiCD = 0.;
  double QTCU = 0.;
  double QTCD = 0.;

  std::vector<double> eta_iCU;
  std::vector<double> eta_iCD;
  std::vector<double> rap_iCU;
  std::vector<double> rap_iCD;
  std::vector<double> phi_iCU;
  std::vector<double> phi_iCD;
  std::vector<double> qT_iCU;
  std::vector<double> qT_iCD;

  double NumBroadCU = 0.;
  double NumBroadCD = 0.;

  double rhoCU = -999.;
  double rhoCD = -999.;
  double rhoCS = -999.;
  double rhoCH = -999.;

  double etaCU = -999.;
  double rapCU = -999.;
  double phiCU = -999.;
  double etaCD = -999.;
  double rapCD = -999.;
  double phiCD = -999.;

  double broadCU = -999.;
  double broadCD = -999.;
  double TotBroad = -999.;
  double WideBroad = -999.;

  TMatrixD TransMomTensorCU(2,2);
  TMatrixD TransMomTensorCD(2,2);

  double SuperSpherocity = -999.;

  double lambda_U2 = -999.;
  double lambda_D2 = -999.;

  // Perform the calculations
  // ------------------------

  // Loop over the input particles
  // .............................

  std::vector<TLorentzVector>::const_iterator iBeg = input.begin();
  std::vector<TLorentzVector>::const_iterator iEnd = input.end();
  for (std::vector<TLorentzVector>::const_iterator i = iBeg; i != iEnd; ++i){

    // Get transverse momentum vector for each input particles
    // .......................................................

    const TLorentzVector &itr = *i;
    TLorentzVector v((itr).Px(), (itr).Py(), (itr).Pz(), (itr).E());
    TVector3 pTi(v.Vect());
    pTi.SetZ(0.);


    // Calculate the sign of the projection of pTi on the thrust axis
    // ..............................................................

    float projection_i = pTi.Dot(ThrustAxis_c);
    bool CU = false;
    bool CD = false;

    if (projection_i > 0.) {CU = true;}
    else if (projection_i < 0.) {CD = true;}

    // Calculate quantites defined in each separate hemisphere
    // .......................................................

    // Note: Need to get  sum of 4 vectors, sum qT, weighted rapidities and weighted azimuthal
    //       angles for each hemisphere

    if (CU) {
      Num4vecCU += v;
      NumEtaCU += (v.Pt())*(v.Eta());
      NumRapCU += (v.Pt())*(v.Rapidity());
      double vphi_cu = v.Phi();
      if (vphi_cu < 0.) vphi_cu = v.Phi()+6.283185307;
      NumPhiCU += (v.Pt())*(vphi_cu);
      QTCU += v.Pt();
      eta_iCU.push_back(v.Eta());
      rap_iCU.push_back(v.Rapidity());
      phi_iCU.push_back(vphi_cu);
      qT_iCU.push_back(v.Pt());
    }

    if (CD) {
      Num4vecCD += v;
      NumEtaCD += (v.Pt())*(v.Eta());
      NumRapCD += (v.Pt())*(v.Rapidity());
      double vphi_cd = v.Phi();
      if (vphi_cd < 0.) vphi_cd = v.Phi()+6.283185307;
      NumPhiCD += (v.Pt())*(vphi_cd);
      QTCD += v.Pt();
      eta_iCD.push_back(v.Eta());
      rap_iCD.push_back(v.Rapidity());
      phi_iCD.push_back(vphi_cd);
      qT_iCD.push_back(v.Pt());
    }

  } // End loop on input particles

  // Calculate centreal sum and heavy masses
  // .......................................

  if (QTC > 0.) {
    rhoCU = ( Num4vecCU.Mag2() ) / (QTC*QTC);
    rhoCD = ( Num4vecCD.Mag2() ) / (QTC*QTC);
  }

  rhoCS = rhoCU + rhoCD;

  if (rhoCU > rhoCD) {rhoCH = rhoCU;}
  else if (rhoCD > rhoCU) {rhoCH = rhoCD;}

  // Calculate jet broadning and supersphero tensor for each hemisphere
  // ...............................................................

  if (QTCU > 0.) {
    etaCU = ( NumEtaCU ) / (QTCU);
    rapCU = ( NumRapCU ) / (QTCU);
    phiCU = ( NumPhiCU ) / (QTCU);
  }

  if (QTCD > 0.) {
    etaCD = ( NumEtaCD ) / (QTCU);
    rapCD = ( NumRapCD ) / (QTCU);
    phiCD = ( NumPhiCD ) / (QTCU);
  }

  for (int ju = 0; ju< eta_iCU.size(); ju++) {
    double delta_eta_ju = eta_iCU[ju] - etaCU;
    double delta_rap_ju = rap_iCU[ju] - rapCU;
    double delta_phi_ju = delta_phi(phi_iCU[ju],phiCU);

    NumBroadCU += qT_iCU[ju]*sqrt( delta_eta_ju*delta_eta_ju + delta_phi_ju*delta_phi_ju );

    TransMomTensorCU(0,0) += qT_iCU[ju]*delta_rap_ju*delta_rap_ju;
    TransMomTensorCU(0,1) += qT_iCU[ju]*delta_rap_ju*delta_phi_ju;
    TransMomTensorCU(1,0) += qT_iCU[ju]*delta_phi_ju*delta_rap_ju;
    TransMomTensorCU(1,1) += qT_iCU[ju]*delta_phi_ju*delta_phi_ju;
  }

  for (int jd = 0; jd< eta_iCD.size(); jd++) {
    double delta_eta_jd = eta_iCD[jd] - etaCD;
    double delta_rap_jd = rap_iCD[jd] - rapCD;
    double delta_phi_jd = delta_phi(phi_iCD[jd], phiCD);

    NumBroadCD += qT_iCD[jd]*sqrt( delta_eta_jd*delta_eta_jd + delta_phi_jd*delta_phi_jd );

    TransMomTensorCD(0,0) += qT_iCD[jd]*delta_rap_jd*delta_rap_jd;
    TransMomTensorCD(0,1) += qT_iCD[jd]*delta_rap_jd*delta_phi_jd;
    TransMomTensorCD(1,0) += qT_iCD[jd]*delta_phi_jd*delta_rap_jd;
    TransMomTensorCD(1,1) += qT_iCD[jd]*delta_phi_jd*delta_phi_jd;
  }

  // Calculate centreal jet broadenings
  // ..................................

  if (QTC > 0.){
    broadCU = ( NumBroadCU ) / (2.*QTC);
    broadCD = ( NumBroadCD ) / (2.*QTC);
  }

  TotBroad = broadCU + broadCD;

  if (broadCU > broadCD) {WideBroad = broadCU;}
  else if (broadCD > broadCU) {WideBroad = broadCD;}

  // Calculate super spherocity
  // ..........................

  if (QTC > 0.) {
    const double inv_QTC = 1. / QTC;

    for (int k1=0; k1<2; k1++) {
      for (int k2=0; k2<2; k2++) {
        TransMomTensorCU(k1,k2) *= inv_QTC;
        TransMomTensorCD(k1,k2) *= inv_QTC;
      }
    }

    TDecompSVD * aTSVD_CU = new TDecompSVD(TransMomTensorCU);
    TVectorD TLambda_CU = aTSVD_CU->GetSig();
    TDecompSVD * aTSVD_CD = new TDecompSVD(TransMomTensorCD);
    TVectorD TLambda_CD = aTSVD_CD->GetSig();

    lambda_U2 = TLambda_CU[1];
    lambda_D2 = TLambda_CD[1];

    delete aTSVD_CU;
    delete aTSVD_CD;
  } // end if QTC>0

  if (lambda_U2 < lambda_D2) {SuperSpherocity = lambda_U2;}
  else if (lambda_D2 < lambda_U2) {SuperSpherocity = lambda_D2;}

  // Return quantities
  // -----------------

  Variables["CentralMass"] = rhoCS;
  Variables["HeavyMass"] = rhoCH;
  Variables["TotBroadenings"] = TotBroad;
  Variables["WideBroadenings"] = WideBroad;
  Variables["SuperSpherocity"] = SuperSpherocity;

  return Variables;

}
//-------------------------------------------------------------------------------------------------------------





//*************************************************************************************************************
//*************************************************************************************************************
//
//    Truth Level Particle Reconstruction Functions:
//
//         -Fill_TruthPart: Fill a TruthParticle object.
//
//         -Get_Tops: Fill a collection of TruthParticle with the top and anti-tops in the event.
//
//         -Get_VectorBosons: Fill a collection of TruthParticle with all prompt vector bosons in an event.
//
//         -Bare_Welectron4Vec: Fill a Collection of TLorentzVector with bare level electrons from W
//
//         -Bare_WelectronTruePart: Fill a Collection of TruthParticle with bare level electrons from W
//
//         -Get_BarePromptLepton: Fill a collection of TruthParticle with bare prompt leptons and another one for neutrinos.
//
//         -TrueJetsReco: Reconstruct jets, and fill a Collection of TruthJets, and deal with b-tagging.
//
//         -Get_BottomQuarks: Fill a collection of TruthParticle with the bottom and anti-bottoms in the event.
//
//         -Get_BottomHadrons: Fill a collection of TruthParticle with the B-hadrons from promt b-quarks.
//
//
//*************************************************************************************************************
//*************************************************************************************************************


// ============================================================================================================
void ANA_utils::Fill_TruthPart(Pythia8::Event event, int index, TruthPart* p_TruthPart) {
//
//  Fill a TruthPart object from the particle "index" in the Pythia event record.
//  
//     Note: The object is filled via the pointer set in the function which calls Fill_TruthPart
//  
// ============================================================================================================

  // Fill a temporary 4-vector from the particle kinematic
  // -----------------------------------------------------

  double partPt  = event[index].pT();
  double partEta = event[index].eta();
  double partPhi = event[index].phi();
  double partM = event[index].m();
    
  TLorentzVector temp_4vec;
  temp_4vec.SetPtEtaPhiM(partPt,partEta,partPhi,partM);

  // Get the internal Pythia information on this particle
  // ----------------------------------------------------
  
  double charge = event[index].charge();
  int status = event[index].status();
  int statusHepMC = event[index].statusHepMC();
  int pdgid = event[index].id();
  int part_index = index;

  // Fill the TruthPart object
  // -------------------------
  
  p_TruthPart->Set(temp_4vec.E(), temp_4vec.Px(), temp_4vec.Py(), temp_4vec.Pz(), partM, partPt, partEta, partPhi, charge, part_index, status, statusHepMC, pdgid);
  
}
// ============================================================================================================



// ============================================================================================================
void ANA_utils::Get_Tops(Pythia8::Event event, std::vector<TruthPart>* p_Top_Coll) {
//
// Find the top and the anti-top in an event using the last entry in the event record before the top decay. Fill
// a TruthPart object with each of the top and anti-top and store them in a collection. Make sure that if there
// is one single top, only that particle will be stored in the collection and that the function won't crash.  
// 
//
// ============================================================================================================
  
  // Find the top particles
  // ----------------------
  // Note: Start from the last entry in the event record and move up the list until a top is found. That is
  //       the last one before it decays.
  
  int top = -1;
  int antitop = -1;

  for (int i = event.size() - 1; i > 0; i--) {
    if (event[i].id() == 6 && top == -1) top = i;
    if (event[i].id() == -6 && antitop == -1) antitop = i;
  }

  // Fill truth particle object and store in the collection for each top
  // -------------------------------------------------------------------
  // Note: Need to define a pointer to the top TruthPart in order for the function Fill_TruthPart to
  //       fill it properly.
  
  TruthPart temp_top;
  TruthPart* p_temp_top;
  TruthPart temp_antitop;
  TruthPart* p_temp_antitop;

  p_temp_top = &temp_top;   
  p_temp_antitop = &temp_antitop; 

  if (top != -1) Fill_TruthPart(event, top, p_temp_top);
  if (antitop != -1) Fill_TruthPart(event, antitop, p_temp_antitop);

  // Store the top particles in the top collection
  // .............................................

  if (top != -1 ) p_Top_Coll->push_back(temp_top);
  if (antitop != -1 )p_Top_Coll->push_back(temp_antitop);
  
}
// ============================================================================================================




// ============================================================================================================
void ANA_utils::Get_VectorBosons(Pythia8::Event event, std::vector<TruthPart>* p_VecBoson_Coll) {
//
// Find the W and Z bosons, fill a TruthPart object for each vector bosons and store them in a collection. Note
// that carbon copies and chains of intermediate steps are also stored in the even record, giving many more W
// or Z in the event record than there is in the event process. The last copy is the one that carry the right
// energy and momentum. To find it, the Pythia function Particle::iBotCopyID() looks for copy of the same
// particle down stream. If there is no further copy, the function returns the index of the particle tested.
// The vector bosons to keep are those for which there is no further copy in the "descendence".
//
//    Note 1: This is a different approach for doing the same thing as what is done to find the top and the
//            anti-top. It has been shown that they yield the same result, although the current strategy seems
//            more robust. Both are kept to provide different examples.
//
//    Note 2: For further comments on how this function works, see Get_Tops.
//
// ============================================================================================================

  // Find the W or Z bosons
  // ----------------------

  std::vector<int> vecboson_index;
  
  for (int i = event.size() - 1; i > 0; i--) {
    if (event[i].idAbs() == 24 || event[i].idAbs() == 23) {
      if (event[i].iBotCopyId()==i) vecboson_index.push_back(i);
    }
  }

  // Fill truth particle object and store in the collection for each vector boson 
  // ----------------------------------------------------------------------------

  if (vecboson_index.size() == 0) return;
  
  for (int i_vb = 0; i_vb < vecboson_index.size(); i_vb++)
    {

      TruthPart temp_vb;
      TruthPart* p_temp_vb;
      
      p_temp_vb = &temp_vb;
      
      Fill_TruthPart(event, vecboson_index[i_vb], p_temp_vb);
    
      // Store the truth particle the vector boson collection
      // ....................................................

      p_VecBoson_Coll->push_back(temp_vb);

    }

} 
// ============================================================================================================



// ============================================================================================================
void ANA_utils::Get_PromptPhotons(Pythia8::Event event, std::vector<TruthPart>* p_PromptPhotons_Coll) {
//
// Find the W and Z bosons, fill a TruthPart object for each vector bosons and store them in a collection. Note
// that carbon copies and chains of intermediate steps are also stored in the even record, giving many more W
// or Z in the event record than there is in the event process. The last copy is the one that carry the right
// energy and momentum. To find it, the Pythia function Particle::iBotCopyID() looks for copy of the same
// particle down stream. If there is no further copy, the function returns the index of the particle tested.
// The vector bosons to keep are those for which there is no further copy in the "descendence".
//
//    Note 1: This is a different approach for doing the same thing as what is done to find the top and the
//            anti-top. It has been shown that they yield the same result, although the current strategy seems
//            more robust. Both are kept to provide different examples.
//
//    Note 2: For further comments on how this function works, see Get_Tops.
//
// ============================================================================================================

  // Find the prompt photons
  // -----------------------

     // Note: To get this we first look at the photon involved in the hard process and follow its carbon copy
     //       to the bottom of the chain. We verified with the daughters history that it is the right particle.
     //       The index of this particle is stored. 

  std::vector<int> promptphoton_index;


  for (int i = 0; i < event.size(); ++i)
    {
      if (event[i].status() == -23 and event[i].id() == 22)
	{
	  promptphoton_index.push_back(event[i].iBotCopy());
	}

    }

  

  
  // Fill truth particle object and store in the collection for each vector boson 
  // ----------------------------------------------------------------------------

     // Note: We verified that this is a stable photon with the momentum coming from the prompt production
     //       with recoiling modifications. 
  

  if (promptphoton_index.size() == 0) return;
  
  for (int i_phot = 0; i_phot < promptphoton_index.size(); i_phot++)
    {

      TruthPart temp_phot;
      TruthPart* p_temp_phot;
      
      p_temp_phot = &temp_phot;
      
      Fill_TruthPart(event, promptphoton_index[i_phot], p_temp_phot);
      
    
      // Store the truth particle the vector boson collection
      // ....................................................

      p_PromptPhotons_Coll->push_back(temp_phot);

    }
  
} 
// ============================================================================================================



// ============================================================================================================
void ANA_utils::Bare_Welectron4Vec(Pythia8::Event event, std::vector<TLorentzVector>* p_Born_Coll) {
//
// Find the electron at Bare level from the decay of a W, and store the 4-vector in the empty collection.
//
//   Note: This function is very simple and works only for events with a single W decaying to an electron. It
//         is kept as a simple example that can be used. Get_BarePromptLepton is a much more general function.
//
// ============================================================================================================

  // Find the W boson
  // ----------------
  
  int idxW = -1;

  for (int i = event.size() - 1; i > 0; i--) {
    if (event[i].idAbs() == 24) {
      idxW = i;
      break;
    }
  }
  if (idxW == -1) {
    cout << "Error: Could not find W" << endl;
    return;
  }

  // Find the bare electron from the W decay
  // ---------------------------------------

  int idxElec = idxW;
  while(true) {
    int daughter = event[idxElec].daughter1();
    if   (daughter == 0) break;
    else                 idxElec = daughter;
  }

  // Check if the last W "descendent" is a stable electron
  // -----------------------------------------------------
  
  if (event[idxElec].idAbs() != 11 || !event[idxElec].isFinal()) {
    cout << "Error: Found incorrect decay product of the W" << endl;
    return;
  }

  // Fill 4-vector
  // -------------
  
  double elecPt  = event[idxElec].pT();
  double elecEta = event[idxElec].eta();
  double elecPhi = event[idxElec].phi();
  double elecM = event[idxElec].m();
  
  TLorentzVector temp_elec;
  temp_elec.SetPtEtaPhiM(elecPt,elecEta,elecPhi,elecM);
    
  // Store the 4 vector in the bare electron collection
  // ..................................................
    
  p_Born_Coll->push_back(temp_elec);

} 
// ============================================================================================================



// ============================================================================================================
void ANA_utils::Bare_WelectronTruePart(Pythia8::Event event, std::vector<TruthPart>* p_BornElec_Coll) {
//
// Find the electron at Bare level from the decay of a W, fill a true particle objects with the information about
// this electron and store it in the empty collection.
//
//    Note: While the function is more sophisticated than Bare_Welectron4Vec, it is still less general than
//          Get_BarePromptLepton which should be the favored function to find bare leptons. This one is
//          nevertheless kept for reference and for testing of other functions.   
//
// ============================================================================================================

  // Find the W or Z bosons
  // ----------------------

  std::vector<int> vecboson_index;
  
  for (int i = event.size() - 1; i > 0; i--) {
    if (event[i].idAbs() == 24 || event[i].idAbs() == 23) {
      if (event[i].iBotCopyId()==i) vecboson_index.push_back(i);
    }
  }
  
  // Find the electron from the W decay
  // ----------------------------------
  
  int idxElec = vecboson_index[0];
  while(true) {
    int daughter = event[idxElec].daughter1();
    if   (daughter == 0) break;
    else                 idxElec = daughter;
  }

  // Check if the last W "descendent" is a stable electron
  // -----------------------------------------------------

  if (event[idxElec].idAbs() != 11 || !event[idxElec].isFinal()) {
    cout << "Error: Found incorrect decay product of the W" << endl;
    return;
  }
      
  // Fill truth particle object and store in the collection
  // ------------------------------------------------------

  double elecPt  = event[idxElec].pT();
  double elecEta = event[idxElec].eta();
  double elecPhi = event[idxElec].phi();
  double elecM = event[idxElec].m();
    
  TLorentzVector temp_elec;
  temp_elec.SetPtEtaPhiM(elecPt,elecEta,elecPhi,elecM);

  double charge = event[idxElec].charge();
  int status = event[idxElec].status();
  int statusHepMC = event[idxElec].statusHepMC();
  int pdgid = event[idxElec].id();
  int part_index = idxElec;

  TruthPart temp_truth;
  
  temp_truth.Set(temp_elec.E(), temp_elec.Px(), temp_elec.Py(), temp_elec.Pz(), elecM, elecPt, elecEta, elecPhi, charge, part_index, status, statusHepMC, pdgid);

  // Store the 4 vector in the bare electron collection
  // ..................................................

  p_BornElec_Coll->push_back(temp_truth);

} 
// ============================================================================================================



// ============================================================================================================
void ANA_utils::Get_BarePromptLepton(Pythia8::Event event, std::vector<int> vecboson_index, std::vector<TruthPart>* p_PromptBareLept_Coll, std::vector<TruthPart>* p_ConversionLept_Coll, std::vector<TruthPart>* p_PhotonFSR_Coll, std::vector<TruthPart>* p_Neutrino_Coll) {
//
// Find all leptons at Bare level coming from the decay of a W or a Z. Fill a true particle objects collection
// with the information about these leptons. This function is valid for both electrons and muons, and for events
// with single or multiple vector bosons (e.g. diboson, ttbar, etc.). The StableDaughterList vector keeps all
// the stable particles descending from a vector boson. It therefore has all FSR photons for dressing of leptons
// or jet correction if needed. Neutrinos are also stored.
//
//    Note 1: Because of the QED emission of the charged leptons, and of the carbon copies of the particles stored
//            in the Pythia event record, the immediate daughters of a vector bosons immediately before it decays
//            is only rarely the stable particles to be used at bare level. It is therefore important to follow
//            the entire decay chain until such stable particles are reached. The number of generation to cover
//            varies for different vector bosons in one event, and even for different leptons. For example, if one
//            Z lepton emits two photons, but not the other, the second lepton lepton will be stable from the first
//            generation of the Z daugthers, while the first one will only be stable at the third generation, and
//            the two photons will respectively be stable at the second and third generations. To this end, a
//            complete enquiry of each daughters at each generation is required. After a stress test of 1000 events,
//            there was no need to go further than the 6th generation. The code therefore stops at the 7th generation.
//  
//    Note 2: It happens in a few events that there are two more charged leptons than expected. This is the result
//            of a photon conversion into a pair of lepton-anti-lepton. These are kept as stable leptons assigned
//            to the decay of a vector boson. A cut on the number of charged leptons expected in the event will
//            take care of these. It happens both in the predictions and in the data.
//  
// ============================================================================================================
    
  // Find the stable lepton from the vector boson decays
  // ---------------------------------------------------

  std::vector<int> lepton_index;
  std::vector<int> fsrphoton_index;
  std::vector<int> neutrino_index;
  
  // Loop over all vector bosons found in the event
  // ..............................................
  // Note: This require using the function Get_VectorBosons, and passing the index of each vector boson as
  //       input to this function
  
  for (int i=0; i< vecboson_index.size(); i++) {
    int idxVb = vecboson_index[i];

    // Define vectors to keep the indices of each particles in each generations
    // ........................................................................
    
    std::vector<int> StableDaughterList;

    std::vector<int> gen1;
    std::vector<int> gen2;
    std::vector<int> gen3;
    std::vector<int> gen4;
    std::vector<int> gen5;
    std::vector<int> gen6;
    std::vector<int> gen7;

    // The first generation is straighforwardly the daughters of the vector boson
    // ..........................................................................
    
    gen1.push_back(event[idxVb].daughter1());
    gen1.push_back(event[idxVb].daughter2());

    // If the vector bosons have not been forced to be stable, checked if the daughters are stable
    // ...........................................................................................
    // Note: Any unstable daughter will have their own daughters populating the next generation
    
    if (gen1.size() != 0) {

      if ( event[gen1[0]].isFinal()  )  StableDaughterList.push_back(gen1[0]);
      else {
        for (int d1 = 0; d1 <  (event[gen1[0]].daughterList()).size(); d1++) {
    	  gen2.push_back((event[gen1[0]].daughterList())[d1]);
    	}
      } 
      if ( event[gen1[1]].isFinal()  )  StableDaughterList.push_back(gen1[1]);
      else {
        for (int d2 = 0; d2 <  (event[gen1[1]].daughterList()).size(); d2++){
    	  gen2.push_back((event[gen1[1]].daughterList())[d2]);
    	}
      }

      // Check if the particles in the second generation are stable, if not moved to the next
      // ....................................................................................  

      if (gen2.size() != 0) {
        for (int i_gen2 = 0; i_gen2 < gen2.size(); i_gen2++) {
    	  if ( event[gen2[i_gen2]].isFinal()  ) StableDaughterList.push_back(gen2[i_gen2]);
    	  else {
    	    for (int d3 = 0; d3 <  (event[gen2[i_gen2]].daughterList()).size(); d3++) {
              gen3.push_back((event[gen2[i_gen2]].daughterList())[d3]);
            }
    	  }
    	}
          
        // Check if the particles in the third generation are stable, if not moved to the next
        // ....................................................................................  

        if (gen3.size() != 0) {
    	  for (int i_gen3 = 0; i_gen3 < gen3.size(); i_gen3++) {
    	    if ( event[gen3[i_gen3]].isFinal()  ) StableDaughterList.push_back(gen3[i_gen3]);
    	    else {
              for (int d4 = 0; d4 <  (event[gen3[i_gen3]].daughterList()).size(); d4++) {
    		gen4.push_back((event[gen3[i_gen3]].daughterList())[d4]);
    	      }
	    }
    	  }

          // Check if the particles in the forth generation are stable, if not moved to the next
          // ....................................................................................  
    	  
    	  if (gen4.size() != 0) {
    	    for (int i_gen4 = 0; i_gen4 < gen4.size(); i_gen4++) {
    	      if ( event[gen4[i_gen4]].isFinal()  ) StableDaughterList.push_back(gen4[i_gen4]);
    	      else {
    	        for (int d5 = 0; d5 <  (event[gen4[i_gen4]].daughterList()).size(); d5++) {
    	          gen5.push_back((event[gen4[i_gen4]].daughterList())[d5]);
		}
    	      }
    	    }

            // Check if the particles in the fith generation are stable, if not moved to the next
            // ....................................................................................  
    	      
    	    if (gen5.size() != 0) {
	      for (int i_gen5 = 0; i_gen5 < gen5.size(); i_gen5++) {
                if ( event[gen5[i_gen5]].isFinal()  ) StableDaughterList.push_back(gen5[i_gen5]);
    	        else {
    	          for (int d6 = 0; d6 <  (event[gen5[i_gen5]].daughterList()).size(); d6++) {
    	            gen6.push_back((event[gen5[i_gen5]].daughterList())[d6]);
    	 	  }
                }
    	      }

              // Check if the particles in the sixth generation are stable, if not moved to the next
              // ....................................................................................  

              if (gen6.size() != 0) {
                for (int i_gen6 = 0; i_gen6 < gen6.size(); i_gen6++) {
 	          if ( event[gen6[i_gen6]].isFinal()  ) StableDaughterList.push_back(gen6[i_gen6]);
    		  else {
    		    for (int d7 = 0; d7 <  (event[gen6[i_gen6]].daughterList()).size(); d7++) {
		      gen7.push_back((event[gen6[i_gen6]].daughterList())[d7]);
		    }
                  }
	        }
              } //end gen6!=0
    	    } //end gen5!=0
    	  } // end if gen4!=0
    	} // end if gen3!=0
      } // end if gen2!=0
    } // end if gen1!=0

    // Of all the stable particles, find those that are leptons (electrons or muons)
    // .............................................................................
    
    for (int stab_i = 0; stab_i < StableDaughterList.size(); stab_i++) {
      int temp_id = StableDaughterList[stab_i];
      if ( event[temp_id].idAbs() == 11 || event[temp_id].idAbs() == 13 ) lepton_index.push_back(StableDaughterList[stab_i]);
      if ( event[temp_id].idAbs() == 12 || event[temp_id].idAbs() == 14  || event[temp_id].idAbs() == 16 ) neutrino_index.push_back(StableDaughterList[stab_i]);
      if ( event[temp_id].id() == 22) fsrphoton_index.push_back(StableDaughterList[stab_i]);
    }

  } // end loop over vector bosons

  // Fill truth particle object and store in the collection for each stable charged lepton
  // -------------------------------------------------------------------------------------

  if (lepton_index.size() != 0) {
    for (int i_lep = 0; i_lep < lepton_index.size(); i_lep++) {

      TruthPart temp_lep;
      TruthPart* p_temp_lep;

      p_temp_lep = &temp_lep;

      Fill_TruthPart(event, lepton_index[i_lep], p_temp_lep);
    
      // Store the 4 vector in the bare electron collection
      // ..................................................

      p_PromptBareLept_Coll->push_back(temp_lep);

    }
  }

  // Fill truth particle object and store in the collection for each FSR photon
  // --------------------------------------------------------------------------

  if (fsrphoton_index.size() != 0) {
    for (int i_fsr = 0; i_fsr < fsrphoton_index.size(); i_fsr++) {

      TruthPart temp_fsr;
      TruthPart* p_temp_fsr;

      p_temp_fsr = &temp_fsr;

      Fill_TruthPart(event, fsrphoton_index[i_fsr], p_temp_fsr);

      // Store the 4 vector in the bare electron collection
      // ..................................................

      p_PhotonFSR_Coll->push_back(temp_fsr);

    }
  }
  
  
  // Fill truth particle object and store in the collection for each neutrino
  // ------------------------------------------------------------------------

  if (neutrino_index.size() == 0) return;
  
  for (int i_neu = 0; i_neu < neutrino_index.size(); i_neu++){

    TruthPart temp_neu;
    TruthPart* p_temp_neu;

    p_temp_neu = &temp_neu;

    Fill_TruthPart(event, neutrino_index[i_neu], p_temp_neu);

    // Store the 4 vector in the bare electron collection
    // ..................................................

    p_Neutrino_Coll->push_back(temp_neu);

  }

}
// ============================================================================================================


// ============================================================================================================
void ANA_utils::Get_DressPromptLepton(std::vector<TruthPart> LeptonBare_Coll, std::vector<TruthPart> PhotonFSR_Coll, std::vector<TruthPart>* p_DressLep_Coll, std::vector<TruthPart> *p_PhotonDress_Coll) {
  //
  // Use the bare leptons and the FSR photons collections, already filled, and fill the dress lepton collection
  // and the dressing photon collection. The dressing consists in finding all the photons in a cone of 0.1 to
  // the bare lepton and to sum their 4-vectors to the 4-vector of the bare lepton. Make sure to not double count
  // the dressing in case there are more than one bare leptons. To do this, loop first on the FSR leptons (they all
  // come from prompt leptons), find the closest lepton, and if it is closer than dR=0.1, merge the 4-vector. This
  // way, any FSR photon won't be merged more than once, while multiple merging can happen for each prompt lepton
  // as it should. Keep only the photons that were used in the dressing and store them in the PhotonDress collection.
  // They will be used later to calculate dress jets. However, there might be more than one photon per bare lepton,
  // so the above strategy for matching one photon to its bare lepton needs to already be included in a loop over
  // lepton. The overall logic is therefore the following:
  //
  //         1- Loop over bare lepton. Each bare lepton initialize a dress lepton
  //         2- Inside the loop over leptons, loop over FSR photons given as input
  //         3- Inside the loop on FSR, loop over bare leptons again. This time to find to which lepton the
  //            FSR photon must be macthed
  //         4- If that lepton is macthed with dR<0.1 and if that lepton is the same as the initializing one,
  //            then add the 4-vector of the FSR to the bare to dress it, and keep the TruthPart corresponding
  //            to that FSR photon into the DressCollection. It will be a subset of the FSR collection, because
  //            each FSR are matched to one and only one bare lepton, and the dressing happens only if that
  //            lepton that matches this frs photon happens to be the one of the first loop.
  //         5- After the loop on FSR photon, the dressing will be over so fill the dress lepton collection with
  //            the 4-vector of the dressed lepton and the particle info of the original bare lepton.     
  //
  // ============================================================================================================  
    
  
    // Loop over bare lepton and initialize the dress lepton as the bare lepton
    // .........................................................................
  
    for (int i_lep = 0; i_lep < LeptonBare_Coll.size(); i_lep++)
      {
        double ptlep_i = (LeptonBare_Coll[i_lep]).Pt();
        double etalep_i = (LeptonBare_Coll[i_lep]).Eta();
        double philep_i = (LeptonBare_Coll[i_lep]).Phi();
        if (philep_i < 0.) philep_i = (LeptonBare_Coll[i_lep]).Phi() + 6.283185307;	      
        double elep_i = (LeptonBare_Coll[i_lep]).E();
  
        double charge = (LeptonBare_Coll[i_lep]).Charge();
        int status = (LeptonBare_Coll[i_lep]).Status();
        int statusHepMC = (LeptonBare_Coll[i_lep]).StatusHepMC();
        int pdgid = (LeptonBare_Coll[i_lep]).Pdgid();
        int part_index = (LeptonBare_Coll[i_lep]).Index();
  
        TLorentzVector barelep_i;
        barelep_i.SetPtEtaPhiE(ptlep_i,etalep_i,philep_i,elep_i);
  
        TLorentzVector dresslep_i;
        dresslep_i = barelep_i;
  
  
  
    // Loop over FSR photons and get their 4-vector
    // ............................................
  
        for (int i_fsr = 0; i_fsr < PhotonFSR_Coll.size(); i_fsr++)
    {
      double ptfsr_i = (PhotonFSR_Coll[i_fsr]).Pt();
      double etafsr_i = (PhotonFSR_Coll[i_fsr]).Eta();
      double phifsr_i = (PhotonFSR_Coll[i_fsr]).Phi();
      if (phifsr_i < 0.) phifsr_i = (PhotonFSR_Coll[i_fsr]).Phi() + 6.283185307;	      
      double efsr_i = (PhotonFSR_Coll[i_fsr]).E();
      
      TLorentzVector fsr;
      fsr.SetPtEtaPhiE(ptfsr_i,etafsr_i,phifsr_i,efsr_i);
  
  
    // Find to which lepton this FSR photon is the closest for best match
    // ..................................................................
      
      float dphimin = 10;
      int matchindex = -999;
        
  
      for (int j_lep = 0; j_lep < LeptonBare_Coll.size(); j_lep++)
        {
          double etalep_j = (LeptonBare_Coll[j_lep]).Eta();
          double philep_j = (LeptonBare_Coll[j_lep]).Phi();
          if (philep_j < 0.) philep_j = (LeptonBare_Coll[j_lep]).Phi() + 6.283185307;	      
  
          double dphi_j = delta_phi(phifsr_i,philep_j);
          double deta_j = etafsr_i - etalep_j;
          double dr_j = sqrt(dphi_j*dphi_j + deta_j*deta_j);
          
          if (dr_j<dphimin)
      {
        dphimin = dr_j;
        matchindex=j_lep;
      }
      
        }
  
  
      
    // Dress the charged leptons with matched fsr photon if within dR of 0.1 and Strore the dressed photon particles
    // .............................................................................................................
  
        
      if (dphimin<0.1 && matchindex == i_lep)
        {
          dresslep_i += fsr;
  
          p_PhotonDress_Coll->push_back(PhotonFSR_Coll[i_fsr]);
   
        }
      
    } // end loop on fsr
  
  
    // Fill the particle info of the dressed lepton and store it in the collection
    // ...........................................................................
  
        
        TruthPart temp_dressedlep;
    
        temp_dressedlep.Set(dresslep_i.E(), dresslep_i.Px(), dresslep_i.Py(), dresslep_i.Pz(), dresslep_i.M(), dresslep_i.Pt(), dresslep_i.Eta(), dresslep_i.Phi(), charge, part_index, status, statusHepMC, pdgid);
  
        p_DressLep_Coll->push_back(temp_dressedlep);
  
      } // end loop on bare lepton
  
  
  
  }
  // ============================================================================================================
  




// ============================================================================================================
void ANA_utils::Get_BornPromptLepton(Pythia8::Event event, std::vector<int> vecbosindex, std::vector<TruthPart>*  p_LeptonBorn_Coll) {
//
// Find all leptons at Born level coming from the decay of a W or a Z. Fill a true particle objects collection
// with the information about these leptons. This function is valid for both electrons and muons, and for events
// with single or multiple vector bosons (e.g. diboson, ttbar, etc.). The FSR photons are already stored from
// the function Get_BarePromptLepton.
//
// ============================================================================================================

  // Find the Born leptons from the vector boson decays
  // --------------------------------------------------

  std::vector<int> lepton_index;

  // Loop over all vector bosons found in the event
  // ..............................................
  // Note: This require using the function Get_VectorBosons, and passing the index of each vector boson as
  //       input to this function

  for (int i=0; i< vecbosindex.size(); i++)
    {
      int idxVb = vecbosindex[i];

      // Keep the daugthers of the vector bosons if they are Born charged leptons
      // ........................................................................
      // Note: The Born leptons are the immediate daugthers of the last vector boson in the event record,
      //       i.e. the particles input to this function.

      int idxDaugther1 = event[idxVb].daughter1();
      int idxDaugther2 = event[idxVb].daughter2();

      if ( event[idxDaugther1].idAbs() == 11 ||  event[idxDaugther1].idAbs() == 13) lepton_index.push_back(idxDaugther1);
      if ( event[idxDaugther2].idAbs() == 11 ||  event[idxDaugther2].idAbs() == 13) lepton_index.push_back(idxDaugther2);
    }

  // Fill truth particle object and store in the collection for each Born charged lepton
  // -----------------------------------------------------------------------------------

  if (lepton_index.size() != 0)
    {
      for (int i_lep = 0; i_lep < lepton_index.size(); i_lep++)
        {

          TruthPart temp_lep;
          TruthPart* p_temp_lep;

          p_temp_lep = &temp_lep;

          Fill_TruthPart(event, lepton_index[i_lep], p_temp_lep);

          // Store the 4 vector in the bare electron collection
          // ..................................................

          p_LeptonBorn_Coll->push_back(temp_lep);

        }
    }

}
// ============================================================================================================




// ============================================================================================================
void ANA_utils::TrueJetsReco(Pythia8::Event event, std::vector<int> partskipped, std::vector<TruthJets>* p_TruthJets_Coll, float ptcut, bool doKtSplitting, std::vector<double>* p_KtSplittingScale1, std::vector<double>* p_KtSplittingScale2, std::vector<double>* p_KtSplittingScale3, double R) {
//
// Reconstruct the truth-level jets using FastJets, sort them in pT, fill a TruthJets objects with the relevant
// variables for each jet, and store it in the empty collection. A b-tagger is also implemented in this function.
// To this end, first, the b-quarks or B-hadrons are found, then they are used to create ghost particles. These
// ghosts are clustered with the jets, but they do not impact jet reconstruction. In a last step, the constituents
// of each jets are investigated to see if they have a ghost or not. If they do, they will be b-tagged.
//
//    Note 1: A dR b-quark/b-hadron assignment is also implemented for comparison, but the ghost approach is prefered.
//
//    Note 2: The "partskipped" vectors given the particle event index of the stable particles not to be clustered
//            in the jets, such as final state leptons coming from W and Z for example.
//
//    The code contains further notes below.   
//
// ============================================================================================================

  // Find b-quarks and b-hadrons in the event
  // ----------------------------------------

     // Note: There is a very strong correlation between the kinematic of the b-quark and the kinematic of the b-hadron
     //       following the b-quark hadronization. When printing the pT and eta of the b-quarks and the b-hadrons, one
     //       would see that these quantities have very similar values for b-quarks and b-hadrons. The correlation is
     //       however much weaker with the bjets.

  bool find_bquarks_bhadrons = true;
  if(R==1.0) find_bquarks_bhadrons = false;
  
  std::vector<TruthPart> BottomQuark_Coll;
  std::vector<TruthPart>* p_BottomQuark_Coll = &BottomQuark_Coll;

  if(find_bquarks_bhadrons) Get_BottomQuarks(event, p_BottomQuark_Coll);

  std::vector<TruthPart> BottomHadron_Coll;
  std::vector<TruthPart>* p_BottomHadron_Coll = &BottomHadron_Coll;

  if(find_bquarks_bhadrons) Get_BottomHadrons(event, BottomQuark_Coll, p_BottomHadron_Coll);

  
  // Create ghost particles with each b-quark or b-hadron
  // ----------------------------------------------------

     // Note: A ghost particle is a particle that exactly has the direction of a given particle, but has its energy to almost 0.
     //       Such ghost particle are useful for b-quarks or b-hadrons association to jets, because the ghost will be clustered
     //       in the jet, but the energy and direction of the jet will not be affected at all by the ghost. By looking to which
     //       jet the ghost has been clustered, we'll know to which jet the b-quark or b-hadron can be assigned.

  std::vector<TLorentzVector> Ghost_Bquarks;
  std::vector<TLorentzVector> Ghost_Bhadrons;

  
     // Create ghosts from b-quarks
     // ...........................
      
  if(find_bquarks_bhadrons)
    {
  
      for (int i_bq = 0; i_bq<BottomQuark_Coll.size(); i_bq++)
	{
	  TLorentzVector ghost_bqrk;
	  ghost_bqrk.SetPtEtaPhiM( ((BottomQuark_Coll[i_bq]).Pt()/1000000.), (BottomQuark_Coll[i_bq]).Eta(), (BottomQuark_Coll[i_bq]).Phi(), ((BottomQuark_Coll[i_bq]).M()/1000000.));
	  
	  Ghost_Bquarks.push_back(ghost_bqrk);
	}


     // Create ghosts from b-hadrons
     // ............................
    
      for (int i_bh = 0; i_bh<BottomHadron_Coll.size(); i_bh++)
	{
	  TLorentzVector ghost_bhad;
	  ghost_bhad.SetPtEtaPhiM( ((BottomHadron_Coll[i_bh]).Pt()/1000000.), (BottomHadron_Coll[i_bh]).Eta(), (BottomHadron_Coll[i_bh]).Phi(), ((BottomHadron_Coll[i_bh]).M()/1000000.));
        
	  Ghost_Bhadrons.push_back(ghost_bhad);
	}

    }


  // Build Inclusive Jets
  // --------------------

     // FastJet Initialization
     // ......................

  fastjet::Strategy             strategy = fastjet::Best;
  fastjet::RecombinationScheme  recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition       *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, R, recombScheme, strategy);
  

     // Fastjet input
     // .............
  
  std::vector <fastjet::PseudoJet> fjInputs;

  
     // Reset Fastjet input
     // ...................
  
  fjInputs.resize(0);
  int index=0;

  
     // Loop over event record and specify which particles to pass to FastJet
     // .....................................................................
  
  for (int i = 0; i < event.size(); ++i)
    {

      // Final state only
      if (!event[i].isFinal())        continue;


      // No neutrinos
      if (event[i].idAbs() == 12 || event[i].idAbs() == 14 || event[i].idAbs() == 16)     continue;


      // Only |eta| < 4.9
      if (fabs(event[i].eta()) > 4.9) continue;

      // Do not include the stable particles listed in the partskipped vector
         // Note: This is used to not include the decay product of the W, Z, H and top for example

      bool skip_i = false;
      for (int j = 0; j < partskipped.size(); j++)
	{
	  if (partskipped[j] == i) skip_i = true;
	}
      if (skip_i == true) continue;

      
     // Store as input to Fastjet and set a unique identifier for each input particle
     // .............................................................................

      fastjet::PseudoJet particle( event[i].px(),event[i].py(), event[i].pz(), event[i].e() );

      particle.set_user_index(index);
      fjInputs.push_back( particle );
      index++;
    }

  
     // Add the ghost to the input particles
     // ....................................

        // Note: In contrary to the stable particles constituting the jets, the index value used for
        //       ghost particles is negative, with 0 < ghost_index < -10 for b-quark ghosts, and
        //       with -9 < ghost_index < -100 for b-hadron ghosts.
  
  if(find_bquarks_bhadrons)
    {

      int ghost_index = -1;
      for (int i_gbq = 0; i_gbq<Ghost_Bquarks.size(); i_gbq++)
	{
	  fastjet::PseudoJet ghost_particle( (Ghost_Bquarks[i_gbq]).Px(), (Ghost_Bquarks[i_gbq]).Py(), (Ghost_Bquarks[i_gbq]).Pz(), (Ghost_Bquarks[i_gbq]).E() );
	  
	  ghost_particle.set_user_index(ghost_index);
	  fjInputs.push_back( ghost_particle );
	  ghost_index--;
	}

      for (int i_gbh = 0; i_gbh<Ghost_Bhadrons.size(); i_gbh++)
	{
	  fastjet::PseudoJet ghost_particle( (Ghost_Bhadrons[i_gbh]).Px(), (Ghost_Bhadrons[i_gbh]).Py(), (Ghost_Bhadrons[i_gbh]).Pz(), (Ghost_Bhadrons[i_gbh]).E() );
	  
	  ghost_index = ghost_index - 10;
	  ghost_particle.set_user_index(ghost_index);
	  fjInputs.push_back( ghost_particle );
	}

    }


     // Check that there is some input to Fastjet
     // .........................................
  
  if (fjInputs.size() == 0)
    {
      cout << "Error: event with no final state particles" << endl;
      return;
    }
  

     // Run Fastjet algorithm
     // .....................
  
  vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
  fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

  
     // Extract inclusive jets sorted by pT (note minimum pT of ptcut GeV)
     // .................................................................
  
  inclusiveJets = clustSeq.inclusive_jets(ptcut);
  sortedJets    = sorted_by_pt(inclusiveJets);



  // Fill the jet object and store it in the true jet collection
  // -----------------------------------------------------------

  std::vector<double> temp_scale0, temp_scale1, temp_scale2;
  
  fastjet::JetDefinition kTjetDef;
  if(doKtSplitting) kTjetDef = fastjet::JetDefinition(fastjet::kt_algorithm, R, recombScheme, strategy);

  
     // Loop over all jets produced by Fastjet
     // ......................................

  for (int i_jet = 0; i_jet < sortedJets.size(); i_jet++)
    {

     // Get the constituents of the jet
     // ...............................

      std::vector<fastjet::PseudoJet> constituents = (sortedJets[i_jet]).constituents();

     // Get kT splitting scales
     // ...............................

      if(doKtSplitting)
	{
	  if(constituents.size() == 0)
	    {
	      p_KtSplittingScale1->push_back(0);
	      p_KtSplittingScale2->push_back(0);
	      p_KtSplittingScale3->push_back(0);
	    }
	  else
	    {
	      // get kT splitting scales
	      fastjet::ClusterSequence kTclustSeq(constituents, kTjetDef);
	      fastjet::PseudoJet kTjet = fastjet::sorted_by_pt(kTclustSeq.inclusive_jets()).front();
	      p_KtSplittingScale1->push_back(R*sqrt(kTclustSeq.exclusive_subdmerge(kTjet,1)));
	      p_KtSplittingScale2->push_back(R*sqrt(kTclustSeq.exclusive_subdmerge(kTjet,2)));
	      p_KtSplittingScale3->push_back(R*sqrt(kTclustSeq.exclusive_subdmerge(kTjet,3)));
	    }
	} // End if doKtSplitting
      
            
      // Check if the jet is a B-jet
      // ...........................

          // Note: We have two methods for btagging:
          //
          //            1- Assign a b-quark or a b-hadron to a jet if its angular distance to this jet is smaller than 0.4
          //
          //            2- Include a b-quark or b-hadron ghost (particle pointing in the same direction but with almost no energy and mass)
          //               in the jet clustering and determine in which jets are the ghosts.
          //
          //       These methods are only approximative because the b-quarks hadronize and the b-hadrons decay, and the product could be far
          //       from each others, have low pT or end up in different jets. It is frequent that a b-quark or a b-hadron doesn't get assigned
          //       to any jet. If the jet threshold is low or the jet size is large, the b-quark or b-hadron can be assigned to many jets, all
          //       close in distance compared to the jet size paramater. These problems do not affect the ghost approach which, by construction
          //       always assign each quark or each hadron to a single jet. It is however possible that two b-quarks or b-hadrons get assigned
          //       to the same jet. It is more frequent for larger jet size parameter (about <0.1% for jets of R=0.1, 0.5% for jets of R=0.4 and
          //       about 5% for jets of R=0.8). Also  the assigned jet could be the ghost itself (so no real jet assignment), or it could be
          //       assign to a jet that does not contain the bulk of the original b-quark or b-hadron energy.

      
          // dR Approach
          // .  .  .  .

      bool bqtag = false;
      bool bhtag = false;

      for (auto bottom : BottomQuark_Coll)
	{
	  if (sqrt(pow(bottom.Eta() - (sortedJets[i_jet]).eta(),2) + 
		   pow(delta_phi(bottom.Phi(),(sortedJets[i_jet]).phi()),2)) <= 0.4)
	    {
	      bqtag = true;
	    }
	}
      
      for (auto bottom : BottomHadron_Coll)
	{
	  if (sqrt(pow(bottom.Eta() - (sortedJets[i_jet]).eta(),2) + 
		   pow(delta_phi(bottom.Phi(),(sortedJets[i_jet]).phi()),2)) <= 0.4)
	    {
	      bhtag = true;
	    }
	}

      
         // Ghost Approach
         // .  .  .  .  .
      
      bool bq_ghosttag = false;
      bool bh_ghosttag = false;
      int ghost_mult = 0;
      
      if(find_bquarks_bhadrons)
	{
	  for (int i_const = 0; i_const<constituents.size(); i_const++)
	    {
	      if ( (constituents[i_const]).user_index() < 0 && (constituents[i_const]).user_index() > -10 )
		{
		  bq_ghosttag = true;
		  ghost_mult++;
		}

	      if ( (constituents[i_const]).user_index() < -9  )
		{
		  bh_ghosttag = true;
		  ghost_mult++;
		}
	    }
	}
      
      int npart_injet = constituents.size() - ghost_mult;

      
      // Fill a TruthJets object
      // .......................
      
      TruthJets temp_truthjet;
      temp_truthjet.Set((sortedJets[i_jet]).e(),   (sortedJets[i_jet]).px(),  (sortedJets[i_jet]).py(), 
		        (sortedJets[i_jet]).pz(),  (sortedJets[i_jet]).m(),   (sortedJets[i_jet]).pt(), 
		        (sortedJets[i_jet]).eta(), (sortedJets[i_jet]).rap(), (sortedJets[i_jet]).phi(), 
		         bq_ghosttag, bh_ghosttag, npart_injet);

      
      // Put in the jet collection
      // .........................

      p_TruthJets_Coll->push_back(temp_truthjet);
      
  }

  // Close objects
  // -------------
  
  delete jetDef;
  
}
// ============================================================================================================




// ============================================================================================================
void ANA_utils::PartonJetsReco(Pythia8::Event event, Pythia8::Event partonevent, std::vector<TruthJets>* p_PartonJets_Coll, float ptcut) {
//
// Reconstruct the parton-level jets using FastJets, sort them in pT, fill a TruthJets objects with the relevant
// variables for each jet, and store it in the empty collection. A b-tagger is also implemented in this function.
// To this end, the entire event record is used to find the b-quarks after all QCD radiation. It is then used to
// create ghost particles. Then the parton-only event record is used to build the jets. The ghosts are also included
// in this clustering without impacting the final parton-level jets, just before hadronization. As for the Truth jets
// defined in TruthJetsReco, in the last step the constituents of each jets are investigated to see if they have a
// ghost or not. If they do, they will be b-tagged.
//
//
//    Note: The "partskipped" vectors used in the TruthJetsReco function is useless here because the event record
//          is not the same. It is however easier to get the lepton to ignore here because they come from the
//          hard interaction and therefore have a status value between 21 and 29. We can simply ignore them. 
//
//    More explanations can be found in the function TrueJetsReco. 
//
// ============================================================================================================
  
  // Find b-quarks in the total event
  // --------------------------------

  std::vector<TruthPart> BottomQuark_Coll;
  std::vector<TruthPart>* p_BottomQuark_Coll = &BottomQuark_Coll;

  Get_BottomQuarks(event, p_BottomQuark_Coll);

  
  // Create ghost particles with each b-quark
  // ----------------------------------------

  std::vector<TLorentzVector> Ghost_Bquarks;
  
  for (int i_bq = 0; i_bq<BottomQuark_Coll.size(); i_bq++)
    {
      TLorentzVector ghost_bqrk;
      ghost_bqrk.SetPtEtaPhiM( ((BottomQuark_Coll[i_bq]).Pt()/1000000.), (BottomQuark_Coll[i_bq]).Eta(), (BottomQuark_Coll[i_bq]).Phi(), ((BottomQuark_Coll[i_bq]).M()/1000000.));
      
      Ghost_Bquarks.push_back(ghost_bqrk);
    }

  
  // Build Inclusive Jets with FastJet
  // ---------------------------------

     // FastJet Initialization
     // ......................

  double Rparam = 0.4;
  fastjet::Strategy             strategy = fastjet::Best;
  fastjet::RecombinationScheme  recombScheme = fastjet::E_scheme;
  fastjet::JetDefinition       *jetDef = NULL;
  jetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, Rparam, recombScheme, strategy);

  
     // Fastjet input
     // .............
  
  std::vector <fastjet::PseudoJet> fjInputs;

  
     // Reset Fastjet input
     // ...................
  
  fjInputs.resize(0);
  int index=0;

  
     // Loop over the parton-event record and set the fastjet input particles
     // .....................................................................

        // Note: There is no need to remove unwanted particles because they have already be removed from the parton event record
  
  for (int i = 0; i < partonevent.size(); ++i)
    {

     // Store as input to Fastjet and set a unique identifier for each input particle
     // .............................................................................

      fastjet::PseudoJet particle( partonevent[i].px(), partonevent[i].py(), partonevent[i].pz(), partonevent[i].e() );

      particle.set_user_index(index);
      fjInputs.push_back( particle );
      index++;
    }

     // Add the ghost to the input particles
     // ....................................

        // Note: In contrary to the stable particles constituting the jets, the index value used for
        //       ghost particles is negative, with 0 < ghost_index < -10 for b-quark ghosts.
  
  int ghost_index = -1;
  
  for (int i_gbq = 0; i_gbq<Ghost_Bquarks.size(); i_gbq++)
    {
      fastjet::PseudoJet ghost_particle( (Ghost_Bquarks[i_gbq]).Px(), (Ghost_Bquarks[i_gbq]).Py(), (Ghost_Bquarks[i_gbq]).Pz(), (Ghost_Bquarks[i_gbq]).E() );

      ghost_particle.set_user_index(ghost_index);
      fjInputs.push_back( ghost_particle );
      ghost_index--;
    }


     // Check that there is some input to Fastjet
     // .........................................
  
  if (fjInputs.size() == 0)
    {
      cout << "Error: event with no final state particles" << endl;
      return;
    }

  
     // Run Fastjet algorithm
     // .....................
  
  vector <fastjet::PseudoJet> inclusiveJets, sortedJets;
  fastjet::ClusterSequence clustSeq(fjInputs, *jetDef);

  
     // Extract inclusive jets sorted by pT (note minimum pT of ptcut GeV)
     // .................................................................
  
  inclusiveJets = clustSeq.inclusive_jets(ptcut);
  sortedJets    = sorted_by_pt(inclusiveJets);

  
  // Fill the jet object and store it in the true jet collection
  // -----------------------------------------------------------

     // Loop over all jets produced by Fastjet
     // ......................................

  for (int i_jet = 0; i_jet < sortedJets.size(); i_jet++)
    {
      
     // Get the constituents of the jet
     // ...............................

      std::vector<fastjet::PseudoJet> constituents = (sortedJets[i_jet]).constituents();

      
      // Check if the jet is a B-jet
      // ...........................

      bool bq_ghosttag = false;
      int ghost_mult = 0;
      
      for (int i_const = 0; i_const<constituents.size(); i_const++)
	{
	  if ( (constituents[i_const]).user_index() < 0 && (constituents[i_const]).user_index() > -10 )
	    {
	      bq_ghosttag = true;
	      ghost_mult++;
	    }
	}
      
      int npart_injet = constituents.size() - ghost_mult;

      
      // Fill a TruthJets object
      // .......................
      
      TruthJets temp_truthjet;
      temp_truthjet.Set((sortedJets[i_jet]).e(),   (sortedJets[i_jet]).px(),  (sortedJets[i_jet]).py(), 
			(sortedJets[i_jet]).pz(),  (sortedJets[i_jet]).m(),   (sortedJets[i_jet]).pt(), 
			(sortedJets[i_jet]).eta(), (sortedJets[i_jet]).rap(), (sortedJets[i_jet]).phi(), 
			bq_ghosttag, false, npart_injet);

      
      // Put in the jet collection
      // .........................
      
      p_PartonJets_Coll->push_back(temp_truthjet);
      
  }

  
  // Close objects
  // -------------
  
  delete jetDef;
  
}
// ============================================================================================================




// ============================================================================================================
void ANA_utils::Get_BottomQuarks(Pythia8::Event event, std::vector<TruthPart>* p_BQuarks_Coll) {
//
// Find the bottom quarks that have come from a top, a Z, a W, or a Higgs. Fill a TruthPart object for each of
// the b-quarks that come from the above particles, just before it decays or hadronizes. Store it in the
// b-qaurk collection.
//  
//    Note: The b-quark must absolutely come from one of the aforementioned particles, but it should be the
//          the last in line before decaying or hadronizing, i.e. it must be after any gluon or photon
//          radiation, as well as not a carbon copy of the right b-quark. The strategy is therefore to find
//          the b-quarks immediately following the top/Z/W/H, and following their generations of descendants,
//          until there is no more b-quark in the next generation. A test showed that about 0.1% of the events
//          had up to the 13th generation, on a 1000 events test. The code therefore follow up to the 15th
//          generation of gluon/photon FSR or carbon copies, to make sure that there essentially no chances to
//          pick the wrong b-quarks. For most events, there are between 3 and 7 generations.   
//
// ============================================================================================================

  // Find the indices of the b-quarks directly coming from tops, W, Z, or H
  // ----------------------------------------------------------------------
  
  std::vector<int> bquark_uplevelindices;

  for (int i = 1; i < event.size(); i++)
    {
      if ((event[i].idAbs() == 5) &&
	  (event[event[i].mother1()].idAbs() == 6 || event[event[i].mother1()].idAbs() == 24 || event[event[i].mother1()].idAbs() == 23 || event[event[i].mother1()].idAbs() == 25) )
	{
	  bquark_uplevelindices.push_back(i);
	}
    }

  // Find the indices of the b-quarks just before they decay or hadronize
  // --------------------------------------------------------------------

  std::vector<int> bquark_indices;
  
  // Loop over all b-quarks directly coming from top/W/Z/H
  // .....................................................
  // Note: An integer will be used to keep the particle index of the last b-quark before decay or hadronization

  for (int i=0; i< bquark_uplevelindices.size(); i++)
    {
      
      int lastbquark_index=bquark_uplevelindices[i];

      // Get all the daughters of the considered b-quark and loop over them
      // ..................................................................
      // Note: Define an identifier to keep track of the generation that follows the last b-quark before decay and hadronization. That
      //       generatio is the first one that follows a top/W/Z/H decaying to b-quarks that do not contain any b-quark in it.

      std::vector<int> gen2 = event[lastbquark_index].daughterList();
      int topgen =1;
      if (gen2.size()>0) topgen = 2;
      for (int i2 = 0; i2<gen2.size();i2++)
	{

	// Check for b-quarks in the generation
	// ....................................
	  
	    // Note: If there are no b-quarks in this generation, there is nothing else to do because we already know which is the last b-quark index.
	    //       If there is a b-quark, move to the next generation by looking at its daughters.

	    // Note: This process is repeated over 15 generations, to make sure that we really get the last b-quark of the chain. The constraint on the
	    //       for loop (i<genX.size() )) guarantees that it is not possible to get further than this last generation.
	  
	  if (event[gen2[i2]].idAbs() != 5) continue;
	  else
	    {
	      lastbquark_index=gen2[i2];
	      std::vector<int> gen3 = event[lastbquark_index].daughterList();
	      if (gen3.size()>0) topgen = 3;
	      for (int i3 = 0; i3<gen3.size();i3++)
		{
		  
		  if (event[gen3[i3]].idAbs() !=5) continue;
		  else
		    {
		      lastbquark_index=gen3[i3];
		      std::vector<int> gen4 = event[lastbquark_index].daughterList();
		      if (gen4.size()>0) topgen = 4;
		      for (int i4 = 0; i4<gen4.size();i4++)
			{
			  if (event[gen4[i4]].idAbs() !=5) continue;
			  else
			    {
			      lastbquark_index=gen4[i4];
			      std::vector<int> gen5 = event[lastbquark_index].daughterList();
			      if (gen5.size()>0) topgen = 5;
			      for (int i5 = 0; i5<gen5.size();i5++)
				{
				  if (event[gen5[i5]].idAbs() !=5) continue;
				  else
				    {
				      lastbquark_index=gen5[i5];
				      std::vector<int> gen6 = event[lastbquark_index].daughterList();
				      if (gen6.size()>0) topgen = 6;
				      for (int i6 = 0; i6<gen6.size();i6++)
					{
					  if (event[gen6[i6]].idAbs() !=5) continue;
					  else
					    {
					      lastbquark_index=gen6[i6];
					      std::vector<int> gen7 = event[lastbquark_index].daughterList();
					      if (gen7.size()>0) topgen = 7;
					      for (int i7 = 0; i7<gen7.size();i7++)
						{
						  if (event[gen7[i7]].idAbs() !=5) continue;
						  else
						    {
						      lastbquark_index=gen7[i7];
						      std::vector<int> gen8 = event[lastbquark_index].daughterList();
						      if (gen8.size()>0) topgen = 8;
						      for (int i8 = 0; i8<gen8.size();i8++)
							{
							  if (event[gen8[i8]].idAbs() !=5) continue;
							  else
							    {
							      lastbquark_index=gen8[i8];
							      std::vector<int> gen9 = event[lastbquark_index].daughterList();
							      if (gen9.size()>0) topgen = 9;
							      for (int i9 = 0; i9<gen9.size();i9++)
								{
								  if (event[gen9[i9]].idAbs() !=5) continue;
								  else
								    {
								      lastbquark_index=gen9[i9];
								      std::vector<int> gen10 = event[lastbquark_index].daughterList();
								      if (gen10.size()>0) topgen = 10;
								      for (int i10 = 0; i10<gen10.size();i10++)
									{
									  if (event[gen10[i10]].idAbs() !=5) continue;
									  else
									    {
									      lastbquark_index=gen10[i10];
									      std::vector<int> gen11 = event[lastbquark_index].daughterList();
									      if (gen11.size()>0) topgen = 11;
									      for (int i11 = 0; i11<gen11.size();i11++)
										{
										  if (event[gen11[i11]].idAbs() !=5) continue;
										  else
										    {
										      lastbquark_index=gen11[i11];
										      std::vector<int> gen12 = event[lastbquark_index].daughterList();
										      if (gen12.size()>0) topgen = 12;
										      for (int i12 = 0; i12<gen12.size();i12++)
											{
											  if (event[gen12[i12]].idAbs() !=5) continue;
											  else
											    {
											      lastbquark_index=gen12[i12];
											      std::vector<int> gen13 = event[lastbquark_index].daughterList();
											      if (gen13.size()>0) topgen = 13;
											      for (int i13 = 0; i13<gen13.size();i13++)
												{
												  if (event[gen13[i13]].idAbs() !=5) continue;
												  else
												    {
												      lastbquark_index=gen13[i13];
												      std::vector<int> gen14 = event[lastbquark_index].daughterList();
												      if (gen14.size()>0) topgen = 14;
												      for (int i14 = 0; i14<gen14.size();i14++)
													{
													  if (event[gen14[i14]].idAbs() !=5) continue;
													  else
													    {
													      lastbquark_index=gen14[i14];
													      std::vector<int> gen15 = event[lastbquark_index].daughterList();
													      if (gen15.size()>0) topgen = 15;
													    }
													} // End Generation 14
												    }
												} // End Generation 13
											    }
											} // End Generation 12
										    }
										} // End Generation 11
									    }
									} // End Generation 10
								    }
								} // End Generation 9
							    }
							} // End Generation 8
						    }
						} // End Generation 7
					    }
					} // End Generation 6
				    }
				} // End Generation 5
			    }
			} // End Generation 4
		    }
		} // End Generation 3
	    }
	} // End Generation 2
      

      //  Push the last index in the vector of b-quarks indices for which the b-quark will now decay or hadronize
      //  .......................................................................................................
      
      bquark_indices.push_back(lastbquark_index);

    }

  // Fill a TruthPart object for each b-quark, and push it back into the b-quark collection
  // -------------------------------------------------------------------------------------
  
  for (auto i : bquark_indices)
    {
      
      TruthPart temp_b;
      TruthPart* p_temp_b = &temp_b;

      Fill_TruthPart(event, i, p_temp_b);
      p_BQuarks_Coll->push_back(temp_b);
    }

}
// ============================================================================================================




// ============================================================================================================
void ANA_utils::Get_BottomHadrons(Pythia8::Event event, std::vector<TruthPart> BQuark_Coll, std::vector<TruthPart>* p_BHadrons_Coll) {
//
// Starting from a collection of b-quarks, this function returns a collection of B-hadrons that come from these
// b-quarks.  
//
// ============================================================================================================

  // Create a list of all B-hadrons PdgIds
  // -------------------------------------
  
  int BhadronsPdgId[88] = {511, 521, 10511, 10521, 513, 523, 10513, 10523, 20513, 20523, 515, 525, 531, 10531, 533, 10533, 20553, 535, 541, 10541, 543, 10543, 20543, 545, 551, 10551, 100551, 110551, 200551, 210551, 553, 10553, 20553, 30553, 100553, 110553, 120553, 130553, 200553, 210553, 220553, 300553, 9000553, 9010553, 555, 10555, 20555, 100555, 110555, 120555, 200555, 557, 100557, 5122, 5112, 5212, 5222, 5114, 5214, 5224, 5132, 5232, 5312, 5322, 5314, 5324, 5332, 5334, 5142, 5242, 5412, 5422, 5414, 5424, 5342, 5432, 5434, 5442, 5444, 5512, 5522, 5514, 5524, 5532, 5534, 5542, 5544, 5554};


  // Loop over all the bquarks coming from top/W/Z/H
  // -----------------------------------------------

  for (int i_bq=0; i_bq<BQuark_Coll.size(); i_bq++)
    {

      // Find the list of daughters of this b-quark
      // ------------------------------------------

      int bquark_index = (BQuark_Coll[i_bq]).Index();
      std::vector<int> bquark_daughterlist = event[bquark_index].daughterList();

      // Find the particle index of the B-hadrons following the b-quark hadronization
      // ----------------------------------------------------------------------------

      // Loop over the daughters of the b-quark
      // ......................................

      for (int i_bd = 0; i_bd<bquark_daughterlist.size(); i_bd++)
	{
	  int temp_index = bquark_daughterlist[i_bd];

          // Check if the particle is in the list of possible B-hadrons and if yes, push this particle back in the B-Hadron Collection
          // .........................................................................................................................

	  int temp_id = event[temp_index].id();

	  for (auto i_list : BhadronsPdgId)
	    {

	      if ( i_list == fabs(temp_id) )
		{
		  TruthPart temp_bhad;
		  TruthPart* p_temp_bhad = &temp_bhad;

		  Fill_TruthPart(event, temp_index, p_temp_bhad);
		  p_BHadrons_Coll->push_back(temp_bhad);
		}
	    }
	}
    }
}
// ============================================================================================================
