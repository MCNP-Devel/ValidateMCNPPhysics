//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// neutron_hp -- source file
// J.P. Wellisch, Nov-1996
// A prototype of the low energy neutron transport model.
//
#include "../include/VFDNeutronHPFissionBaseFS.hh"
#include "G4SystemOfUnits.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4Proton.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4NeutronHPDataUsed.hh"

  void VFDNeutronHPFissionBaseFS::Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & bit)
  {
    G4String tString = dirName;
    G4bool dbool;
    G4NeutronHPDataUsed aFile = theNames.GetName(static_cast<G4int>(A), static_cast<G4int>(Z), M, tString, bit, dbool);
    G4String filename = aFile.GetName();
    SetAZMs( A, Z, M, aFile );
    //theBaseA = aFile.GetA();
    //theBaseZ = aFile.GetZ();
    //if(!dbool  || ( Z<2.5 && ( std::abs(theBaseZ - Z)>0.0001 || std::abs(theBaseA - A)>0.0001) ) )
    if ( !dbool || ( Z<2.5 && ( std::abs(theNDLDataZ - Z)>0.0001 || std::abs(theNDLDataA - A)>0.0001)) )
    {
      hasAnyData = false;
      hasFSData = false;
      hasXsec = false;
      return; // no data for exactly this isotope.
    }

    std::ifstream theData(filename, std::ios::in);
    G4int dummy;
    if(!(theData))
    {
      theData.close();
      hasFSData = false;
      hasXsec = false;
      hasAnyData = false;
      return; // no data for this FS for this isotope
    }
    theData >> dummy>>dummy;
    G4int total;
    theData >> total;
    theXsection->Init(theData, total, eV);
    theXsection->ThinOut(0.02);
    if (!(theData >> dummy))
    {
      hasFSData = false;
      theData.close();
      return;
    }
    theData >> dummy;

    theAngularDistribution.Init(theData);

    theData >> dummy >> dummy;

    theEnergyDistribution.Init(theData);
    theData.close();

  }

G4HadFinalState * VFDNeutronHPFissionBaseFS::ApplyYourself(const G4HadProjectile &aTrack)
  {
    theResult.Clear();
    theNeutron.SetDefinition(const_cast<G4ParticleDefinition *>(aTrack.GetDefinition()));
    theNeutron.SetMomentum( aTrack.Get4Momentum().vect() );
    theNeutron.SetKineticEnergy( aTrack.GetKineticEnergy() );

    G4Nucleus aNucleus;
   theTarget = aNucleus.GetBiasedThermalNucleus( targetMass, G4ThreeVector(1.,0.,0.), 0.);

    G4int nPrompt=3;
// if therere were no data for this isotope, break out.
    if(!HasFSData()) { return 0; }

    G4int i;
    G4ReactionProduct boosted;
    boosted.Lorentz(theNeutron, theTarget);
    G4double eKinetic = boosted.GetKineticEnergy();

    theAngularDistribution.SetNeutron(theNeutron);
    theAngularDistribution.SetTarget(theTarget);

// Build neutrons
    G4ReactionProduct * theNeutrons = new G4ReactionProduct[nPrompt];
    for(i=0; i<nPrompt; i++)
    {
      theNeutrons[i].SetDefinition(G4Neutron::Neutron());
    }

// sample energies
   G4int dummy;
   for(i=0; i<nPrompt; i++)
   {
     // always in the lab system (if file-5)
     theNeutrons[i].SetKineticEnergy(theEnergyDistribution.Sample(eKinetic, dummy));
   }

// sample neutron angular distribution
   for(i=0; i<nPrompt; i++)
   {
     theAngularDistribution.SampleAndUpdate(theNeutrons[i]);
   }

// already in lab. Add neutrons to dynamic particle vector
   for(i=0; i<nPrompt; i++)
   {
      G4DynamicParticle * it = new G4DynamicParticle;
      it->SetDefinition(theNeutrons[i].GetDefinition());
      it->SetMomentum(theNeutrons[i].GetMomentum());
      theResult.AddSecondary(it);
   }
   delete [] theNeutrons;

//    for(int x=0; x<theResult.GetNumberOfSecondaries(); x++)
//    {
//        delete theResult.GetSecondary(x)->GetParticle();
//    }
//    theResult.Clear();
   return &theResult;
  }
