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
//
//
#ifndef VFDNeutronHPFissionBaseFS_h
#define VFDNeutronHPFissionBaseFS_h 1

#include "globals.hh"
#include "G4ReactionProduct.hh"
#include "G4DynamicParticleVector.hh"
#include "VFDNeutronHPFinalState.hh"
#include "G4NeutronHPNames.hh"
#include "G4NeutronHPVector.hh"
#include "VFDNeutronHPEnergyDistribution.hh"
#include "G4NeutronHPAngular.hh"
#include "G4Track.hh"
#include "G4IonTable.hh"
#include "VFDNeutronHPFSFissionFS.hh"

class VFDNeutronHPFissionBaseFS : public VFDNeutronHPFinalState
{
    public:

    VFDNeutronHPFissionBaseFS()
    {
        hasXsec = true;
        theXsection = new G4NeutronHPVector;
    }
    virtual ~VFDNeutronHPFissionBaseFS()
    {
        delete theXsection;
    }

    void Init (G4double A, G4double Z, G4int M, G4String & dirName, G4String & bit);

    G4HadFinalState * ApplyYourself(const G4HadProjectile &aTrack);

    virtual G4double GetXsec(G4double anEnergy)
    {
        return std::max(0., theXsection->GetY(anEnergy));
    }
    virtual G4NeutronHPVector * GetXsec() { return theXsection; }

    inline void SetNeutron(const G4ReactionProduct & aNeutron)
    {
      theNeutron = aNeutron;
      theAngularDistribution.SetNeutron(aNeutron);
    }

    inline void SetTarget(const G4ReactionProduct & aTarget)
    {
      theTarget = aTarget;
      theAngularDistribution.SetTarget(aTarget);
    }
    void SetMass(double mass)
    {
        targetMass = mass;
    }

    VFDNeutronHPFinalState * New()
    {
        VFDNeutronHPFissionBaseFS * theNew = new VFDNeutronHPFissionBaseFS;
        return theNew;
    }

    G4double GetMinEnergy()
    {
        if(theXsection->GetVectorLength()==0)
            return 0.;
        else
            return theXsection->GetEnergy(0);
    }

    G4double GetMaxEnergy()
    {
        if(theXsection->GetVectorLength()==0)
            return 20.;
        else
        {
            int count=1;
            while(theXsection->GetXsec(theXsection->GetVectorLength()-count)==0.0)
            {
                count++;
            }
            return theXsection->GetEnergy(theXsection->GetVectorLength()-count);
        }
    }

  private:


  G4NeutronHPVector * theXsection;
  VFDNeutronHPEnergyDistribution theEnergyDistribution;
  G4NeutronHPAngular theAngularDistribution;
  VFDNeutronHPFSFissionFS theYieldData;

  G4ReactionProduct theNeutron;
  G4ReactionProduct theTarget;
  G4double targetMass;
  G4int offset;

  private:

};
#endif
