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
// $Id: VFDIonTable.cc 79333 2014-02-24 10:36:17Z gcosmo $
//
//
// --------------------------------------------------------------
//	GEANT 4 class implementation file
//
//	History: first implementation, based on object model of
//	27 June 1998  H.Kurashige
// ---------------------------------------------------------------
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige
//      use G4NucleiPropoerties to get nuceli Mass 17  Nov.,98 H.Kurashige
//      use G4GenericIon for process List
//      modify fomula of Ion mass       09 Dec., 98 H.Kurashige
//          -----
//      Modified GetIon methods         17 Aug. 99 H.Kurashige
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige
//      Modified Element Name for Z>103  06 Apr. 01 H.Kurashige
//      Remove test of cuts in SetCuts   16 Jan  03 V.Ivanchenko
//      Add G4IsomerTable                        5 May. 2013  H.Kurashige

#include <iostream>
#include <iomanip>
#include <sstream>

#include "G4ios.hh"
#include "G4Threading.hh"

#include "VFDIonTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4StateManager.hh"
#include "G4Ions.hh"
#include "G4UImanager.hh"
#include "G4NucleiProperties.hh"
#include "G4HyperNucleiProperties.hh"

#include "G4IsotopeProperty.hh"
#include "G4VIsotopeTable.hh"
#include "G4IsomerTable.hh"
#include "G4NuclideTable.hh"

// It is very important for multithreaded Geant4 to keep only one copy of the
// particle table pointer and the ion table pointer. However, we try to let
// each worker thread hold its own copy of the particle dictionary and the
// ion list. This implementation is equivalent to make the ion table thread
// private. The two shadow ponters are used by each worker thread to copy the
// content from the master thread.
//


////////////////////
// can't use the copy constructor of G4IonTable since it has not been defined in the GEANT4 class
VFDIonTable::VFDIonTable(const G4IonTable *copyIonTab)
 : G4IonTable::G4IonTable(), pNuclideTable(0)
{
    pDefList=0;
    isIsomerCreated=0;
}


////////////////////
VFDIonTable::~VFDIonTable()
{
  if((pNuclideTable!=0)&&(pDefList!=0))
  {
    for (int i = 0; i< int(pNuclideTable->entries()); ++i)
    {
      if( pDefList[i] != 0 )
        delete pDefList[i];
    }
    delete [] pDefList;
  }
}

////////////////////
void VFDIonTable::PreloadNuclideNew()
{
  if (isIsomerCreated) return;

  if (pNuclideTable==0)
  {
    pNuclideTable = G4NuclideTable::GetNuclideTable();
    pNuclideTable->GenerateNuclide();
    RegisterIsotopeTable(pNuclideTable);
  }

    pDefList = new G4ParticleDefinition *[pNuclideTable->entries()];
  for ( size_t i = 0 ; i != pNuclideTable->entries() ; i++ ) {
     const G4IsotopeProperty*  fProperty = pNuclideTable->GetIsotopeByIndex( i );
     G4int Z  = fProperty->GetAtomicNumber();
     G4int A  = fProperty->GetAtomicMass();
     G4double Eex  = fProperty->GetEnergy();
     pDefList[i]=GetIon(Z,A,Eex);
  }

  isIsomerCreated = true;
}


