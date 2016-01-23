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
// $Id: VFDIonTable.hh 79333 2014-02-24 10:36:17Z gcosmo $
//
//
// ------------------------------------------------------------
//	GEANT 4 class header file
//
//	History: first implementation,
//      based on object model of June 27, 98 H.Kurashige
// ------------------------------------------------------------
//      added clear()			20 Mar., 08 H.Kurashige
//      modified GetIon                 02 Aug., 98 H.Kurashige
//      added Remove()                  06 Nov.,98 H.Kurashige
//      add GetNucleusMass              15 Mar. 99  H.Kurashige
//          -----
//      Modified GetIon methods                  17 Aug. 99 H.Kurashige
//      New design using G4VIsotopeTable          5 Oct. 99 H.Kurashige
//      Add GetNucleusEncoding according PDG 2006 9 Oct. 2006 H.Kurashige
//      Use STL map                              30 Jul. 2009 H.Kurashige
//      Add G4IsomerTable                        5 May. 2013  H.Kurashige
//      Add GetIsomerMass                       25 July 2013  H.Kurashige
//
#ifndef VFDIonTable_h
#define VFDIonTable_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"

#include <cmath>
#include <vector>
#include <map>

class G4ParticleTable;
class G4VIsotopeTable;
class G4IsotopeProperty;
class G4IsomerTable;
class G4NuclideTable;

class VFDIonTable : public G4IonTable
{
 // Class Description
 //   VFDIonTable is the table of pointer to G4ParticleDefinition
 //   In VFDIonTable, each G4ParticleDefinition pointer is stored
 //

 public:
   // Use STL map as list of ions
   typedef  std::multimap<G4int, const G4ParticleDefinition*> G4IonList;
   typedef  std::multimap<G4int, const G4ParticleDefinition*>::iterator G4IonListIterator;

 public:
  // constructor
   VFDIonTable(const G4IonTable *copyIonTab);

 protected:

 public:
  // destructor
   virtual ~VFDIonTable();

   void PreloadNuclideNew();
   // All nuclide with a life time longer than certain value will be created
   // prior to the event loop.


 private:
    G4ParticleDefinition **pDefList;
   G4NuclideTable* pNuclideTable;
   G4bool         isIsomerCreated;
   // Isomer table and flag of creation
};
#endif
