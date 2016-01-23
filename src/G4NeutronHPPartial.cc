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
// 13-Jan-06 fix in Sample by T. Koi
//
#include "G4NeutronHPPartial.hh"
#include "G4NeutronHPInterpolator.hh"
#include "Randomize.hh"

G4NeutronHPVector * G4NeutronHPPartial::GetY(G4double e1)
  {
    G4NeutronHPVector * aBuffer = new G4NeutronHPVector();
    G4int i;
    if(nData==1)
    {
      for(i=0; i<data[0].GetVectorLength(); i++)
      {
        aBuffer->SetInterpolationManager(data[0].GetInterpolationManager());
        aBuffer->SetData(i , data[0].GetX(i), data[0].GetY(i));
      }
      return aBuffer;
    }
    for (i=0; i<nData; i++)
    {
       if(X[i]>e1) break;
    }
    if(i==nData) i--;
    if(0==i) i=1;
    data[i-1].GetLabel();
    data[i].GetLabel();
    aBuffer->Merge(theManager.GetScheme(i),e1,data+i-1,data+i);
    return aBuffer;
  }

  G4double G4NeutronHPPartial::Sample(G4double x)
  {
    G4NeutronHPVector * aBuffer = GetY(x);
    G4double result = aBuffer->Sample();
    delete aBuffer;
    return result;
  }
