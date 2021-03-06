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
// 070523 bug fix for G4FPE_DEBUG on by A. Howard ( and T. Koi)
// 080808 bug fix in Sample() and GetXsec() by T. Koi
//
#include "G4NeutronHPVector.hh"
#include "G4SystemOfUnits.hh"

  // if the ranges do not match, constant extrapolation is used.
  G4NeutronHPVector & operator + (G4NeutronHPVector & left, G4NeutronHPVector & right)
  {
    G4NeutronHPVector * result = new G4NeutronHPVector;
    G4int j=0;
    G4double x;
    G4double y;
    G4int running = 0;
    for(G4int i=0; i<left.GetVectorLength(); i++)
    {
      while(j<right.GetVectorLength())
      {
        if(right.GetX(j)<left.GetX(i)*1.001)
        {
          x = right.GetX(j);
          y = right.GetY(j)+left.GetY(x);
          result->SetData(running++, x, y);
          j++;
        }
        //else if(std::abs((right.GetX(j)-left.GetX(i))/(left.GetX(i)+right.GetX(j)))>0.001)
        else if( left.GetX(i)+right.GetX(j) == 0
              || std::abs((right.GetX(j)-left.GetX(i))/(left.GetX(i)+right.GetX(j))) > 0.001 )
        {
          x = left.GetX(i);
          y = left.GetY(i)+right.GetY(x);
          result->SetData(running++, x, y);
          break;
        }
        else
        {
          break;
        }
      }
      if(j==right.GetVectorLength())
      {
        x = left.GetX(i);
        y = left.GetY(i)+right.GetY(x);
        result->SetData(running++, x, y);
      }
    }
    result->ThinOut(0.02);
    return *result;
  }

  G4NeutronHPVector::G4NeutronHPVector()
  {
    theData = new G4NeutronHPDataPoint[20];
    nPoints=20;
    nEntries=0;
    Verbose=0;
    theIntegral=0;
    totalIntegral=-1;
    isFreed = 0;
    maxValue = -DBL_MAX;
    the15percentBorderCash = -DBL_MAX;
    the50percentBorderCash = -DBL_MAX;
    label = -DBL_MAX;

  }

  G4NeutronHPVector::G4NeutronHPVector(G4int n)
  {
    nPoints=std::max(n, 20);
    theData = new G4NeutronHPDataPoint[nPoints];
    nEntries=0;
    Verbose=0;
    theIntegral=0;
    totalIntegral=-1;
    isFreed = 0;
    maxValue = -DBL_MAX;
    the15percentBorderCash = -DBL_MAX;
    the50percentBorderCash = -DBL_MAX;
  }

  G4NeutronHPVector::~G4NeutronHPVector()
  {
//    if(Verbose==1)G4cout <<"G4NeutronHPVector::~G4NeutronHPVector"<<G4endl;
      delete [] theData;
//    if(Verbose==1)G4cout <<"Vector: delete theData"<<G4endl;
      delete [] theIntegral;
//    if(Verbose==1)G4cout <<"Vector: delete theIntegral"<<G4endl;
    theHash.Clear();
    isFreed = 1;
  }

  G4NeutronHPVector & G4NeutronHPVector::
  operator = (const G4NeutronHPVector & right)
  {
    if(&right == this) return *this;

    G4int i;

    totalIntegral = right.totalIntegral;
    if(right.theIntegral!=0) theIntegral = new G4double[right.nEntries];
    for(i=0; i<right.nEntries; i++)
    {
      SetPoint(i, right.GetPoint(i)); // copy theData
      if(right.theIntegral!=0) theIntegral[i] = right.theIntegral[i];
    }
    theManager = right.theManager;
    label = right.label;

    Verbose = right.Verbose;
    the15percentBorderCash = right.the15percentBorderCash;
    the50percentBorderCash = right.the50percentBorderCash;
    theHash = right.theHash;
   return *this;
  }


  G4double G4NeutronHPVector::GetXsec(G4double e)
  {
    if(nEntries == 0) return 0;
    if(!theHash.Prepared()) Hash();
    G4int min = theHash.GetMinIndex(e);
    G4int i;
    for(i=min ; i<nEntries; i++)
    {
      //if(theData[i].GetX()>e) break;
      if(theData[i].GetX() >= e) break;
    }
    G4int low = i-1;
    G4int high = i;
    if(i==0)
    {
      low = 0;
      high = 1;
    }
    else if(i==nEntries)
    {
      low = nEntries-2;
      high = nEntries-1;
    }
    G4double y;
    if(e<theData[nEntries-1].GetX())
    {
      // Protect against doubled-up x values
      //if( (theData[high].GetX()-theData[low].GetX())/theData[high].GetX() < 0.000001)
      if ( theData[high].GetX() !=0
       //080808 TKDB
       //&&( theData[high].GetX()-theData[low].GetX())/theData[high].GetX() < 0.000001)
       &&( std::abs( (theData[high].GetX()-theData[low].GetX())/theData[high].GetX() ) < 0.000001 ) )
      {
        y = theData[low].GetY();
      }
      else
      {
        y = theInt.Interpolate(theManager.GetScheme(high), e,
                               theData[low].GetX(), theData[high].GetX(),
		  	       theData[low].GetY(), theData[high].GetY());
      }
    }
    else
    {
      y=theData[nEntries-1].GetY();
    }
    return y;
  }

  void G4NeutronHPVector::Dump()
  {
    G4cout << nEntries<<G4endl;
    for(G4int i=0; i<nEntries; i++)
    {
      G4cout << theData[i].GetX()<<" ";
      G4cout << theData[i].GetY()<<" ";
//      if (i!=1&&i==5*(i/5)) G4cout << G4endl;
      G4cout << G4endl;
    }
    G4cout << G4endl;
  }

  void G4NeutronHPVector::Check(G4int i)
  {
    if(i>nEntries) throw G4HadronicException(__FILE__, __LINE__, "Skipped some index numbers in G4NeutronHPVector");
    if(i==nPoints)
    {
      nPoints = static_cast<G4int>(1.2*nPoints);
      G4NeutronHPDataPoint * buff = new G4NeutronHPDataPoint[nPoints];
      for (G4int j=0; j<nEntries; j++) buff[j] = theData[j];
      delete [] theData;
      theData = buff;
    }
    if(i==nEntries) nEntries=i+1;
  }

  void G4NeutronHPVector::
  Merge(G4InterpolationScheme aScheme, G4double aValue,
        G4NeutronHPVector * active, G4NeutronHPVector * passive)
  {
    // interpolate between labels according to aScheme, cut at aValue,
    // continue in unknown areas by substraction of the last difference.

    CleanUp();
    G4int s_tmp = 0, n=0, m_tmp=0;
    G4NeutronHPVector * tmp;
    G4int a = s_tmp, p = n, t;
    while ( (a<active->GetVectorLength()) && (p<passive->GetVectorLength()) )
    {
      if(active->GetEnergy(a) < passive->GetEnergy(p))
      {
        G4double xa  = active->GetEnergy(a);
        G4double yy = theInt.Interpolate(aScheme, aValue, active->GetLabel(), passive->GetLabel(),
                                                          active->GetXsec(a), passive->GetXsec(xa));
        SetData(m_tmp, xa, yy);
        theManager.AppendScheme(m_tmp, active->GetScheme(a));
        m_tmp++;
        a++;
      }
      else if(active->GetEnergy(a) == passive->GetEnergy(p))
      {
        G4double xa  = passive->GetEnergy(p);
        G4double yy = theInt.Interpolate(aScheme, aValue, active->GetLabel(), passive->GetLabel(),
                                                          active->GetXsec(a), passive->GetXsec(p));
        SetData(m_tmp, xa, yy);
        theManager.AppendScheme(m_tmp, passive->GetScheme(p));
        m_tmp++;
        p++;
        a++;
      }
      else
      {
        G4double xa  = passive->GetEnergy(p);
        G4double yy = theInt.Interpolate(aScheme, aValue, active->GetLabel(), passive->GetLabel(),
                                                          active->GetXsec(xa), passive->GetXsec(p));
        SetData(m_tmp, xa, yy);
        theManager.AppendScheme(m_tmp, passive->GetScheme(p));
        m_tmp++;
        p++;
      }
    }
    while (a<active->GetVectorLength())
    {
        G4double xa  = active->GetEnergy(a);
        G4double yy = theInt.Interpolate(aScheme, aValue, active->GetLabel(), passive->GetLabel(),
                                                          active->GetXsec(a), passive->GetXsec(xa));
        SetData(m_tmp, xa, yy);
        theManager.AppendScheme(m_tmp, active->GetScheme(a));
        m_tmp++;
        a++;
    }
    while ( p<passive->GetVectorLength() )
    {
        G4double xa  = passive->GetEnergy(p);
        G4double yy = theInt.Interpolate(aScheme, aValue, active->GetLabel(), passive->GetLabel(),
                                                          active->GetXsec(xa), passive->GetXsec(p));
        SetData(m_tmp, xa, yy);
        theManager.AppendScheme(m_tmp, passive->GetScheme(p));
        m_tmp++;
        p++;
    }

    // Rebuild the Hash;
    if(theHash.Prepared())
    {
      ReHash();
    }
  }

  void G4NeutronHPVector::ThinOut(G4double precision)
  {
    // anything in there?
    if(GetVectorLength()==0) return;
    // make the new vector
    G4NeutronHPDataPoint * aBuff = new G4NeutronHPDataPoint[nPoints];
    G4double x, x1, x2, y, y1, y2;
    G4int count = 0, current = 2, start = 1;

    // First element always goes and is never tested.
    aBuff[0] = theData[0];

    // Find the rest
    while(current < GetVectorLength())
    {
      x1=aBuff[count].GetX();
      y1=aBuff[count].GetY();
      x2=theData[current].GetX();
      y2=theData[current].GetY();
      for(G4int j=start; j<current; j++)
      {
	x = theData[j].GetX();
	if(x1-x2 == 0) y = (y2+y1)/2.;
	else y = theInt.Lin(x, x1, x2, y1, y2);
	if (std::abs(y-theData[j].GetY())>precision*y)
	{
	  aBuff[++count] = theData[current-1]; // for this one, everything was fine
          start = current; // the next candidate
	  break;
	}
      }
      current++ ;
    }
    // The last one also always goes, and is never tested.
    aBuff[++count] = theData[GetVectorLength()-1];
    delete [] theData;
    theData = aBuff;
    nEntries = count+1;

    // Rebuild the Hash;
    if(theHash.Prepared())
    {
      ReHash();
    }
  }

  G4bool G4NeutronHPVector::IsBlocked(G4double aX)
  {
    G4bool result = false;
    std::vector<G4double>::iterator i;
    for(i=theBlocked.begin(); i!=theBlocked.end(); i++)
    {
      G4double aBlock = *i;
      if(std::abs(aX-aBlock) < 0.1*MeV)
      {
        result = true;
	theBlocked.erase(i);
	break;
      }
    }
    return result;
  }

  G4double G4NeutronHPVector::Sample() // Samples X according to distribution Y
  {
    G4double result;
    G4int j;
    for(j=0; j<GetVectorLength(); j++)
    {
      if(GetY(j)<0) SetY(j, 0);
    }

    if(theBuffered.size() !=0 && G4UniformRand()<0.5)
    {
      result = theBuffered[0];
      theBuffered.erase(theBuffered.begin());
      if(result < GetX(GetVectorLength()-1) ) return result;
    }
    if(GetVectorLength()==1)
    {
      result = theData[0].GetX();
    }
    else
    {
      if(theIntegral==0) { IntegrateAndNormalise(); }
      do
      {
//080808
/*
        G4double rand;
        G4double value, test, baseline;
        baseline = theData[GetVectorLength()-1].GetX()-theData[0].GetX();
        do
        {
          value = baseline*G4UniformRand();
          value += theData[0].GetX();
          test = GetY(value)/maxValue;
          rand = G4UniformRand();
        }
        //while(test<rand);
        while( test < rand && test > 0 );
        result = value;
*/
        G4double rand;
        G4double value, test;
        do
        {
           rand = G4UniformRand();
           G4int ibin = -1;
           for ( G4int i = 0 ; i < GetVectorLength() ; i++ )
           {
              if ( rand < theIntegral[i] )
              {
                 ibin = i;
                 break;
              }
           }
           if ( ibin < 0 ) G4cout << "TKDB 080807 " << rand << G4endl;
           // result
           rand = G4UniformRand();
           G4double x1, x2;
           if ( ibin == 0 )
           {
              x1 = theData[ ibin ].GetX();
              value = x1;
              break;
           }
           else
           {
              x1 = theData[ ibin-1 ].GetX();
           }

           x2 = theData[ ibin ].GetX();
           value = rand * ( x2 - x1 ) + x1;
	   //***********************************************************************
	   /*
           test = GetY ( value ) / std::max ( GetY( ibin-1 ) , GetY ( ibin ) );
	   */
	   //***********************************************************************
	   //EMendoza - Always linear interpolation:
	   G4double y1=theData[ ibin-1 ].GetY();
           G4double y2=theData[ ibin ].GetY();
	   G4double mval=(y2-y1)/(x2-x1);
	   G4double bval=y1-mval*x1;
	   test =(mval*value+bval)/std::max ( GetY( ibin-1 ) , GetY ( ibin ) );
	   //***********************************************************************
        }
        while ( G4UniformRand() > test );
        result = value;
//080807
      }
      while(IsBlocked(result));
    }
    return result;
  }

  G4double G4NeutronHPVector::Get15percentBorder()
  {
    if(the15percentBorderCash>-DBL_MAX/2.) return the15percentBorderCash;
    G4double result;
    if(GetVectorLength()==1)
    {
      result = theData[0].GetX();
      the15percentBorderCash = result;
    }
    else
    {
      if(theIntegral==0) { IntegrateAndNormalise(); }
      G4int i;
      result = theData[GetVectorLength()-1].GetX();
      for(i=0;i<GetVectorLength();i++)
      {
	if(theIntegral[i]/theIntegral[GetVectorLength()-1]>0.15)
	{
	  result = theData[std::min(i+1, GetVectorLength()-1)].GetX();
          the15percentBorderCash = result;
	  break;
	}
      }
      the15percentBorderCash = result;
    }
    return result;
  }

  G4double G4NeutronHPVector::Get50percentBorder()
  {
    if(the50percentBorderCash>-DBL_MAX/2.) return the50percentBorderCash;
    G4double result;
    if(GetVectorLength()==1)
    {
      result = theData[0].GetX();
      the50percentBorderCash = result;
    }
    else
    {
      if(theIntegral==0) { IntegrateAndNormalise(); }
      G4int i;
      G4double x = 0.5;
      result = theData[GetVectorLength()-1].GetX();
      for(i=0;i<GetVectorLength();i++)
      {
	if(theIntegral[i]/theIntegral[GetVectorLength()-1]>x)
	{
	  G4int it;
	  it = i;
	  if(it == GetVectorLength()-1)
	  {
	    result = theData[GetVectorLength()-1].GetX();
	  }
	  else
	  {
	    G4double x1, x2, y1, y2;
	    x1 = theIntegral[i-1]/theIntegral[GetVectorLength()-1];
	    x2 = theIntegral[i]/theIntegral[GetVectorLength()-1];
	    y1 = theData[i-1].GetX();
	    y2 = theData[i].GetX();
	    result = theLin.Lin(x, x1, x2, y1, y2);
	  }
          the50percentBorderCash = result;
	  break;
	}
      }
      the50percentBorderCash = result;
    }
    return result;
  }
