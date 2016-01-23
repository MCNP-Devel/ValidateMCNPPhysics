#include "../include/FSSpectrumData.hh"


FSSpectrumData::FSSpectrumData()
{
    //ctor
}

FSSpectrumData::~FSSpectrumData()
{
    //dtor
}

void FSSpectrumData::AddData(G4HadFinalState* result)
{
    G4DynamicParticle *particle;
    int numSec, pSecY, nSecY, nDelY;
    pSecY = nSecY = nDelY = 0;

    numSec=result->GetNumberOfSecondaries();
    AddPrimary(result[0]);

    for(int s=0; s<numSec; s++)
    {
        particle=result->GetSecondary(s)->GetParticle();
        if(particle->GetParticleDefinition()->GetParticleName()=="gamma")
        {
            pSecY++;
            AddPhoton(particle[0]);
        }
        else if(particle->GetParticleDefinition()->GetParticleName()=="neutron")
        {
            if(result->GetSecondary(s)->GetTime()>0.)
            {
                nDelY++;
                AddDelayed(particle[0]);
            }
            else
            {
                nSecY++;
                AddSecondary(particle[0]);
            }
        }
        delete particle;
    }
    SetSecYield(nSecY);
    SetDelYield(nDelY);
    SetPhYield(pSecY);

    result->Clear();
}

void FSSpectrumData::AddSecondary(G4DynamicParticle &nSec)
{
    nSecMomAngle.push_back(std::sqrt(std::pow(nSec.GetMomentumDirection().getPhi(),2)+std::pow(nSec.GetMomentumDirection().getTheta(),2)));
    nSecKEn.push_back(std::log10(nSec.GetKineticEnergy()));
}

void FSSpectrumData::AddDelayed(G4DynamicParticle &nDel)
{
    nDelMomAngle.push_back(std::sqrt(std::pow(nDel.GetMomentumDirection().getPhi(),2)+std::pow(nDel.GetMomentumDirection().getTheta(),2)));
    nDelKEn.push_back(std::log10(nDel.GetKineticEnergy()));
}

void FSSpectrumData::AddPhoton(G4DynamicParticle &pSec)
{
    pSecMomAngle.push_back(std::sqrt(std::pow(pSec.GetMomentumDirection().getPhi(),2)+std::pow(pSec.GetMomentumDirection().getTheta(),2)));
    pSecKEn.push_back(std::log10(pSec.GetKineticEnergy()));
}

double FSSpectrumData::CompareFSData(std::string &outFileName, FSSpectrumData &mcnpData, bool *relevantData)
{
    std::stringstream stream;
    double totalDiff=0.;

    if(relevantData[2])
    {
        stream << "Comparing the secondary neutron deflection angle" << std::endl;
        totalDiff+=CompareHist(stream, nSecMomAngle, mcnpData.GetNSecMomAngle());
    }

    if(relevantData[4])
    {
        stream << "Comparing the delayed neutron deflection angle" << std::endl;
        totalDiff+=CompareHist(stream, nDelMomAngle, mcnpData.GetNDelMomAngle());
    }

    if(relevantData[6])
    {
        stream << "Comparing the secondary photon deflection angle" << std::endl;
        totalDiff+=CompareHist(stream, pSecMomAngle, mcnpData.GetPSecMomAngle());
    }

    if(relevantData[9])
    {
        stream << "Comparing the secondary neutron kinetic energy" << std::endl;
        totalDiff+=CompareHist(stream, nSecKEn, mcnpData.GetNSecKEn());
    }

    if(relevantData[10])
    {
        stream << "Comparing the delayed neutron kinetic energy" << std::endl;
        totalDiff+=CompareHist(stream, nDelKEn, mcnpData.GetNDelKEn());
    }

    if(relevantData[11])
    {
        stream << "Comparing the secondary photon kinetic energy" << std::endl;
        totalDiff+=CompareHist(stream, pSecKEn, mcnpData.GetPSecKEn());
    }

    if(relevantData[12])
    {
        stream << "Comparing the secondary neutron yield" << std::endl;
        totalDiff+=CompareHist(stream, nSecYield, mcnpData.GetNSecYield());
    }

    if(relevantData[13])
    {
        stream << "Comparing the delayed neutron yield" << std::endl;
        totalDiff+=CompareHist(stream, nDelYield, mcnpData.GetNDelYield());
    }

    if(relevantData[14])
    {
        stream << "Comparing the secondary photon yield" << std::endl;
        totalDiff+=CompareHist(stream, pSecYield, mcnpData.GetPSecYield());
    }

    SetDataStream( outFileName, stream, false );
    return totalDiff;
}

double FSSpectrumData::CompareFSData(std::string &outFileName, FSSpectrumData &mcnpData, int dataTypeIndex, double *binBounds, int binVecSize)
{
    std::stringstream stream;
    double totalDiff=0.;

    else if(dataTypeIndex==2)
    {
        totalDiff+=CompareHist(stream, nSecMomAngle, mcnpData.GetNSecMomAngle(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==4)
    {
        totalDiff+=CompareHist(stream, nDelMomAngle, mcnpData.GetNDelMomAngle(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==6)
    {
        totalDiff+=CompareHist(stream, pSecMomAnglei, mcnpData.GetPSecMomAnglei(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==9)
    {
        totalDiff+=CompareHist(stream, nSecKEn, mcnpData.GetNSecKEn(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==10)
    {
        totalDiff+=CompareHist(stream, nDelKEn, mcnpData.GetNDelKEn(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==11)
    {
        totalDiff+=CompareHist(stream, pSecKEn, mcnpData.GetPSecKEn(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==12)
    {
        totalDiff+=CompareHist(stream, nSecYield, mcnpData.GetNSecYield(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==13)
    {
        totalDiff+=CompareHist(stream, nDelYield, mcnpData.GetNDelYield(), binBounds, binVecSize);
    }

    else if(dataTypeIndex==14)
    {
        totalDiff+=CompareHist(stream, pSecYield, mcnpData.GetPSecYield(), binBounds, binVecSize);
    }

    SetDataStream( outFileName, stream, false );
    return totalDiff;
}

double FSSpectrumData::CompareHist(std::stringstream &stream, std::vector<double> &g4ndlData, std::vector<double> &mcnpData, double *binLimits, int binVecSize)
{
    double binBounds[numBins+1];
    double *binDim;
    double maxNum, minNum, sumDiff=0.;
    int numBinBound;

    stream.fill(' ');
    stream.precision(6);

    if(binVecSize==-1)
    {
        if(g4ndlData.size()>0)
            minNum=maxNum=g4ndlData[0];
        else if(mcnpData.size()>0)
            minNum=maxNum=mcnpData[0];
        else
        {
            return 0.;
        }

        stream << std::endl;

        for(int i=0; i<int(g4ndlData.size()); i++)
        {
            if(maxNum<g4ndlData[i])
            {
                maxNum=g4ndlData[i];
            }
            if(minNum>g4ndlData[i])
            {
                minNum=g4ndlData[i];
            }
        }

        for(int i=0; i<int(mcnpData.size()); i++)
        {
            if(maxNum<mcnpData[i])
            {
                maxNum=mcnpData[i];
            }
            if(minNum>mcnpData[i])
            {
                minNum=mcnpData[i];
            }
        }

        for(int i=0; i<numBins+1; i++)
        {
            binBounds[i]=(maxNum-minNum)*i/numBins+minNum;
        }
        binDim=binBounds;
        numBinBound=numBins+1;
    }
    else
    {
        binDim=binLimits;
        numBinBound=binVecSize;

        if((g4ndlData.size()==0)&&(mcnpData.size()==0))
        {
            for(int i=0; i<3; i++)
            {
                if(i==0)
                {
                    stream << "G4NDL Hist" << std::endl;
                }
                if(i==1)
                {
                    stream << "MCNP Hist" << std::endl;
                }
                if(i==2)
                {
                    stream << "Diff Sq Hist" << std::endl;
                }
                for(int j=0; j<numBinBound-1; j++)
                {
                    stream << std::setw(14) << std::right << 0.;
                    if(((j+1)%6==0)||(j==numBinBound-2))
                        stream << '\n';
                }
            }
            stream << '\n';
        }
    }

    double *g4ndlHist = new double [numBinBound-1];
    double *mcnpHist = new double [numBinBound-1];
    double *diffHist = new double [numBinBound-1];
    for(int i=0; i<numBinBound-1; i++)
    {
        g4ndlHist[i]=0.;
        mcnpHist[i]=0.;
        diffHist[i]=0.;
    }

    for(int i=0; i<int(g4ndlData.size()); i++)
    {
        for(int j=1; j<numBinBound-1; j++)
        {
            if(g4ndlData[i]<=binDim[j])
            {
                g4ndlHist[j-1]++;
                break;
            }
        }
    }

    if(g4ndlData.size()>0)
    {
        for(int j=0; j<numBins; j++)
        {
            g4ndlHist[j]/=g4ndlData.size();
        }
    }

    for(int i=0; i<int(mcnpData.size()); i++)
    {
        for(int j=1; j<numBinBound-1; j++)
        {
            if(mcnpData[i]<=binDim[j])
            {
                mcnpHist[j-1]++;
                break;
            }
        }
    }

    if(mcnpData.size()>0)
    {
        for(int j=0; j<numBinBound-1; j++)
        {
            mcnpHist[j]/=mcnpData.size();
        }
    }

    if(binVecSize==-1)
    {
        for(int j=0; j<numBinBound; j++)
        {
            stream << std::setw(14) << std::right << binBounds[j];
            if(((j+1)%6==0)||(j==numBinBound-1))
                stream << '\n';
        }
    }

    stream << "G4NDL Hist" << std::endl;
    for(int j=0; j<numBinBound-1; j++)
    {
        stream << std::setw(14) << std::right << g4ndlHist[j];
        if(((j+1)%6==0)||(j==numBinBound-2))
            stream << '\n';
    }

    stream << "MCNP Hist" << std::endl;
    for(int j=0; j<numBinBound-1; j++)
    {
        stream << std::setw(14) << std::right << mcnpHist[j];
        if(((j+1)%6==0)||(j==numBinBound-2))
            stream << '\n';
    }

    stream << "Diff Sq Hist" << std::endl;
    for(int j=0; j<numBinBound-1; j++)
    {
        diffHist[j] = std::pow(mcnpHist[j]-g4ndlHist[j],2);
        sumDiff+=diffHist[j];
        stream << std::setw(14) << std::right << diffHist[j];
        if(((j+1)%6==0)||(j==numBinBound-2))
            stream << '\n';
    }

    if(g4ndlHist)
        delete [] g4ndlHist;
    if(mcnpHist)
        delete [] mcnpHist;
    if(diffHist)
        delete [] diffHist;

    stream << '\n';
    return sumDiff;
}

void FSSpectrumData::GetBinLimits(double &minVal, double &maxVal, int dataTypeIndex, bool &hasData)
{
    if(dataTypeIndex==0)
    {
        minVal=GetMin(nPrimMomDirPhi, hasData);
        maxVal=GetMax(nPrimMomDirPhi);
    }

    else if(dataTypeIndex==1)
    {
        minVal=GetMin(nPrimMomDirTheta, hasData);
        maxVal=GetMax(nPrimMomDirTheta);
    }

    else if(dataTypeIndex==2)
    {
        minVal=GetMin(nSecMomDirPhi, hasData);
        maxVal=GetMax(nSecMomDirPhi);
    }

    else if(dataTypeIndex==3)
    {
        minVal=GetMin(nSecMomDirTheta, hasData);
        maxVal=GetMax(nSecMomDirTheta);
    }

    else if(dataTypeIndex==4)
    {
        minVal=GetMin(nDelMomDirPhi, hasData);
        maxVal=GetMax(nDelMomDirPhi);
    }

    else if(dataTypeIndex==5)
    {
        minVal=GetMin(nDelMomDirTheta, hasData);
        maxVal=GetMax(nDelMomDirTheta);
    }

    else if(dataTypeIndex==6)
    {
        minVal=GetMin(pSecMomDirPhi, hasData);
        maxVal=GetMax(pSecMomDirPhi);
    }

    else if(dataTypeIndex==7)
    {
        minVal=GetMin(pSecMomDirTheta, hasData);
        maxVal=GetMax(pSecMomDirTheta);
    }

    else if(dataTypeIndex==8)
    {
        minVal=GetMin(nPrimKEn, hasData);
        maxVal=GetMax(nPrimKEn);
    }

    else if(dataTypeIndex==9)
    {
        minVal=GetMin(nSecKEn, hasData);
        maxVal=GetMax(nSecKEn);
    }

    else if(dataTypeIndex==10)
    {
        minVal=GetMin(nDelKEn, hasData);
        maxVal=GetMax(nDelKEn);
    }

    else if(dataTypeIndex==11)
    {
        minVal=GetMin(pSecKEn, hasData);
        maxVal=GetMax(pSecKEn);
    }

    else if(dataTypeIndex==12)
    {
        minVal=GetMin(nSecYield, hasData);
        maxVal=GetMax(nSecYield);
    }

    else if(dataTypeIndex==13)
    {
        minVal=GetMin(nDelYield, hasData);
        maxVal=GetMax(nDelYield);
    }

    else if(dataTypeIndex==14)
    {
        minVal=GetMin(pSecYield, hasData);
        maxVal=GetMax(pSecYield);
    }
}

double FSSpectrumData::GetMin(std::vector<double> &valVec, bool &hasData)
{
    if(valVec.size()==0)
    {
        return 0.;
    }

    hasData=true;

    double minVal = valVec[0];
    for(int i=1; i<int(valVec.size()); i++)
    {
        if(valVec[i]<minVal)
            minVal=valVec[i];
    }

    return minVal;
}

double FSSpectrumData::GetMax(std::vector<double> &valVec)
{
    if(valVec.size()==0)
    {
        return 0.;
    }

    double maxVal = valVec[0];
    for(int i=1; i<int(valVec.size()); i++)
    {
        if(valVec[i]>maxVal)
            maxVal=valVec[i];
    }

    return maxVal;
}

void FSSpectrumData::SetDataStream( std::string filename , std::stringstream& ss, bool overWrite )
{
    // Use regular text file
    std::string compfilename(filename);

    if(compfilename.substr((compfilename.length()-2),2)==".z")
    {
        compfilename.erase(compfilename.size()-2, 2);
    }

    std::ofstream *out;
    if(overWrite)
        out = new std::ofstream ( compfilename.c_str() , std::ios::out | std::ios::trunc );
    else
        out = new std::ofstream ( compfilename.c_str() , std::ios::out | std::ios::app );
    if ( ss.good() )
    {
         ss.seekg( 0 , std::ios::end );
         int file_size = ss.tellg();
         ss.seekg( 0 , std::ios::beg );
         char* filedata = new char[ file_size ];
         while ( ss ) {
            ss.read( filedata , file_size );
            if(!file_size)
            {
                std::cout << "\n #### Error the size of the stringstream is invalid ###" << std::endl;
                break;
            }
         }
         out->write(filedata, file_size);
        if (out->fail())
        {
            std::cout << std::endl << "writing the ascii data to the output file " << compfilename << " failed" << std::endl
                 << " may not have permission to delete an older version of the file" << std::endl;
        }

         delete [] filedata;
    }
    else
    {
    // found no data file
    //                 set error bit to the stream
     ss.setstate( std::ios::badbit );

     std::cout << std::endl << "### failed to write to ascii file " << compfilename << " ###" << std::endl;
    }
    out->close();
    delete out;
   ss.str("");
}
