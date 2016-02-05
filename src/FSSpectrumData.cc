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
    AddSecondary(result[0]);

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

void FSSpectrumData::AddSecondary(G4HadFinalState &nSec)
{
    if(nSec.GetEnergyChange()>0.)
    {
        nSecMomAngle.push_back(std::cos(std::sqrt(std::pow(nSec.GetMomentumChange().getPhi(),2)+std::pow(nSec.GetMomentumChange().getTheta(),2))));
        nSecKEn.push_back(std::log10(nSec.GetEnergyChange()));
    }
}

void FSSpectrumData::AddSecondary(G4DynamicParticle &nSec)
{
    if(nSec.GetKineticEnergy()>0.)
    {
        nSecMomAngle.push_back(std::cos(std::sqrt(std::pow(nSec.GetMomentumDirection().getPhi(),2)+std::pow(nSec.GetMomentumDirection().getTheta(),2))));
        nSecKEn.push_back(std::log10(nSec.GetKineticEnergy()));
    }
}

void FSSpectrumData::AddDelayed(G4DynamicParticle &nDel)
{
    if(nDel.GetKineticEnergy()>0.)
    {
        nDelMomAngle.push_back(std::cos(std::sqrt(std::pow(nDel.GetMomentumDirection().getPhi(),2)+std::pow(nDel.GetMomentumDirection().getTheta(),2))));
        nDelKEn.push_back(std::log10(nDel.GetKineticEnergy()));
    }
}

void FSSpectrumData::AddPhoton(G4DynamicParticle &pSec)
{
    if(pSec.GetKineticEnergy()>0.)
    {
        pSecMomAngle.push_back(std::cos(std::sqrt(std::pow(pSec.GetMomentumDirection().getPhi(),2)+std::pow(pSec.GetMomentumDirection().getTheta(),2))));
        pSecKEn.push_back(std::log10(pSec.GetKineticEnergy()));
    }
}

double FSSpectrumData::CompareFSData(std::string &outFileName, double **mcnpIsoData, int dataTypeIndex, double *binBounds, int binVecSize)
{
    std::stringstream stream;
    double totalDiff=0.;

    if(dataTypeIndex==0)
    {
        totalDiff+=CompareHist(stream, nSecMomAngle, mcnpIsoData[0], binBounds, binVecSize);
    }

    else if(dataTypeIndex==1)
    {
        totalDiff+=CompareHist(stream, nDelMomAngle, mcnpIsoData[1], binBounds, binVecSize);
    }

    else if(dataTypeIndex==2)
    {
        totalDiff+=CompareHist(stream, pSecMomAngle, mcnpIsoData[2], binBounds, binVecSize);
    }

    else if(dataTypeIndex==3)
    {
        totalDiff+=CompareHist(stream, nSecKEn, mcnpIsoData[3], binBounds, binVecSize);
    }

    else if(dataTypeIndex==4)
    {
        totalDiff+=CompareHist(stream, nDelKEn, mcnpIsoData[4], binBounds, binVecSize);
    }

    else if(dataTypeIndex==5)
    {
        totalDiff+=CompareHist(stream, pSecKEn, mcnpIsoData[5], binBounds, binVecSize);
    }

    else if(dataTypeIndex==6)
    {
        totalDiff+=CompareHist(stream, nSecYield, mcnpIsoData[6], binBounds, binVecSize);
    }

    else if(dataTypeIndex==7)
    {
        totalDiff+=CompareHist(stream, nDelYield, mcnpIsoData[7], binBounds, binVecSize);
    }

    else if(dataTypeIndex==8)
    {
        totalDiff+=CompareHist(stream, pSecYield, mcnpIsoData[8], binBounds, binVecSize);
    }

    SetDataStream( outFileName, stream, false );
    return totalDiff;
}

double FSSpectrumData::CompareHist(std::stringstream &stream, std::vector<double> &g4ndlData, double *mcnpIsoData, double *binLimits, int binVecSize)
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

        if(g4ndlData.size()==0)
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
                for(int j=0; j<numBinBound; j++)
                {
                    stream << std::setw(14) << std::right << 0.;
                    if(((j+1)%6==0)||(j==numBinBound-1))
                        stream << '\n';
                }
            }
            stream << '\n';
        }
    }

    double *g4ndlHist = new double [numBinBound];
    double *diffHist = new double [numBinBound];
    for(int i=0; i<numBinBound; i++)
    {
        g4ndlHist[i]=0.;
        diffHist[i]=0.;
    }

    for(int i=0; i<int(g4ndlData.size()); i++)
    {
        for(int j=0; j<numBinBound; j++)
        {
            if((g4ndlData[i]<=binDim[j])||(j==numBinBound-1))
            {
                g4ndlHist[j]++;
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
    for(int j=0; j<numBinBound; j++)
    {
        stream << std::setw(14) << std::right << g4ndlHist[j];
        if(((j+1)%6==0)||(j==numBinBound-1))
            stream << '\n';
    }

    stream << "MCNP Hist" << std::endl;
    for(int j=0; j<numBinBound; j++)
    {
        stream << std::setw(14) << std::right << mcnpIsoData[j];
        if(((j+1)%6==0)||(j==numBinBound-1))
            stream << '\n';
    }

    stream << "Diff Sq Hist" << std::endl;
    for(int j=0; j<numBinBound; j++)
    {
        diffHist[j] = std::pow(mcnpIsoData[j]-g4ndlHist[j],2);
        sumDiff+=diffHist[j];
        stream << std::setw(14) << std::right << diffHist[j];
        if(((j+1)%6==0)||(j==numBinBound-1))
            stream << '\n';
    }

    if(g4ndlHist)
        delete [] g4ndlHist;
    if(diffHist)
        delete [] diffHist;

    stream << '\n';
    return sumDiff;
}

bool FSSpectrumData::HasData(int dataTypeIndex)
{
    bool hasData=false;
    if(dataTypeIndex==0)
    {
        (nSecMomAngle.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==1)
    {
        (nDelMomAngle.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==2)
    {
        (pSecMomAngle.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==3)
    {
        (nSecKEn.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==4)
    {
        (nDelKEn.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==5)
    {
        (pSecKEn.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==6)
    {
        (nSecYield.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==7)
    {
        (nDelYield.size()==0) ? (hasData=false) : (hasData=true);
    }

    else if(dataTypeIndex==8)
    {
        (pSecYield.size()==0) ? (hasData=false) : (hasData=true);
    }
    return hasData;
}

void FSSpectrumData::GetBinLimits(double &minVal, double &maxVal, int dataTypeIndex, bool &hasData)
{

    if(dataTypeIndex==0)
    {
        minVal=GetMin(nSecMomAngle, hasData);
        maxVal=GetMax(nSecMomAngle);
    }

    else if(dataTypeIndex==1)
    {
        minVal=GetMin(nDelMomAngle, hasData);
        maxVal=GetMax(nDelMomAngle);
    }

    else if(dataTypeIndex==2)
    {
        minVal=GetMin(pSecMomAngle, hasData);
        maxVal=GetMax(pSecMomAngle);
    }

    else if(dataTypeIndex==3)
    {
        minVal=GetMin(nSecKEn, hasData);
        maxVal=GetMax(nSecKEn);
    }

    else if(dataTypeIndex==4)
    {
        minVal=GetMin(nDelKEn, hasData);
        maxVal=GetMax(nDelKEn);
    }

    else if(dataTypeIndex==5)
    {
        minVal=GetMin(pSecKEn, hasData);
        maxVal=GetMax(pSecKEn);
    }

    else if(dataTypeIndex==6)
    {
        minVal=GetMin(nSecYield, hasData);
        maxVal=GetMax(nSecYield);
    }

    else if(dataTypeIndex==7)
    {
        minVal=GetMin(nDelYield, hasData);
        maxVal=GetMax(nDelYield);
    }

    else if(dataTypeIndex==8)
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
