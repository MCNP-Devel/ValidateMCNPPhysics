#ifndef FSSPECTRUMDATA_HH
#define FSSPECTRUMDATA_HH

#include "G4HadFinalState.hh"
#include "G4DynamicParticle.hh"
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#define numBins 100

class FSSpectrumData
{
    public:
        FSSpectrumData();
        virtual ~FSSpectrumData();
        void AddData(G4HadFinalState* result);
        void AddSecondary(G4DynamicParticle &nSec);
        void AddDelayed(G4DynamicParticle &nDel);
        void AddPhoton(G4DynamicParticle &pSec);
        void AddPrimary(G4HadFinalState &nPrim);

        std::vector<double>& GetNSecMomAngle() {return nSecMomAnglei;}
        std::vector<double>& GetNDelMomAngle() {return nDelMomAngle;}
        std::vector<double>& GetPSecMomAngle() {return pSecMomAngle;}
        std::vector<double>& GetNSecKEn() {return nSecKEn;}
        std::vector<double>& GetNDelKEn() {return nDelKEn;}
        std::vector<double>& GetPSecKEn() {return pSecKEn;}
        std::vector<double>& GetNSecYield() {return nSecYield;}
        std::vector<double>& GetNDelYield() {return nDelYield;}
        std::vector<double>& GetPSecYield() {return pSecYield;}

        void SetSecYield(int yield)
        {
            nSecYield.push_back(yield);
        }
        void SetDelYield(int yield)
        {
            nDelYield.push_back(yield);
        }
        void SetPhYield(int yield)
        {
            pSecYield.push_back(yield);
        }
        double CompareFSData(std::string &outFileName, FSSpectrumData &mcnpData, bool *relevantData);
        double CompareFSData(std::string &outFileName, FSSpectrumData &mcnpData, int dataTypeIndex, double *binBounds, int binVecSize);
        double CompareHist(std::stringstream &stream, std::vector<double> &g4ndlData, std::vector<double> &mcnpData, double *binLimits=NULL, int binVecSize=-1);
        void GetBinLimits(double &minVal, double &maxVal, int dataTypeIndex, bool &hasData);
        double GetMin(std::vector<double> &valVec, bool &hasData);
        double GetMax(std::vector<double> &valVec);
        void SetDataStream( std::string filename , std::stringstream& ss, bool overWrite );
        void Clear()
        {
            nSecMomAngle.clear();
            nDelMomAngle.clear();
            pSecMomAngle.clear();
            nSecKEn.clear();
            nDelKEn.clear();
            pSecKEn.clear();
            nSecYield.clear();
            nDelYield.clear();
            pSecYield.clear();
        }
    protected:

    private:
    std::vector<double> nSecMomAngle, nDelMomAngle, pSecMomAngle;
    std::vector<double> nSecKEn, nDelKEn, pSecKEn;
    std::vector<double> nSecYield, nDelYield, pSecYield;
};

#endif // FSSPECTRUMDATA_HH
