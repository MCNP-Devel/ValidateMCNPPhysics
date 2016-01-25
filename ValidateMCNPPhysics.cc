#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <dirent.h>
#include <iomanip>

//Models
#include "globals.hh"
#include "G4NeutronHPChannel.hh"
#include "G4HadronicInteraction.hh"
#include "G4NeutronHPChannelList.hh"
#include "include/VFDNeutronHPElasticFS.hh"
#include "include/VFDNeutronHP2AInelasticFS.hh"
#include "include/VFDNeutronHP2N2AInelasticFS.hh"
#include "include/VFDNeutronHP2NAInelasticFS.hh"
#include "include/VFDNeutronHP2NDInelasticFS.hh"
#include "include/VFDNeutronHP2NInelasticFS.hh"
#include "include/VFDNeutronHP2NPInelasticFS.hh"
#include "include/VFDNeutronHP2PInelasticFS.hh"
#include "include/VFDNeutronHP3AInelasticFS.hh"
#include "include/VFDNeutronHP3NAInelasticFS.hh"
#include "include/VFDNeutronHP3NInelasticFS.hh"
#include "include/VFDNeutronHP3NPInelasticFS.hh"
#include "include/VFDNeutronHP4NInelasticFS.hh"
#include "include/VFDNeutronHPAInelasticFS.hh"
#include "include/VFDNeutronHPD2AInelasticFS.hh"
#include "include/VFDNeutronHPDAInelasticFS.hh"
#include "include/VFDNeutronHPDInelasticFS.hh"
#include "include/VFDNeutronHPHe3InelasticFS.hh"
#include "include/VFDNeutronHPN2AInelasticFS.hh"
#include "include/VFDNeutronHPN2PInelasticFS.hh"
#include "include/VFDNeutronHPN3AInelasticFS.hh"
#include "include/VFDNeutronHPNAInelasticFS.hh"
#include "include/VFDNeutronHPND2AInelasticFS.hh"
#include "include/VFDNeutronHPNDInelasticFS.hh"
#include "include/VFDNeutronHPNHe3InelasticFS.hh"
#include "include/VFDNeutronHPNInelasticFS.hh"
#include "include/VFDNeutronHPNPAInelasticFS.hh"
#include "include/VFDNeutronHPNPInelasticFS.hh"
#include "include/VFDNeutronHPNT2AInelasticFS.hh"
#include "include/VFDNeutronHPNTInelasticFS.hh"
#include "include/VFDNeutronHPNXInelasticFS.hh"
#include "include/VFDNeutronHPPAInelasticFS.hh"
#include "include/VFDNeutronHPPDInelasticFS.hh"
#include "include/VFDNeutronHPPInelasticFS.hh"
#include "include/VFDNeutronHPPTInelasticFS.hh"
#include "include/VFDNeutronHPT2AInelasticFS.hh"
#include "include/VFDNeutronHPTInelasticFS.hh"
#include "include/VFDNeutronHPCaptureFS.hh"
#include "include/VFDNeutronHPFSFissionFS.hh"
#include "include/VFDNeutronHPFissionBaseFS.hh"
#include "include/VFDNeutronHPFinalState.hh"

//Data Containers
#include "include/IsotopeMass.hh"
#include "include/FSSpectrumData.hh"
#include "include/ElementNames.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTable.hh"
#include "topc.h"
#include "include/TaskInput.hh"
#include "include/TaskOutput.hh"
#include "include/MarshaledTaskInput.h"
#include "include/MarshaledTaskOutput.h"
#include "include/MarshaledObj.h"
#include "include/VFDIonTable.hh"

#define SetNumEnerFlag 1
#define numIncEner 261
#define numPointPerOrder 20
#define numDataTypes 9
#define numFS 36
#define numSamples 600

#ifndef TOPC_USE
#define TOPC_USE 0
#endif

// need to repeat FS application multiple times for the same incoming energy, for processes that use sampling and random techniques (ie everything except elastic) and then compare
// distributions, this can easily be done by setting the yield to 50 for these reactions, while seperately extracting the real yield using new function
// or just reapeat the process multiple timesS
// could get rid of unnessary code in FS classes to speed up calculation,  no need to calc dopp broad nuc or use so many data containers, could make our own data container that stores all of the data we care about
// add in fission/ff/ and capture/FSMF6

// we changed the max tries from 1000 to 50 in VFDNeutronHPPhotonDist.cc:372 to improve efficiency since if it doesn't work after 50 tries it probably won't work for a 1000
// we took out the recoiling isotope as a secondary for efficiency sake since we don't look at it, to add back in search for it using GetIon key word
// we filter out photons with zero kinetic energy from the results
// we got rid of photon and delayed data for fission and inelastic reaction since we want to just focus on checking that the neutronics is correct for now

// 080808 Something unexpected is happen in G4NeutronHPLabAngularEnergy is caused by the energy regime of the angular-energy not extending down to the minimum energy supplied in the cross-section data, this is not our fault
// this is the fault of the G4NDL data, and seems to mainly occur for the isotopes of berylium
using namespace std;

void ExtractDir(G4String dirName, VFDNeutronHPFinalState ***isoFSData, string **isoNameList, int &vecSize, VFDNeutronHPFinalState *fsType, std::pair <double,double> **enerBound);
#if TOPC_USE
    double CompareData(string outDir, int fsSizeG4NDL, int fsSizeMCNP, bool sampling);
#else
    double CompareData(string outDir, string *isoNameListG4NDL, VFDNeutronHPFinalState **isoFSDataG4NDL, int fsSizeG4NDL, std::pair <double,double> *enerBoundG4NDL, string *isoNameListMCNP, VFDNeutronHPFinalState **isoFSDataMCNP, int fsSizeMCNP, std::pair <double,double> *enerBoundMCNP, bool sampling, bool *relevantData);
#endif
double CalcDiff(string outFileName, VFDNeutronHPFinalState *isoFSDataG4NDL, std::pair <double,double> enerBoundG4NDL, VFDNeutronHPFinalState *isoFSDataMCNP, std::pair <double,double> enerBoundMCNP, bool sampling, bool *relevantData);
void ExtractZA(string fileName, int &Z, int &A);
void ExtractZA(string fileName, string &Z, string &A);
bool DirectoryExists( const char* pzPath );
void GetDataStream( string filename, std::stringstream& ss);
void SetDataStream( string filename , std::stringstream& ss, bool overWrite );

double incEnerVec[numIncEner];
int incEnerSize;
//G4Track *nTrack;
G4DynamicParticle *nParticle;
bool master;

string dataTypeHeaders[numDataTypes] = {"Comparing the secondary neutron deflection angle", "Comparing the delayed neutron deflection angle", "Comparing the secondary photon deflection angle",
                                    "Comparing the secondary neutron kinetic energy", "Comparing the delayed neutron kinetic energy", "Comparing the secondary photon kinetic energy",
                                    "Comparing the secondary neutron yield", "Comparing the delayed neutron yield", "Comparing the secondary photon yield"};


#if TOPC_USE
TaskOutput *taskOut;
MarshaledTaskOutput *taskOutMarsh;
MarshaledTaskInput *taskInMarsh;
TaskInput *taskIn;

int fsIndexGbl;
string outDirGbl, inMCNPFSDirGbl;
string *isoNameListG4NDLGbl;
VFDNeutronHPFinalState **isoFSDataG4NDLGbl;
std::pair<double,double> *enerBoundG4NDLGbl;
bool *relevantDataGbl[numFS];
double *sumErrorVecGbl;
double binBounds[3][numBins];

string fsDirNameListGbl[numFS] = {"Elastic/FS/", "Capture/FS/", "Fission/FS/", "Fission/FC/", "Fission/SC/", "Fission/TC/", "Fission/LC/", "Inelastic/F01/",
                                    "Inelastic/F02/", "Inelastic/F04/", "Inelastic/F05/", "Inelastic/F06/", "Inelastic/F07/", "Inelastic/F08/", "Inelastic/F09/",
                                    "Inelastic/F10/", "Inelastic/F11/", "Inelastic/F13/", "Inelastic/F14/", "Inelastic/F15/", "Inelastic/F18/", "Inelastic/F19/",
                                    "Inelastic/F21/", "Inelastic/F22/", "Inelastic/F23/", "Inelastic/F24/", "Inelastic/F25/", "Inelastic/F26/", "Inelastic/F27/",
                                    "Inelastic/F28/", "Inelastic/F30/", "Inelastic/F31/", "Inelastic/F33/", "Inelastic/F34/", "Inelastic/F35/","Inelastic/F36/"};

//CalcDiff
//This function is run by the slave when a task is assigned to it
//this function pushes neutrons with varing energies throught the physics list specified by the directory


TOPC_BUF CalcDiff(void *input)
{
    // TaskInput has j, k, sampling
    MarshaledTaskInput inputTaskMarsh(input);
    TaskInput *inputTask = inputTaskMarsh.unmarshal();
    char line[256];
	double dummy;
	double mcnpIsoData[numIncEner][numDataTypes][numBins];
    bool relevantDataMCNP[numDataTypes];

    int j = inputTask->g4ndlIndex;
    bool sampling = inputTask->sampling;

    if(inputTask)
        delete inputTask;

    string outFileName = outDirGbl+fsDirNameListGbl[fsIndexGbl]+isoNameListG4NDLGbl[j]+".txt";
    string inMCNPFileName = inMCNPFSDirGbl+fsDirNameListGbl[fsIndexGbl]+isoNameListG4NDLGbl[j].substr(0,isoNameListG4NDLGbl[j].find_last_of('_'))+".txt";

    stringstream inMCNPFile;
    GetDataStream(inMCNPFileName, inMCNPFile);

    stringstream stream;

    stream.getline(line,256);
    stream.getline(line,256);
    stream.getline(line,256);
    stream.getline(line,256);

    for(int i=0; i<numDataTypes; i++)
    {
        stream >> relevantDataMCNP[i];
    }

	for (int i = 0; i<numIncEner; i++)
	{
		stream >> dummy;
	}

	for (int i = 0; i<numDataTypes; i++)
	{
		stream.getline(line, 256);
		stream.getline(line, 256);
		stream.getline(line, 256);
		for (int k = 0; k < numBins; k++)
		{
			stream >> dummy;
		}

		if (relevantDataMCNP[i] == 1)
		{
			for (int j = 0; j < numIncEner; j++)
			{
				for (int k = 0; k < numBins; k++)
				{
					stream >> dummy;
					mcnpIsoData[j][i][k] = dummy;
				}
			}
		}
	}

    G4HadProjectile *projectile;
    int points=1;
    double totalDiff=0.;
    FSSpectrumData *fsDataG4NDL;
    G4HadFinalState *resultG4NDLVec;

    fsDataG4NDL = new FSSpectrumData [incEnerSize];

    stream << "Comparison of G4NDL Data to the Converted MCNP data using\n" << incEnerSize << " incoming energy points\n" << numBins << " bins per histogram\n" << std::endl;
    for(int i=0; i<numDataTypes; i++)
    {
        stream << relevantDataGbl[fsIndexGbl][i] << '\t';
    }
    stream << '\n' << endl;

    stream << std::setprecision(6);
    for(int i=0; i<incEnerSize; i++)
    {
        stream << std::setw(14) << std::right << incEnerVec[i];
        if(((i+1)%6==0)||(i==incEnerSize-1))
            stream << '\n';
    }

    SetDataStream( outFileName, stream, true);


    for(int i=0; i<incEnerSize; i++)
    {
        nParticle->SetKineticEnergy(incEnerVec[i]);
        projectile = new G4HadProjectile(*nParticle);

        if(sampling)
        {
            points=numSamples;
        }
        for(int l=0; l<points; l++)
        {
            resultG4NDLVec = isoFSDataG4NDLGbl[j]->ApplyYourself(*projectile);
            if((resultG4NDLVec==0))
            {
                stream.clear();
                stream.str("");
                stream << "No final state data to compare!" << std::endl;
                SetDataStream( outFileName, stream, false);
                delete projectile;

                if(taskOut)
                    delete taskOut;
                taskOut = new TaskOutput(totalDiff/incEnerSize);

                if(taskOutMarsh)
                    delete taskOutMarsh;
                taskOutMarsh = new MarshaledTaskOutput(taskOut);

                return TOPC_MSG( taskOutMarsh->getBuffer(), taskOutMarsh->getBufferSize());
            }

            fsDataG4NDL[i].AddData(resultG4NDLVec);
        }
        delete projectile;
    }
    bool hasData;
    int boundType;
    for(int h=0; h<numDataTypes; h++)
    {
        if((relevantDataGbl[fsIndexGbl][h])&&(relevantDataMCNP[h]))
        {
            hasData=false;
            stream.clear();
            stream.str("");
            stream << '\n' << dataTypeHeaders[h] << '\n' << std::endl;

            if(!hasData)
            {
                continue;
            }

            boundType=floor(h/3);

            for(int i=0; i<numBins; i++)
            {
                stream << std::setw(14) << std::right << binBounds[boundType][i];
                if(((i+1)%6==0)||(i==numBins-1))
                    stream << '\n';
            }
            SetDataStream( outFileName, stream, false);

            for(int i=0; i<incEnerSize; i++)
            {
                totalDiff+=fsDataG4NDL[i].CompareFSData(outFileName, mcnpIsoData[i], h, binBounds[boundType], numBins);
            }
        }
    }

    delete [] fsDataG4NDL;
    delete isoFSDataG4NDLGbl[j];
    isoFSDataG4NDLGbl[j]=NULL;

    if(taskOut)
        delete taskOut;
    taskOut = new TaskOutput(totalDiff/incEnerSize);

    if(taskOutMarsh)
        delete taskOutMarsh;
    taskOutMarsh = new MarshaledTaskOutput(taskOut);

    return TOPC_MSG( taskOutMarsh->getBuffer(), taskOutMarsh->getBufferSize());
}

//CheckProgress
//this function is run by the master process after a slave has finished running DoBroadening()
//This takes in the results of the task assigned to slave and stores it in the master
TOPC_ACTION CheckProgress(void* input, void *output)
{
    MarshaledTaskInput inputTaskMarsh(input);
    TaskInput* inputTask = inputTaskMarsh.unmarshal();
    int j = inputTask->g4ndlIndex;

    MarshaledTaskOutput outputTaskMarsh(output);
    TaskOutput* outputTask = outputTaskMarsh.unmarshal();
    sumErrorVecGbl[j]=outputTask->totalDiff;

    if(inputTask)
        delete inputTask;

    if(outputTask)
        delete outputTask;

    return NO_ACTION;
}
#endif


int main (int argc, char** argv)
{

#if TOPC_USE
    TOPC_OPT_trace = 0;
	TOPC_OPT_slave_timeout=6048000;
    TOPC_init(&argc, &argv);
    master = bool(TOPC_is_master());
    //int processID = TOPC_rank();
#else
    master=true;
#endif

    CLHEP::HepRandom *CLHEPRand = new CLHEP::HepRandom();

    CLHEP::RanecuEngine *theEngine = new CLHEP::RanecuEngine(2013092304);

    CLHEPRand->setTheEngine(theEngine);

    CLHEPRand->showEngineStatus();
    CLHEPRand->saveEngineStatus();

    G4Neutron::NeutronDefinition();
    G4Geantino::GeantinoDefinition();
    G4Deuteron::DeuteronDefinition();
    G4Triton::TritonDefinition();
    G4Gamma::GammaDefinition();
    G4He3::He3Definition();
    G4Proton::ProtonDefinition();
    G4AntiProton::AntiProtonDefinition();
    G4PionPlus::PionPlusDefinition();
    G4PionZero::PionZeroDefinition();
    G4PionMinus::PionMinusDefinition();
    G4KaonPlus::KaonPlusDefinition();
    G4KaonZero::KaonZeroDefinition();
    G4KaonMinus::KaonMinusDefinition();
    G4KaonZeroShort::KaonZeroShortDefinition();
    G4KaonZeroLong::KaonZeroLongDefinition();
    G4Alpha::AlphaDefinition();
    G4Electron::ElectronDefinition();
    G4Positron::PositronDefinition();
    G4GenericIon *gion = G4GenericIon::Definition();

    G4ParticleTable* pTable = G4ParticleTable::GetParticleTable();
    G4IonTable *oldIonTable=pTable->fIonTable;
    //VFDIonTable newIonTable = VFDIonTable(oldIonTable);
    //pTable->fIonTable= dynamic_cast<G4IonTable*>(&newIonTable);

    G4ProcessManager gProcessMan(gion);
    gion->SetProcessManager(&gProcessMan);
    oldIonTable->InitializeLightIons();
    oldIonTable->PreloadNuclide();
    pTable->SetReadiness(true);
    G4ParticleDefinition *neutron = pTable->FindParticle("neutron");

    nParticle = new G4DynamicParticle(neutron,G4ThreeVector(0.,0.,1.), 10.);

    //nTrack = new G4Track(apValueDynamicParticle, 0., G4ThreeVector(0.,0.,0.));

    int fsSizeG4NDL=0, fsSizeMCNP=0;
    double tempDiff, binInc, incNEner, energyStep;
    vector<double> sumDiffVec;
    string G4NDLDir, compareCS="false", dirName, temp;
    ElementNames* elementNames;
    IsotopeMass* isoMasses;
    elementNames->SetElementNames();
    isoMasses->SetIsotopeMass();

#if TOPC_USE==0
    string fsDirNameListGbl[numFS] = {"Elastic/FS/", "Capture/FS/", "Fission/FS/", "Fission/FC/", "Fission/SC/", "Fission/TC/", "Fission/LC/", "Inelastic/F01/",
                                    "Inelastic/F02/", "Inelastic/F04/", "Inelastic/F05/", "Inelastic/F06/", "Inelastic/F07/", "Inelastic/F08/", "Inelastic/F09/",
                                    "Inelastic/F10/", "Inelastic/F11/", "Inelastic/F13/", "Inelastic/F14/", "Inelastic/F15/", "Inelastic/F18/", "Inelastic/F19/",
                                    "Inelastic/F21/", "Inelastic/F22/", "Inelastic/F23/", "Inelastic/F24/", "Inelastic/F25/", "Inelastic/F26/", "Inelastic/F27/",
                                    "Inelastic/F28/", "Inelastic/F30/", "Inelastic/F31/", "Inelastic/F33/", "Inelastic/F34/", "Inelastic/F35/","Inelastic/F36/"};
    int fsIndexGbl;
    string outDirGbl;
    string *isoNameListG4NDLGbl, *isoNameListMCNPGbl;
    VFDNeutronHPFinalState **isoFSDataG4NDLGbl, **isoFSDataMCNPGbl;
    std::pair<double,double> *enerBoundG4NDLGbl, *enerBoundMCNPGbl;
    bool *relevantDataGbl[numDataTypes];
#endif

    string fileName;
    bool relDataTemp[5][numDataTypes] = {{1, 0, 0, 1, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {1, 1, 0, 1, 1, 0, 1, 1, 0}, {1, 0, 0, 1, 0, 0, 1, 0, 0},{1, 0, 0, 1, 0, 0, 1, 0, 0}};
    relevantDataGbl[0] = relDataTemp[0];
    relevantDataGbl[1] = relDataTemp[1];
    relevantDataGbl[2] = relDataTemp[2];
    for(int i=3; i<7; i++)
    {
        relevantDataGbl[i] = relDataTemp[3];
    }
    for(int i=7; i<numFS; i++)
    {
        relevantDataGbl[i] = relDataTemp[4];
    }

	incNEner = 1.0E-05;
	for (int i = 0; i < 13; i++)
	{
		if (10 * incNEner > 20.0)
			energyStep = (20.0 - incNEner) / numPointPerOrder;
		else
			energyStep = (9 * incNEner) / numPointPerOrder;
		for (int j = 0; j < numPointPerOrder; j++)
		{
			incEnerVec[(i - 1)*numPointPerOrder + j] = incNEner;
			incNEner = incNEner + energyStep;
		}
	}

	incEnerVec[numIncEner] = incNEner;

    binInc = 2*3.1415927/(numBins);
    binBounds[0][0] = -3.1415927+binInc;
    for(int i=1; i<numBins; i++)
    {
        binBounds[0][i]=binBounds[0][i-1]+binInc;
    }

    binInc = (log10(20.0)+6)/(numBins-1);
    binBounds[0][0] = -6;
    for(int i=1; i<numBins; i++)
    {
        binBounds[0][i]=binBounds[0][i-1]+binInc;
    }

    binInc = 20.0/(numBins-1);
    binBounds[0][0] = 0;
    for(int i=1; i<numBins; i++)
    {
        binBounds[0][i]=binBounds[0][i-1]+binInc;
    }

    VFDNeutronHPFinalState* fsTypeList[numFS] = {new VFDNeutronHPElasticFS, new VFDNeutronHPCaptureFS, new VFDNeutronHPFSFissionFS, new VFDNeutronHPFissionBaseFS, new VFDNeutronHPFissionBaseFS
    , new VFDNeutronHPFissionBaseFS, new VFDNeutronHPFissionBaseFS, new VFDNeutronHPNInelasticFS, new VFDNeutronHPNXInelasticFS
    , new VFDNeutronHP2NInelasticFS, new VFDNeutronHP3NInelasticFS, new VFDNeutronHPNAInelasticFS, new VFDNeutronHPN3AInelasticFS
    , new VFDNeutronHP2NAInelasticFS, new VFDNeutronHP3NAInelasticFS, new VFDNeutronHPNPInelasticFS, new VFDNeutronHPN2AInelasticFS
    , new VFDNeutronHPNDInelasticFS, new VFDNeutronHPNTInelasticFS, new VFDNeutronHPNHe3InelasticFS
    , new VFDNeutronHP4NInelasticFS, new VFDNeutronHP2NPInelasticFS, new VFDNeutronHPN2PInelasticFS, new VFDNeutronHPNPAInelasticFS
    , new VFDNeutronHPPInelasticFS, new VFDNeutronHPDInelasticFS, new VFDNeutronHPTInelasticFS, new VFDNeutronHPHe3InelasticFS, new VFDNeutronHPAInelasticFS
    , new VFDNeutronHP2AInelasticFS, new VFDNeutronHP2PInelasticFS, new VFDNeutronHPPAInelasticFS
    , new VFDNeutronHPT2AInelasticFS, new VFDNeutronHPPDInelasticFS, new VFDNeutronHPPTInelasticFS, new VFDNeutronHPDAInelasticFS };

    //bool compCS=0;
    stringstream stream;

    #if TOPC_USE
        if(argc == 2)
        {
            stream << argv[1];
            stream >> fileName;
            stream.str("");
            stream.clear();
        }
        else
        {
            for(int i=0; i<numFS; i++)
            {
                delete fsTypeList[i];
            }
            delete nParticle;
            delete CLHEPRand;
            delete theEngine;
            TOPC_finalize();
            cout << "Please input a file containing the G4NDL data directory, the directory contianing the converted MCNP data,\n the output directory, and the compare CS data flag" << endl;
            return 0;
        }
        GetDataStream(fileName, stream);
        stream >> G4NDLDir >> inMCNPFSDirGbl >> outDirGbl;
    #else
        if(argc == 5)
        {
            stream << ' ' << argv[1] << ' ' << argv[2] << ' ' << argv[3] << ' ' << argv[4];
            stream >> G4NDLDir >> inMCNPFSDirGbl >> outDirGbl >> compareCS;
        }
        if(argc == 4)
        {
            stream << ' ' << argv[1] << ' ' << argv[2] << ' ' << argv[3];
            stream >> G4NDLDir >> inMCNPFSDirGbl >> outDirGbl;
        }
        else
        {
            cout << "Please input the G4NDL data directory, the directory contianing the converted MCNP data, the output directory, and the compare CS data flag" << endl;
            for(int i=0; i<numFS; i++)
            {
                delete fsTypeList[i];
            }
            delete nParticle;
            delete CLHEPRand;
            delete theEngine;
            return 0;
        }
    #endif

//    if(compareCS=="true"||compareCS=="1"||compareCS=="True")
//        compCS=1;

    if(G4NDLDir[G4NDLDir.size()-1]!='/')
        G4NDLDir.push_back('/');

    if(inMCNPFSDirGbl[inMCNPFSDirGbl.size()-1]!='/')
        inMCNPFSDirGbl.push_back('/');

    if(outDirGbl[outDirGbl.size()-1]!='/')
        outDirGbl.push_back('/');

    for(fsIndexGbl=0; fsIndexGbl<numFS; fsIndexGbl++)
    {
        if(master)
            cout << "\n\nExtracting Directory " << G4NDLDir+fsDirNameListGbl[fsIndexGbl] << endl;
        ExtractDir(G4NDLDir+fsDirNameListGbl[fsIndexGbl], &isoFSDataG4NDLGbl, &isoNameListG4NDLGbl, fsSizeG4NDL, fsTypeList[fsIndexGbl], &enerBoundG4NDLGbl);
        if(master)
            cout << "\nComparing FS Data " << G4NDLDir+fsDirNameListGbl[fsIndexGbl] << endl;

        #if TOPC_USE
        sumDiffVec.push_back(CompareData(outDirGbl+fsDirNameListGbl[fsIndexGbl], fsSizeG4NDL, fsSizeMCNP, true));
        #else
        sumDiffVec.push_back(CompareData(outDirGbl+fsDirNameListGbl[fsIndexGbl], isoNameListG4NDLGbl, isoFSDataG4NDLGbl, fsSizeG4NDL, enerBoundG4NDLGbl,
                                            isoNameListMCNPGbl, isoFSDataMCNPGbl, fsSizeMCNP, enerBoundMCNPGbl, true, relevantDataGbl[fsIndexGbl]));
        #endif

        if(master)
        {
            string command = "chmod -R 755 "+outDirGbl+fsDirNameListGbl[fsIndexGbl];
            system(command.c_str());
        }

        for(int j=0; j<fsSizeG4NDL; j++)
        {
            if(isoFSDataG4NDLGbl[j])
                delete isoFSDataG4NDLGbl[j];
        }
        if(fsSizeG4NDL>0)
        {
            delete [] isoFSDataG4NDLGbl;
            delete [] isoNameListG4NDLGbl;
            delete [] enerBoundG4NDLGbl;
        }
    }
    for(int i=0; i<numFS; i++)
    {
        delete fsTypeList[i];
    }
    //delete nTrack;
    delete nParticle;

    for(int i=0; i<int(sumDiffVec.size()); i++)
    {
        for(int j=i; j<int(sumDiffVec.size()); j++)
        {
            if(sumDiffVec[i]<sumDiffVec[j])
            {
                tempDiff=sumDiffVec[i];
                temp=fsDirNameListGbl[i];
                sumDiffVec[i]=sumDiffVec[j];
                fsDirNameListGbl[i]=fsDirNameListGbl[j];
                sumDiffVec[j]=tempDiff;
                fsDirNameListGbl[j]=temp;
            }
        }
    }
    stream.str("");
    stream.clear();
    stream << "The FS directories, ordered from greatest to lowest discrepancy\n" << std::endl;

    for(int i=0; i<int(sumDiffVec.size()); i++)
    {
        stream << fsDirNameListGbl[i] << " has summed difference squared of " << sumDiffVec[i] << std::endl;
    }
    SetDataStream( outDirGbl+"FileDiffList.txt", stream, true);

    elementNames->ClearStore();
    isoMasses->ClearStore();
    delete CLHEPRand;
    delete theEngine;

    if(master)
    {
        system( ("chmod -R  777 "+outDirGbl).c_str());
    }

#if TOPC_USE
    TOPC_finalize();
#endif
    return 0;
}

void ExtractDir(G4String dirName, VFDNeutronHPFinalState ***isoFSData, string **isoNameList, int &vecSize, VFDNeutronHPFinalState *fsType, std::pair <double,double> **enerBound)
{
    int i, j;
    int Z, A, pos1, pos2;
    bool ahead;
    string isoName;
    vector<string> isoNameVec;
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dirName.c_str())) != NULL)
    {
        while ((ent = readdir (dir)) != NULL)
        {
            if((string(ent->d_name)!="..")&&(string(ent->d_name)!="."))
            {
                isoName = ent->d_name;
                ExtractZA(isoName, Z, A);
                if((Z!=-1)&&(A!=-1))
                {
                    ahead=false;
                    for(i=0; i<int(isoNameVec.size()); i++)
                    {
                        for(j=0; j<int(min(isoName.length(), isoNameVec[i].length())); j++)
                        {
                            if(isoName[j]>isoNameVec[i][j])
                            {
                                break;
                            }
                            if(isoName[j]<isoNameVec[i][j])
                            {
                                ahead=true;
                                break;
                            }
                            if((isoName.length()<isoNameVec[i].length())&&(j==int(isoName.length()-1)))
                            {
                                ahead=true;
                            }
                        }
                        if(ahead)
                        {
                            break;
                        }
                    }
                    isoNameVec.insert(isoNameVec.begin()+i,isoName);
                }

            }
        }
    }

    closedir(dir);
    pos1 = dirName.find_last_of('/');
    pos2 = dirName.find_last_of('/', pos1-1);
    G4String theFSType = dirName.substr(pos2+1,pos1-pos2);
    dirName = dirName.substr(0,pos2+1);

    vecSize=isoNameVec.size();
    isoFSData[0] = new VFDNeutronHPFinalState *[vecSize];
    isoNameList[0] = new string [vecSize];
    enerBound[0] = new std::pair<double,double> [vecSize];

    IsotopeMass* isoMasses;
    int total;
    double minE, maxE, dummy, dummyE;
    stringstream stream;

    for(i=0; i<int(isoNameVec.size()); i++)
    {
        isoFSData[0][i] = fsType->New();
        isoNameList[0][i] = isoNameVec[i];
        ExtractZA(isoNameList[0][i], Z, A);
        isoFSData[0][i]->Init(A, Z, 0, dirName, theFSType);

        if(isoFSData[0][i]->HasXsec())
        {
            enerBound[0][i].first = (1+1.0e-5)*isoFSData[0][i]->GetMinEnergy();
            enerBound[0][i].second = (1-1.0e-5)*isoFSData[0][i]->GetMaxEnergy();
        }
        else
        {

            GetDataStream(dirName+"CrossSection/"+isoNameList[0][i], stream);
            stream >> dummy >> dummy >> total;
            stream >> minE >> dummy;
            maxE=minE;
            for(j=0; j<total-1; j++)
            {
                stream >> dummyE >> dummy;
                if(dummy!=0.)
                {
                    maxE=dummyE;
                }
            }

            enerBound[0][i].first = (1+1.0e-5)*minE/1.0e+06;
            enerBound[0][i].second = (1-1.0e-5)*maxE/1.0e+06;

            stream.str("");
            stream.clear();
        }

        VFDNeutronHPFissionBaseFS *fisBaseFS = dynamic_cast<VFDNeutronHPFissionBaseFS*>(isoFSData[0][i]);
        if(fisBaseFS)
        {
            fisBaseFS->SetMass(isoMasses->GetIsotopeMass(Z,A));
        }
    }
    isoNameVec.clear();
}

#if TOPC_USE
double CompareData(string outDir, int fsSizeG4NDL, bool sampling)
{
    double sumDiff=0.;
    TOPC_raw_begin_master_slave( CalcDiff, CheckProgress, NULL );

    double tempDiff;
    string temp;
    stringstream stream;

    if(master)
    {
        bool match;
        sumErrorVecGbl = new double[fsSizeG4NDL];
        for(int i=0; i<fsSizeG4NDL; i++)
        {
            sumErrorVecGbl[i]=-1;
        }
        string zG4NDL, aG4NDL;
        stream << "The FS files that have been compared in this directory, ordered from greatest to lowest discrepancy\n" << std::endl;

        if(!(DirectoryExists((outDir).c_str())))
        {
            system( ("mkdir -p -m=666 "+outDir).c_str());
            if(DirectoryExists((outDir).c_str()))
            {
                cout << "created directory " << outDir << "\n" << endl;
            }
            else
            {
                cout << "\nError: could not create directory " << outDir << "\n" << endl;
                return -1;
            }
        }

        cout << "[";
        cout.flush();
        int nBar=0;
        for(int j=0; j<fsSizeG4NDL; j++)
        {
//           
			if(taskIn)
				delete taskIn;

			taskIn = new TaskInput(j, k, sampling);

			if (taskInMarsh)
				delete taskInMarsh;

			taskInMarsh = new MarshaledTaskInput(taskIn);

			TOPC_raw_submit_task_input(TOPC_MSG(taskInMarsh->getBuffer(), taskInMarsh->getBufferSize()));

            while(nBar<int(j*80/fsSizeG4NDL))
            {
                nBar++;
                cout << "#";
                cout.flush();
            }
        }
        cout << "#]" << endl;
    }

    TOPC_raw_end_master_slave();

    if(master)
    {
        for(int i=0; i<fsSizeG4NDL; i++)
        {
            for(int j=i; j<fsSizeG4NDL; j++)
            {
                if(sumErrorVecGbl[i]<sumErrorVecGbl[j])
                {
                    tempDiff=sumErrorVecGbl[i];
                    temp=isoNameListG4NDLGbl[i];
                    sumErrorVecGbl[i]=sumErrorVecGbl[j];
                    isoNameListG4NDLGbl[i]=isoNameListG4NDLGbl[j];
                    sumErrorVecGbl[j]=tempDiff;
                    isoNameListG4NDLGbl[j]=temp;
                }
            }
        }

        for(int i=0; i<fsSizeG4NDL; i++)
        {
            if(sumErrorVecGbl[i]>-1)
            {
                stream << isoNameListG4NDLGbl[i] << " has summed difference squared of " << sumErrorVecGbl[i] << std::endl;
                sumDiff+=sumErrorVecGbl[i];
            }
        }
        if(sumErrorVecGbl)
            delete [] sumErrorVecGbl;
        SetDataStream( outDir+"FileDiffList.txt", stream, true);
        if(fsSizeG4NDL!=0)
        {
            sumDiff/=fsSizeG4NDL;
        }
    }

    return sumDiff;
}

#else
double CompareData(string outDir, string *isoNameListG4NDL, VFDNeutronHPFinalState **isoFSDataG4NDL, int fsSizeG4NDL, std::pair <double,double> *enerBoundG4NDL,
                    string *isoNameListMCNP, VFDNeutronHPFinalState **isoFSDataMCNP, int fsSizeMCNP, std::pair <double,double> *enerBoundMCNP, bool sampling, bool *relevantData)
{
    double sumDiff=0.;
    double* sumErrorVecGbl;

    double tempDiff;
    string temp;
    stringstream stream;

    if(master)
    {
        bool match;
        sumErrorVecGbl = new double[fsSizeG4NDL];
        for(int i=0; i<fsSizeG4NDL; i++)
        {
            sumErrorVecGbl[i]=-1;
        }
        string zG4NDL, zMCNP, aG4NDL, aMCNP;
        stream << "The FS files that have been compared in this directory, ordered from greatest to lowest discrepancy\n" << std::endl;

        if(!(DirectoryExists((outDir).c_str())))
        {
            system( ("mkdir -p -m=666 "+outDir).c_str());
            if(DirectoryExists((outDir).c_str()))
            {
                cout << "created directory " << outDir << "\n" << endl;
            }
            else
            {
                cout << "\nError: could not create directory " << outDir << "\n" << endl;
                return -1;
            }
        }

        cout << "[";
        cout.flush();
        int nBar=0;
        for(int j=0; j<fsSizeG4NDL; j++)
        {
            match=false;
            ExtractZA(isoNameListG4NDL[j], zG4NDL, aG4NDL);
            if(j<fsSizeMCNP)
            {
                for(int k=j; k<fsSizeMCNP; k++)
                {
                    ExtractZA(isoNameListMCNP[k], zMCNP, aMCNP);
                    if((zG4NDL==zMCNP)&&(aG4NDL==aMCNP))
                    {
                        sumErrorVecGbl[j]=(CalcDiff(outDir+isoNameListG4NDL[j]+".txt", isoFSDataG4NDL[j], enerBoundG4NDL[j], isoFSDataMCNP[k], enerBoundMCNP[k], sampling, relevantData));
                        match=true;
                        break;
                    }
                }
                if(!match)
                {
                    for(int k=0; k<j; k++)
                    {
                        ExtractZA(isoNameListMCNP[k], zMCNP, aMCNP);
                        if((zG4NDL==zMCNP)&&(aG4NDL==aMCNP))
                        {
                            sumErrorVecGbl[j]=(CalcDiff(outDir+isoNameListG4NDL[j]+".txt", isoFSDataG4NDL[j], enerBoundG4NDL[j], isoFSDataMCNP[k], enerBoundMCNP[k], sampling, relevantData));
                            match=true;
                            break;
                        }
                    }
                }
            }
            else
            {
                for(int k=0; k<fsSizeMCNP; k++)
                {
                    ExtractZA(isoNameListMCNP[k], zMCNP, aMCNP);
                    if((zG4NDL==zMCNP)&&(aG4NDL==aMCNP))
                    {
                        sumErrorVecGbl[j]=(CalcDiff(outDir+isoNameListG4NDL[j]+".txt", isoFSDataG4NDL[j], enerBoundG4NDL[j], isoFSDataMCNP[k], enerBoundMCNP[k], sampling, relevantData));
                        match=true;
                        break;
                    }
                }
            }
            while(nBar<int(j*80/fsSizeG4NDL))
            {
                nBar++;
                cout << "#";
                cout.flush();
            }
        }
        cout << "#]" << endl;
    }

    for(int i=0; i<fsSizeG4NDL; i++)
    {
        for(int j=i; j<fsSizeG4NDL; j++)
        {
            if(sumErrorVecGbl[i]<sumErrorVecGbl[j])
            {
                tempDiff=sumErrorVecGbl[i];
                temp=isoNameListG4NDL[i];
                sumErrorVecGbl[i]=sumErrorVecGbl[j];
                isoNameListG4NDL[i]=isoNameListG4NDL[j];
                sumErrorVecGbl[j]=tempDiff;
                isoNameListG4NDL[j]=temp;
            }
        }
    }

    for(int i=0; i<fsSizeG4NDL; i++)
    {
        if(sumErrorVecGbl[i]>-1)
        {
            stream << isoNameListG4NDL[i] << " has summed difference squared of " << sumErrorVecGbl[i] << std::endl;
            sumDiff+=sumErrorVecGbl[i];
        }
    }
    if(sumErrorVecGbl)
        delete [] sumErrorVecGbl;
    SetDataStream( outDir+"FileDiffList.txt", stream, true);

    return sumDiff;
}

#endif

double CalcDiff(string outFileName, VFDNeutronHPFinalState *isoFSDataG4NDL, std::pair <double,double> enerBoundG4NDL, VFDNeutronHPFinalState *isoFSDataMCNP, std::pair <double,double> enerBoundMCNP, bool sampling, bool *relevantData)
{
    stringstream stream;
    G4HadProjectile *projectile;
    int points=1;
    double totalDiff=0.;
    FSSpectrumData fsDataG4NDL, fsDataMCNP;
    G4HadFinalState *resultG4NDLVec, *resultMCNPVec;

    double emax, emin, enStep;
    if(enerBoundG4NDL.first<enerBoundMCNP.first)
    {
        emin=enerBoundMCNP.first;
    }
    else
    {
        emin=enerBoundG4NDL.first;
    }

    if(enerBoundG4NDL.second>enerBoundMCNP.second)
    {
        emax=enerBoundMCNP.second;
    }
    else
    {
        emax=enerBoundG4NDL.second;
    }

    int numLoops=(std::ceil(std::log10(emax/emin)));
    incEnerSize = numPointPerOrder*numLoops+1;
    incEnerVec = new double[incEnerSize];
    incEnerVec[0]=emin;

    for(int i=0; i<numLoops; i++)
    {
        (incEnerVec[i*numPointPerOrder]*10<emax)? (enStep=9*incEnerVec[i*numPointPerOrder]/numPointPerOrder) : (enStep=(emax-incEnerVec[i*numPointPerOrder])/numPointPerOrder);
        for(int step=0; step<numPointPerOrder; step++)
        {
            incEnerVec[i*numPointPerOrder+step+1] = incEnerVec[i*numPointPerOrder+step] + enStep;
        }
    }

    stream << "Comparison of G4NDL Data to the Converted MCNP data using\n" << incEnerSize << " incoming energy points\n" << numBins << " bins per histogram\n" << std::endl;
    for(int i=0; i<numDataTypes; i++)
    {
        stream << relevantData[i] << '\t';
    }
    stream << '\n' << endl;

    for(int i=0; i<incEnerSize; i++)
    {
        stream << std::setw(14) << std::right << incEnerVec[i];
        if(((i+1)%6==0)||(i==incEnerSize-1))
            stream << '\n';
    }
    stream << '\n' << endl;

    SetDataStream( outFileName, stream, true);


    for(int i=0; i<incEnerSize; i++)
    {
        nParticle->SetKineticEnergy(incEnerVec[i]);
        projectile = new G4HadProjectile(*nParticle);

        fsDataG4NDL.Clear();
        fsDataMCNP.Clear();

        if(sampling)
        {
            points=numSamples;
        }
        for(int l=0; l<points; l++)
        {
            resultG4NDLVec = isoFSDataG4NDL->ApplyYourself(*projectile);
            resultMCNPVec = isoFSDataMCNP->ApplyYourself(*projectile);
            if((resultG4NDLVec==0)||(resultMCNPVec==0))
            {
                stream.clear();
                stream.str("");
                stream << "No final state data to compare!" << std::endl;
                SetDataStream( outFileName, stream, false);
                delete projectile;
                delete [] incEnerVec;
                return 0.;
            }

            fsDataG4NDL.AddData(resultG4NDLVec);
            fsDataMCNP.AddData(resultMCNPVec);
        }
        totalDiff+=fsDataG4NDL.CompareFSData(outFileName, fsDataMCNP, relevantData);

        delete projectile;
    }
    delete [] incEnerVec;
    return totalDiff/points;
}



//ExtractZA
//extracts the Z and the A isotope numbers from the file name
void ExtractZA(string fileName, int &Z, int &A)
{
    std::size_t startPos=0;
    stringstream ss;
    ElementNames* elementNames;
    IsotopeMass* isoMass;
    while(startPos!=fileName.length() && (fileName[startPos]<'0' || fileName[startPos]>'9'))
        startPos++;

    if(startPos==fileName.length())
    {
        //cout << "### File Name Does Not Contian a Z or an A Value " << fileName << " is Invalid for Broadening ###" << std::endl;
        Z=A=-1;
    }
    else
    {
    ////
        std::size_t found1 = fileName.find_first_of('_', startPos);
        if (found1==std::string::npos)
        {
            /*cout << "### File Name Does Not Contian a '_', two are needed, one to seperate the Z and A value, \
            and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << std::endl;*/
            Z=A=-1;
        }
        else
        {
            std::size_t found2 = fileName.find_first_of('_', found1+1);
            if (found2==std::string::npos)
            {
                /*cout << "### File Name Does Not Contian a second '_', two are needed, one to seperate the Z and A value, \
                and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << std::endl;*/
                Z=A=-1;
            }
            else
            {

                ss.str(fileName.substr(startPos, found1));
                ss >> Z;
                ss.str("");
                ss.clear();
                if(((found2-found1-1) > 2) && (fileName[found2-2] == 'm'))
                    ss.str(fileName.substr(found1+1, found2-found1-3));
                else
                    ss.str(fileName.substr(found1+1, found2-found1-1));
                ss >> A;
                ss.str("");
                ss.clear();
                ss.str(fileName.substr(found2+1));
                if (!(elementNames->CheckName(ss.str(), Z)))
                {
                    //cout << "### " << fileName << " does not include the correct element name at the end ###" << std::endl;
                    Z=A=-1;
                }
                ss.str("");
                ss.clear();
                if((Z>0)&&(A==0))
                {
                    isoMass->GetMostNatAbunIso(Z, A);
                }
            }

        }

    }
}

//ExtractZA
//extracts the Z and the A isotope numbers from the file name
void ExtractZA(string fileName, string &Z, string &A)
{
    std::size_t startPos=0;
    stringstream ss;
    while(startPos!=fileName.length() && (fileName[startPos]<'0' || fileName[startPos]>'9'))
        startPos++;

    if(startPos==fileName.length())
    {
        //cout << "### File Name Does Not Contian a Z or an A Value " << fileName << " is Invalid for Broadening ###" << std::endl;
        Z=A="";
    }
    else
    {
    ////
        std::size_t found1 = fileName.find_first_of('_', startPos);
        if (found1==std::string::npos)
        {
            /*cout << "### File Name Does Not Contian a '_', two are needed, one to seperate the Z and A value, \
            and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << std::endl;*/
            Z=A="";
        }
        else
        {
            std::size_t found2 = fileName.find_first_of('_', found1+1);
            if (found2==std::string::npos)
            {
                /*cout << "### File Name Does Not Contian a second '_', two are needed, one to seperate the Z and A value, \
                and one to seperate the A and the Element name " << fileName << " is Invalid for Broadening ###" << std::endl;*/
                Z=A="";
            }
            else
            {

                ss.str(fileName.substr(startPos, found1));
                ss >> Z;
                ss.str("");
                ss.clear();
                ss.str(fileName.substr(found1+1, found2-found1-1));
                ss >> A;
                ss.str("");
                ss.clear();
            }

        }

    }
}

bool DirectoryExists( const char* pzPath )
{
    if ( pzPath == NULL) return false;

    DIR *pDir;
    bool bExists = false;

    pDir = opendir (pzPath);

    if (pDir != NULL)
    {
        bExists = true;
        closedir (pDir);
    }

    return bExists;
}

void GetDataStream( string filename, std::stringstream& ss)
{
   string* data=NULL;
   std::ifstream* in=NULL;
   //string compfilename(filename);

// Use regular text file
  std::ifstream thefData( filename.c_str() , std::ios::in | std::ios::ate );
  if ( thefData.good() )
  {
     int file_size = thefData.tellg();
     thefData.seekg( 0 , std::ios::beg );
     char* filedata = new char[ file_size ];
     while ( thefData )
     {
        thefData.read( filedata , file_size );
     }
     thefData.close();
     data = new string ( filedata , file_size );
     delete [] filedata;
  }
  else
  {
// found no data file
//                 set error bit to the stream
     ss.setstate( std::ios::badbit );
     if(master)
        cout << std::endl << "### failed to open ascii file " << filename << " ###" << std::endl;
  }

   if (data != NULL)
   {
        ss.str(*data);
        if(data[0][data->size()-1]!='\n')
            ss << "\n";
        ss.seekg( 0 , std::ios::beg );
    }

   if(in!=NULL)
   {
        in->close();
        delete in;
   }

    if(data!=NULL)
        delete data;
}


void SetDataStream( string filename , std::stringstream& ss, bool overWrite )
{
    // Use regular text file
    string compfilename(filename);

    if(compfilename.substr((compfilename.length()-2),2)==".z")
    {
        compfilename.erase((compfilename.length()-2),2);
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
                if(master)
                    cout << "\n #### Error the size of the stringstream is invalid ###" << std::endl;
                break;
            }
         }
         out->write(filedata, file_size);
        if (out->fail())
        {
            if(master)
                cout << std::endl << "writing the ascii data to the output file " << compfilename << " failed" << std::endl
                 << " may not have permission to delete an older version of the file" << std::endl;
        }

         delete [] filedata;
    }
    else
    {
    // found no data file
    //                 set error bit to the stream
     ss.setstate( std::ios::badbit );

        if(master)
            cout << std::endl << "### failed to write to ascii file " << compfilename << " ###" << std::endl;
    }
    out->close();
    delete out;
   ss.str("");
}
