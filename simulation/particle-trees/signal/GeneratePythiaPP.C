#include "TennGen.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <string>

#include "Pythia8/Pythia.h"

using namespace Pythia8;
using namespace tenngen;
using namespace std;



const Int_t MAXPARTS = 2000;
const Int_t nPtBins = 26;
const Float_t pthardbin_rhic[nPtBins]= {10.0, 12.5, 15.0, 17.5, 20.0, 22.5, 25.0, 27.5, 30.0, 33.0, 37.0, 40.0, 42.5, 47.0, 50.0, 52.5, 57.5, 
    60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, -1};

const Float_t pthardbin_lhc[nPtBins]= {10.0 , 20.0 , 30.0 , 40.0 , 50.0, 
    60.0, 70.0, 80.0,   90.0,   100.0,  110.0,  120.0,  130.0,  140.0,  150.0,  160.0,  170.0,  180.0,  190.0,  200.0,  220.0,  240.0,  260.0,  280.0, 300.0,   -1};


int GeneratePythiaPP(int collen, int nevent, int ptbin){

    int nparts;
    Float_t particle_px[MAXPARTS], particle_py[MAXPARTS], particle_pz[MAXPARTS], particle_e[MAXPARTS];
    Float_t weight, ptmin, ptmax, xsec_over_eventweight;

    string datadir = "root-files/";
    if (mkdir(datadir.c_str(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    
    TString filename = Form("%s%dGeV_PP_ptbin%d",datadir.c_str(),collen,ptbin);
    filename= filename + ".root";
    TFile *outFile = new TFile(filename.Data(),"RECREATE");
    TTree* outTree = new TTree("tree", "tree");
    TTree* eventInfo = new TTree("eventInfo","eventInfo");

    eventInfo->Branch("nevent",&nevent,"nevent/I");
    eventInfo->Branch("ptmin",&ptmin,"ptmin/F");
    eventInfo->Branch("ptmax",&ptmax,"ptmax/F");
    eventInfo->Branch("ptbinID", &ptbin, "ptbinID/I");
    eventInfo->Branch("xsec_over_eventweight",&xsec_over_eventweight,"xsec_over_eventweight/F");

    outTree->Branch("weight", &weight, "weight/F");
    outTree->Branch("nparts",&nparts,"nparts/I");
    outTree->Branch("particle_px", particle_px, "particle_px[nparts]/F");
    outTree->Branch("particle_py", particle_py, "particle_py[nparts]/F");
    outTree->Branch("particle_pz", particle_pz, "particle_pz[nparts]/F");
    outTree->Branch("particle_e", particle_e, "particle_e[nparts]/F");

    Pythia pythia;
    Pythia8::Settings& settings = pythia.settings;
    const Pythia8::Info& info = pythia.info;
    Event& event = pythia.event;
    if(collen == 200 )pythia.readFile("pythiaSettings_RHIC.cmnd");
    else if(collen == 2760) pythia.readFile("pythiaSettings_LHC.cmnd");
    else{
        cout << "Invalid collision energy" << endl;
        return 1;
    }
    settings.parm("Main:numberOfEvents", nevent);

    Float_t pthardbin[nPtBins];
    if(collen == 200 ) {
        for (int iptbin = 0; iptbin<nPtBins; iptbin++) pthardbin[iptbin] = pthardbin_rhic[iptbin];
    }
    else if(collen == 2760) {
        for (int iptbin = 0; iptbin<nPtBins; iptbin++) pthardbin[iptbin] = pthardbin_lhc[iptbin];
    }
    
    settings.parm("PhaseSpace:pTHatMin", pthardbin[ptbin]);
    settings.parm("PhaseSpace:pTHatMax", pthardbin[ptbin+1]);
    pythia.init();

    for (int j=0; j<nevent; j++){
        if (!pythia.next()) continue;
        weight = info.weight();
        nparts = 0;
        for (int i = 0; i < event.size(); i++){
            if (event[i].isFinal() && event[i].isCharged()){
               particle_px[nparts] = event[i].px();
                particle_py[nparts] = event[i].py();
                particle_pz[nparts] = event[i].pz();
                particle_e[nparts] = event[i].e();
                nparts++;
            }
        }

        outTree->Fill();
        for (int k =0; k<nparts; k++){
            particle_px[k] = 0.0;
            particle_py[k] = 0.0;
            particle_pz[k] = 0.0;
            particle_e[k] = 0.0;
        }



    }

    
    ptmin = pthardbin[ptbin];
    ptmax = pthardbin[ptbin+1];
    xsec_over_eventweight = (info.sigmaGen() / info.weightSum());
    eventInfo->Fill();
    outFile->Write();
    outFile->Close();
    delete outFile;   
    return 0;
}

int main(int argc, char** argv){
    
    int nevent = 100;
    int ptbin = 0;
    if(argc != 4){
        cout << "Usage: ./GeneratePythiaPP <collen> <nevents> <ptbin>" << endl;
        return 1;
    }
    int collen = atoi(argv[1]);
    nevent = atoi(argv[2]);
    ptbin = atoi(argv[3]);
    cout << "Generating " << nevent << " events in pt bin " << ptbin << endl;
    cout << "Collision energy: " << collen << endl;
    return GeneratePythiaPP(collen,nevent, ptbin);
}
