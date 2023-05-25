#include <iostream>
#include <string>

#include "TennGen.h"
#include "Pythia8/Pythia.h"

using namespace tenngen;
using namespace Pythia8;

const int MAXTRACKS = 2000;
float TrackPt[MAXTRACKS];
float TrackEta[MAXTRACKS];
float TrackPhi[MAXTRACKS];
int TrackKF[MAXTRACKS];
int TrackIndex[MAXTRACKS];

void clearBuff();
int main(){


    float pthardmin = 5.0;
    float pthardmax = 60.0;
    int centBin = 0;
    int nEvent = 10;
    float etaRange = 1.1;
    int eventMultiplicity;

    TGEventList mergedEvents;
    TGEvent tmpEvent;
    TGParticle tmpPart;

    string dir = "output-pythia";
    if (mkdir(dir.c_str(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    string TfileOut = dir+"/MergedEvent.root";

    TFile* outFile = new TFile(TfileOut.c_str(),"RECREATE");
    TTree* outTree = new TTree("particleTree", "");
    TTree* eventTree = new TTree("EventInfo", "");

    eventTree->Branch("numEvents",&nEvent,"numEvents/I");
    eventTree->Branch("cent",&centBin,"cent/I");
    eventTree->Branch("etaRange",&etaRange,"etaRange/F");
    eventTree->Branch("ptHardmin",&pthardmin,"ptHardmin/F");
    eventTree->Branch("ptHardmax",&pthardmax,"ptHardmax/F");

    eventTree->Fill();

    outTree->Branch("multiplicity", &eventMultiplicity, "multiplicity/I");
    outTree->Branch("particlePt", TrackPt, "particlePt[multiplicity]/F");
    outTree->Branch("particleEta", TrackEta, "particleEta[multiplicity]/F");
    outTree->Branch("particlePhi", TrackPhi, "particlePhi[multiplicity]/F");
    outTree->Branch("particleKF", TrackKF, "particleKF[multiplicity]/I");
    outTree->Branch("particleTruth", TrackIndex, "particleKF[multiplicity]/I");
    
    Pythia pythia;
    Settings& settings = pythia.settings;
    Event& event = pythia.event;

    pythia.readFile("pythiaSettings.cmnd");

    settings.parm("Main:numberOfEvents",nEvent);
    settings.parm("PhaseSpace:pTHatMin", pthardmin);
    settings.parm("PhaseSpace:pTHatMax", pthardmax);

    pythia.init();

    int nAbort = 5;
    int iAbort = 0;

    TennGen tg;
    tg.defaultSettings200(true);
    tg.setcent(centBin);
    tg.seteta(etaRange);
    tg.set_Stream(true);

    tg.init();
  


    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      tmpEvent = tg.next();
      tmpEvent.setEventUserIndex(0);
      tmpEvent.Eta(etaRange);

      if (!pythia.next()) { /*Check to make sure pythia didn't crash or output junk*/
            if (++iAbort < nAbort) continue;
            cout << " Event generation aborted prematurely, owing to error!\n";
            break;
      }
      
      for (int i = 0; i < event.size(); ++i){
        
        if (event[i].isFinal() && event[i].isCharged() && abs(event[i].eta())<=etaRange){
            
            tmpPart=TGParticle(event[i].pT(),event[i].eta(), event[i].phi()+M_PI, int(event[i].id()));
            tmpPart.setUserIndex(1);
            tmpEvent+=tmpPart;
            tmpPart.clear();
        }


      }

      mergedEvents+=tmpEvent;
      tmpEvent.clear();

    }

    for(int i =0; i< nEvent; i++){
        eventMultiplicity=mergedEvents[i].size();
        for(int j =0; j<eventMultiplicity; j++){
                TrackPt[j] = mergedEvents[i][j].Pt();
                TrackEta[j]= mergedEvents[i][j].Eta();
                TrackPhi[j]= mergedEvents[i][j].Phi();
                TrackKF[j] = mergedEvents[i][j].Kf();
                TrackIndex[j] = mergedEvents[i][j].UserIndex();
        }
        outTree->Fill();
        clearBuff();
    }
    outFile->Write();
    outFile->Close();   
    return 0;

}
void clearBuff(){
    for(int j =0; j<MAXTRACKS; j++){
        TrackPt[j] = 0;
        TrackEta[j]= 0;
        TrackPhi[j]= 0;
        TrackKF[j] = 0;
        TrackIndex[j] = -1;
    }
}