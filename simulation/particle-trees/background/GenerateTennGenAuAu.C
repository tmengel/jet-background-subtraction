#include "TennGen.h"

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"

#include <iostream>
#include <string>

using namespace tenngen;
using namespace std;
const Int_t MAXPARTS = 2000;
int GenerateTennGenAuAu(int nevent, int centbin, float etaRange, int index){

    int nparts;
    Float_t particle_px[MAXPARTS], particle_py[MAXPARTS], particle_pz[MAXPARTS], particle_e[MAXPARTS];

    // Tenngen
    TennGen tg;
  
    tg.setcollen(200);
    tg.setnevent(1000);
    tg.seteta(etaRange);
    tg.setpsiN(1,-1.0); // -1.0 means it will be randomly choosen between 0-2 pi
    tg.setpsiN(2,0.0);
    tg.setpsiN(3,-1.0);
    tg.setpsiN(4,0.0);

    // Output file
    TString datadir = "root-files/AuAu/";
    if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    // if(index == 0){
    //     datadir = datadir + "no-flow/";
    //     if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    //     tg.setvN(1,false);
    //     tg.setvN(2,false);
    //     tg.setvN(3,false);
    //     tg.setvN(4,false);
    // }
    // else if(index==1){
    //     datadir = datadir + "flow/"; 
    //     if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    //     tg.setvN(1,true);
    //     tg.setvN(2,true);
    //     tg.setvN(3,true);
    //     tg.setvN(4,true);
    // }
    if(centbin ==0){
        tg.setcent(0);
        datadir = datadir + "0-10/";
        if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    } 
    else if(centbin ==1){
        tg.setcent(1);
        datadir = datadir + "10-20/";
        if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    } 
    else if(centbin==2){
        tg.setcent(2);
        datadir = datadir + "20-40/";
        if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    }
    else if(centbin==3){
        tg.setcent(3);
        datadir = datadir + "40-60/";
        if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    } 
    else{
        cout << "centrality bin not defined" << endl;
        return 0;
    }

    // create director if needed
    if (mkdir(datadir.Data(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;
    
    tg.do_Histos(false); // will make debug root plots
    tg.do_TTree(false); // will make TTree 
    tg.set_Batch(false);
    tg.set_Stream(true);


    TString filename = datadir+"200GeV_AuAu.root";
    TFile *outFile = new TFile(filename.Data(),"RECREATE");

    TTree* outTree = new TTree("tree", "tree");
    TTree* eventInfo = new TTree("eventInfo","eventInfo");

    eventInfo->Branch("nevent",&nevent,"nevent/I");
    eventInfo->Branch("centbin",&centbin,"centbin/I");
    eventInfo->Branch("etarange",&etaRange,"etarange/F");

    outTree->Branch("nparts",&nparts,"nparts/I");
    outTree->Branch("particle_px", particle_px, "particle_px[nparts]/F");
    outTree->Branch("particle_py", particle_py, "particle_py[nparts]/F");
    outTree->Branch("particle_pz", particle_pz, "particle_pz[nparts]/F");
    outTree->Branch("particle_e", particle_e, "particle_e[nparts]/F");

    tg.init();
    TGEvent tmpEvent;

    for (int i=0; i<nevent; i++){
            tmpEvent.clear();
            tmpEvent = tg.next();
            nparts = 0;
            for (int j=0; j< tmpEvent.size(); j++){
                particle_px[j] = tmpEvent[j].Px();
                particle_py[j] = tmpEvent[j].Py();
                particle_pz[j] = tmpEvent[j].Pz();
                particle_e[j] = tmpEvent[j].E();
                nparts++;
            }
            outTree->Fill();
            for (int k =0; k<nparts; k++){
                particle_px[k] = 0.0;
                particle_py[k] = 0.0;
                particle_pz[k] = 0.0;
                particle_e[k] = 0.0;
            }

        }
    
    
    eventInfo->Fill();
    outFile->Write();
    outFile->Close();
    delete outFile;   
    return 0;
}

int main(int argc, char** argv){
    
    int nevent = 100;
    int centbin = 0;
    float etaRange = 1.1;

    if(argc <= 4){
        cout << "Usage: ./GenerateTennGenAuAu <nevents> <centrality> <etarange> <vnmode>" << endl;
        return 1;
    }
    nevent = atoi(argv[1]);
    centbin = atoi(argv[2]);
    etaRange = atof(argv[3]);
    int index = -1; 
    index = atoi(argv[4]);  
    cout << "Generating " << nevent << " events with " << centbin << " cent bin and " << etaRange << " eta range" << endl;
    return GenerateTennGenAuAu(nevent, centbin, etaRange,index); 
}
