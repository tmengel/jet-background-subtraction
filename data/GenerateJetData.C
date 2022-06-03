#include "TennGen.h"
#include "Pythia8/Pythia.h"
#include "TMath.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"
#include "fastjet/tools/Subtractor.hh" 
#include <fastjet/config.h>


using namespace tenngen;
using namespace Pythia8;
using namespace fastjet;
using namespace std;
using namespace ROOT;


const Int_t MAXPARTS = 6000;
Int_t nParts, nJetParts;
Float_t constituentPt[MAXPARTS] , constituentEta[MAXPARTS] , constituentPhi[MAXPARTS] , constituentEnergy[MAXPARTS], constituentTrueIndex[MAXPARTS];
Float_t KtPt[MAXPARTS] , KtNparts[MAXPARTS] , KtArea[MAXPARTS], KtEta[MAXPARTS], KtPhi[MAXPARTS];
Float_t Area, rhoArea, rhoNpart, JetEta, JetPhi, JetPt, JetTruePt, EventAvgPt, EventMedianPt, rhoEst;
Float_t weight;
Int_t Events, nKtJets, ptbinID;
Float_t ptmin,ptmax, xsec_over_eventweight;

const Int_t nPtBins = 20;
const Float_t pthardbin[nPtBins]= {5.,7.,8.,10.,11.,12.,15.,16.,17.,20.,21.,23.,25.,27.,30.,35.,40.,45.,50.,-1 };

int GenerateJetData(const Int_t nevent, const Int_t ptbin, double R, double etaRange, int centbin){

    double ghost_maxEta = 2*R + etaRange;
    double selectorRap = etaRange;
    double ktR= R + 0.2;
    double ptjetmin = 5.0;
    ptbinID = ptbin;

    ROOT::EnableImplicitMT();
    string datadir = "jet-trees-raw";
    if (mkdir(datadir.c_str(), 0777) == -1) std::cerr << "Output " << strerror(errno) << endl;

    TString centstring, jetparamstring;
    if(centbin == 0) centstring = "0to10";
    else if(centbin == 1) centstring = "20to40";
    if(R == 0.2) jetparamstring = "R02";
    else if(R == 0.4) jetparamstring = "R04";

    TString filename = datadir+"/AuAu_200GeV_"+jetparamstring+"_"+centstring+"_ptbin"+Form("%d",ptbin)+".root";
    TFile *outFile = new TFile(filename.Data(),"RECREATE");

    //string TFOUT = datadir+"/200GeV_R04_0to10cent_ptbin"+std::to_string(ptbin)+".root";
    //TFile* outFile = new TFile(TFOUT.c_str(),"RECREATE"); 
    TTree* outTree = new TTree("JetTree", "JetTree");

    //TTree* ktTree = new TTree("KTJetTree", "KTJetTree");
    
    TTree* eventInfo = new TTree("eventInfo","eventInfo");
    eventInfo->Branch("nevent",&Events,"nevent/I");
    eventInfo->Branch("ptmin",&ptmin,"ptmin/F");
    eventInfo->Branch("ptmax",&ptmax,"ptmax/F");
    eventInfo->Branch("ptbinID", &ptbinID, "ptbinID/I");
    eventInfo->Branch("xsec_over_eventweight",&xsec_over_eventweight,"xsec_over_eventweight/F");
   
    // ktTree->Branch("nKtJets", &nKtJets, "nKtJets/I");
    // ktTree->Branch("KtPt", KtPt, "KtPt[nKtJets]/F");
    // ktTree->Branch("KtNparts", KtNparts, "KtNparts[nKtJets]/F");
    // ktTree->Branch("KtArea", KtArea, "KtArea[nKtJets]/F");
    // ktTree->Branch("KtEta", KtEta, "KtEta[nKtJets]/F");
    // ktTree->Branch("KtPhi", KtPhi, "Kt[nKtJets]/F");

    outTree->Branch("weight", &weight, "weight/F");
    outTree->Branch("ptbinID", &ptbinID, "ptbinID/I");
    outTree->Branch("rhoArea", &rhoArea, "rhoArea/F");
    outTree->Branch("rhoNpart", &rhoNpart, "rhoNpart/F" );
    outTree->Branch("rhoEst", &rhoEst, "rhoEst/F" );
    outTree->Branch("EventAvgPt", &EventAvgPt, "EventAvgPt/F" );
    outTree->Branch("EventMedianPt", &EventMedianPt, "EventMedianPt/F" );
    outTree->Branch("JetEta", &JetEta, "JetEta/F");
    outTree->Branch("JetPhi", &JetPhi, "JetPhi/F");
    outTree->Branch("JetPt", &JetPt, "JetPt/F");
    outTree->Branch("JetTruePt", &JetTruePt, "JetTruePt/F");
    outTree->Branch("Area", &Area, "Area/F");
    outTree->Branch("nJetParts", &nJetParts, "nJetParts/I");
    outTree->Branch("constituentPt", constituentPt, "constituentPt[nJetParts]/F");
    outTree->Branch("constituentEta", constituentEta, "constituentEta[nJetParts]/F");
    outTree->Branch("constituentPhi", constituentPhi, "constituentPhi[nJetParts]/F");
    outTree->Branch("constituentTrueIndex", constituentTrueIndex, "constituentTrueIndex[nJetParts]/F");
    
    Pythia pythia;
    Pythia8::Settings& settings = pythia.settings;
    const Pythia8::Info& info = pythia.info;
    Event& event = pythia.event;
    pythia.readFile("pythiaSettings.cmnd");
    settings.parm("Main:numberOfEvents", nevent);
    settings.parm("PhaseSpace:pTHatMin", pthardbin[ptbin]);
    settings.parm("PhaseSpace:pTHatMax", pthardbin[ptbin+1]);

    pythia.init();

    TennGen tg;
    tg.defaultSettings200(true);
    if(centbin ==0) tg.setcent(0);
    else if(centbin ==1) tg.setcent(2);
    tg.seteta(1.1);
    tg.set_Stream(true);
    TGEvent tmpEvent;
    tg.init();

    vector<PseudoJet> particles;
    vector<PseudoJet> jets;
    vector<PseudoJet> constituents;
    vector<Float_t> calcVec;


 

    for (int ievent = 0; ievent < nevent; ++ievent){

        nParts = 0;
        Float_t Event_pt =0;
        tmpEvent.clear();
        particles.clear();
        jets.clear();
        constituents.clear();

        calcVec.clear();

        particles.resize(0);
        calcVec.resize(0);

        
        if (!pythia.next()) continue;

        tmpEvent = tg.next();

        for (int i = 0; i < event.size(); i++){
            if (event[i].isFinal() && event[i].isCharged()){
                Event_pt+=event[i].pT();
                particles.push_back(PseudoJet(event[i].px(),event[i].py(),event[i].pz(),event[i].e()));
                particles[nParts].set_user_index(1);
                nParts++;
            }
        }

        for (int j = 0; j < tmpEvent.size(); j++){  
                Event_pt+=tmpEvent[j].Pt();
                particles.push_back(PseudoJet(tmpEvent[j].Px(),tmpEvent[j].Py(),tmpEvent[j].Pz(),tmpEvent[j].E()));
                particles[nParts].set_user_index(0);
                nParts++;
        }

        EventAvgPt = Event_pt/nParts;
   
        fastjet::GhostedAreaSpec area_spec(ghost_maxEta);
        fastjet::AreaDefinition area_def(fastjet::active_area, area_spec);
        fastjet::AreaDefinition area_def_bkgd(active_area_explicit_ghosts, GhostedAreaSpec(ghost_maxEta));
        fastjet::Selector selector = SelectorAbsRapMax(ghost_maxEta) * (!SelectorNHardest(2));
        fastjet::JetDefinition jet_def(antikt_algorithm, R);
        fastjet::JetDefinition jet_def_bkgd(kt_algorithm, 0.4);
        fastjet::JetMedianBackgroundEstimator bge(selector, jet_def_bkgd, area_def_bkgd);
        bge.set_particles(particles);
        fastjet::BackgroundEstimate bkgd_estimate = bge.estimate();

        fastjet::ClusterSequenceArea cs(particles, jet_def, area_def);
        fastjet::ClusterSequenceArea cskt(particles, jet_def_bkgd, area_def);

        jets = sorted_by_pt(cskt.inclusive_jets());
        nKtJets = jets.size();
        for(int i =0; i<nKtJets; i++){
            constituents.resize(0);
            constituents.clear();
            constituents =  sorted_by_pt(jets[i].constituents());
            calcVec.push_back(jets[i].pt()/constituents.size());
            // KtNparts[i] = constituents.size();
            // KtArea[i] = jets[i].area();
            // KtPt[i] = jets[i].pt();
            // KtEta[i] = jets[i].eta();
            // KtPhi[i]= jets[i].phi();
        }
      //  ktTree->Fill();
        

        sort(calcVec.begin(), calcVec.end());
        rhoNpart = calcVec[int(calcVec.size()/2)];

        calcVec.resize(0);
        calcVec.clear();
        constituents.resize(0);
        constituents.clear();

        for(int i =0; i<nKtJets; i++){
            calcVec.push_back(jets[i].pt()/jets[i].area());
            // KtNparts[i] = 0;
            // KtArea[i]=0;
            // KtPt[i]=0;
            // KtEta[i] =0;
            // KtPhi[i] =0;
        }
        sort(calcVec.begin(), calcVec.end());
        rhoArea = calcVec[int(calcVec.size()/2)];

        jets.resize(0);
        jets.clear();



        jets = sorted_by_pt(cs.inclusive_jets());

        constituents = sorted_by_pt(particles);
        EventMedianPt = constituents[int(constituents.size()/2)].pt();
        constituents.resize(0);
        constituents.clear();

        for(int i = 0; i < jets.size(); i++){
            if (abs(jets[i].eta())<=etaRange){

                constituents.resize(0);
                constituents.clear();
                constituents =  sorted_by_pt(jets[i].constituents());
                
                Area = jets[i].area();
                rhoEst = bkgd_estimate.rho();
                JetEta = jets[i].eta();
                JetPhi = jets[i].phi();
                JetPt = jets[i].pt();
                nJetParts = constituents.size();
                JetTruePt = 0;
                for(int j = 0; j< nJetParts; j++){
                 // if(constituents[j].user_index()==-1){
                    constituentPt[j] = constituents[j].pt();
                    constituentEta[j] = constituents[j].eta();
                    constituentPhi[j] = constituents[j].phi();
                    constituentTrueIndex[j] = constituents[j].user_index();
                    if(constituents[j].user_index() == 1) JetTruePt += constituents[j].pt();
                }
            
           
               if(JetTruePt ==0) continue;

                weight = info.weight();
                outTree->Fill(); 
                for(int j = 0; j< nJetParts; j++){
                    constituentPt[j] = 0;
                    constituentEta[j] = 0; 
                    constituentPhi[j] = 0;
                    constituentTrueIndex[j] = 0;
                  
                }   
                   
            }    

        }      
        
    }
    
    Events = nevent;
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


     if(argc != 3){
        cout << "Usage: ./GenerateJetData <jetparam> <centrality>" << endl;
        return 1;
    }
    
    int nevent = 10000;
    double etaRange = 1.1;
    double R = atof(argv[1]);
    int centbin = atoi(argv[2]);
    cout << "R = " << R << " centrality = " << centbin << endl;
    

    for(int i =0; i < (nPtBins-1); i++){
    //for(int i =0; i<1; i++){
        GenerateJetData(nevent,i,R,etaRange,centbin);
        gSystem->ProcessEvents();
    } 
    return 0;
}
