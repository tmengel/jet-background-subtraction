#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include <vector>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

using namespace ROOT; // RDataFrame lives in here
using namespace std;

int FormatRawJetTrees(TString filename, TString outdir, double R, int centbin)
{

    TString cent, jetparam;
    if(centbin == 0) cent = "0to10";
    else if(centbin == 1) cent = "20to40";
    if(R == 0.2) jetparam = "R02";
    else if(R == 0.4) jetparam = "R04";
    TString outRootFile = "AuAu_200GeV_"+jetparam+"_"+cent+"_ptbin";
    TString outRootFilePath = outdir+outRootFile;

    string jetTree = "JetTree";
    string ktTree = "KTJetTree";
    string eventTree = "eventInfo";


    TFile f(filename.Data());

    TTreeReader jettr(jetTree.c_str(), &f);
    TTreeReader weighttr(eventTree.c_str(), &f);

    Int_t ptbin;
    Float_t eventxsec;
    std::vector<Float_t> pToverNpart;

    //Jet Tree Varibles//
    TTreeReaderValue<Float_t> jet_weight(jettr, "weight");
    TTreeReaderValue<Int_t> jet_pTbin(jettr, "ptbinID");
    TTreeReaderValue<Float_t> median_pT_over_area(jettr, "rhoArea");
    TTreeReaderValue<Float_t> median_Npart_over_area(jettr, "rhoNpart");
    TTreeReaderValue<Float_t> fastjet_rho_estimate(jettr, "rhoEst");
    TTreeReaderValue<Float_t> event_average_pT(jettr, "EventAvgPt");
    TTreeReaderValue<Float_t> event_median_pT(jettr, "EventMedianPt");
    TTreeReaderValue<Float_t> jet_pT(jettr, "JetPt");
    TTreeReaderValue<Float_t> jet_eta(jettr, "JetEta");
    TTreeReaderValue<Float_t> jet_phi(jettr, "JetPhi");
    TTreeReaderValue<Float_t> jet_pT_from_pythia(jettr, "JetTruePt");
    TTreeReaderValue<Float_t> jet_area(jettr, "Area");
    TTreeReaderValue<Int_t> jet_Nparts(jettr, "nJetParts");
    TTreeReaderArray<Float_t> constituent_pT(jettr, "constituentPt");
    TTreeReaderArray<Float_t> constituent_eta(jettr, "constituentEta");
    TTreeReaderArray<Float_t> constituent_phi(jettr, "constituentPhi");
    TTreeReaderArray<Float_t> constituent_truth_index(jettr, "constituentTrueIndex");

    //Event Tree Varibles//
    TTreeReaderValue<Float_t> xsec_over_totalweight(weighttr, "xsec_over_eventweight");
    TTreeReaderValue<Int_t> event_pTbin(weighttr, "ptbinID");
    
    weighttr.Next(); 
    ptbin = *event_pTbin; 
    eventxsec = *xsec_over_totalweight;
    TString outfile = outRootFilePath+Form("%d",ptbin)+".root";
    TFile *fout = new TFile(outfile.Data(), "RECREATE");
    TTree *outTree = new TTree("outTree", "outTree");


    //string outfile = outputDir+"/AuAu_200GeV_R04_0to10cent_ptbin"+to_string(ptbin)+".root";
    // TFile* fout = new TFile(outfile.c_str(),"RECREATE");
    // TTree* outTree= new TTree("outTree","outTree");

    Float_t weight, scaledweight, median_pt_over_area, median_npart_over_area, fastjet_rho, event_average_pt, event_median_pt;
    Float_t jetpt, jeteta, jetphi, jetpt_pythia, jetarea;
    Int_t jetnparts, PtHardBin;
    Float_t track0pt, track1pt, track2pt, track3pt, track4pt, average_track_pt, median_track_pt, jetangularity, trackptvariance,trackptskewness,trackptkurtosis;
    Float_t jetpt_pythia_fraction, NumberBased, AreaBased;

    outTree->Branch("weight", &weight, "weight/F");
    outTree->Branch("scaledweight", &scaledweight, "scaledweight/F");
    outTree->Branch("median_pt_over_area", &median_pt_over_area, "median_pt_over_area/F");
    outTree->Branch("median_npart_over_area", &median_npart_over_area, "median_npart_over_area/F" );
    outTree->Branch("fastjet_rho", &fastjet_rho, "fastjet_rho/F" );
    outTree->Branch("event_average_pt", &event_average_pt, "event_average_pt/F" );
    outTree->Branch("event_median_pt", &event_median_pt, "event_median_pt/F" );
    outTree->Branch("jetpt", &jetpt, "jetpt/F");
    outTree->Branch("jeteta", &jeteta, "jeteta/F");
    outTree->Branch("jetphi", &jetphi, "jetphi/F");
    outTree->Branch("jetpt_pythia", &jetpt_pythia, "jetpt_pythia/F");
    outTree->Branch("AreaBased", &AreaBased, "AreaBased/F");
    outTree->Branch("NumberBased", &NumberBased, "NumberBased/F");
    outTree->Branch("jetpt_pythia_fraction", &jetpt_pythia_fraction, "jetpt_pythia_fraction/F");
    outTree->Branch("jetarea", &jetarea, "jetarea/F");
    outTree->Branch("jetnparts", &jetnparts, "jetnparts/I" );
    outTree->Branch("PtHardBin", &PtHardBin, "PtHardBin/I" );
    outTree->Branch("track0pt", &track0pt, "track0pt/F" );
    outTree->Branch("track1pt", &track1pt, "track1pt/F" );
    outTree->Branch("track2pt", &track2pt, "track2pt/F" );
    outTree->Branch("track3pt", &track3pt, "track3pt/F" );
    outTree->Branch("track4pt", &track4pt, "track4pt/F" );
    outTree->Branch("average_track_pt", &average_track_pt, "average_track_pt/F");
    outTree->Branch("median_track_pt", &median_track_pt, "median_track_pt/F");
    outTree->Branch("jetangularity", &jetangularity, "jetangularity/F");
    outTree->Branch("trackptvariance", &trackptvariance, "trackptvariance/F");
    outTree->Branch("trackptskewness", &trackptskewness, "trackptskewness/F");
    outTree->Branch("trackptkurtosis", &trackptkurtosis, "trackptkurtosis/F");

    while (jettr.Next()) {

        weight = *jet_weight;
        scaledweight = weight*eventxsec;
        median_pt_over_area = *median_pT_over_area;
        median_npart_over_area = *median_Npart_over_area;
        fastjet_rho = *fastjet_rho_estimate;
        event_average_pt = *event_average_pT;
        event_median_pt = *event_median_pT;
        jetpt = *jet_pT;
        jeteta = *jet_eta;
        jetphi = *jet_phi;
        jetpt_pythia = *jet_pT_from_pythia;
        jetarea = *jet_area;
        jetnparts = *jet_Nparts;
        AreaBased = jetpt - jetarea*median_pt_over_area;
        NumberBased = jetpt - median_npart_over_area*jetnparts;
        PtHardBin = ptbin;
       
        track0pt = constituent_pT[0];
        track1pt = constituent_pT[1];
        track2pt = constituent_pT[2];
        track3pt = constituent_pT[3];
        track4pt = constituent_pT[4];

        average_track_pt = jetpt/jetnparts;

        pToverNpart.clear();
        pToverNpart.resize(0);

        Float_t temp, pttemp;
        jetangularity =0;
        trackptvariance=0;
        trackptskewness=0;
        trackptkurtosis=0;
        pttemp =0;

        for (int i=0;i < constituent_pT.GetSize(); i++) {
            pToverNpart.push_back(constituent_pT[i]); 
            temp = (constituent_pT[i] - average_track_pt);
            jetangularity+= constituent_pT[i]*TMath::Sqrt( (jeteta-constituent_eta[i])*(jeteta-constituent_eta[i])+ (jetphi-constituent_phi[i])*(jetphi-constituent_phi[i]) );
            trackptvariance+= TMath::Power(temp,2.0);
            trackptskewness+= TMath::Power(temp,3.0);
            trackptkurtosis+= TMath::Power(temp,4.0);
            if(constituent_truth_index[i] == 1.0) pttemp +=  constituent_pT[i];
        }

        sort(pToverNpart.begin(), pToverNpart.end());
        median_track_pt = pToverNpart[int(pToverNpart.size()/2)];
        jetpt_pythia_fraction = pttemp/jetpt;
        jetangularity = jetangularity/jetpt;
        trackptvariance = trackptvariance/jetnparts;
        trackptskewness = trackptskewness/jetnparts;
        trackptkurtosis = trackptkurtosis/jetnparts;

        if(TMath::Abs(jeteta) < 1.1  && jetpt > 5.0 &&jetpt_pythia_fraction !=0 )outTree->Fill();
    }
    
    fout->Write();
    fout->Close();
    return 0;
}

int main( int argc, char** argv){

    if(argc != 3){
        cout << "Usage: ./FormatRawJetTrees <jetparam> <centrality>" << endl;
        return 1;
    }

    double R = atof(argv[1]);
    int centbin = atoi(argv[2]);
    cout << "R = " << R << " centrality = " << centbin << endl;

    int numptbins = 19;
    TString cent, jetparam;

    if(centbin == 0) cent = "0to10";
    else if(centbin == 1) cent = "20to40";

    if(R == 0.2) jetparam = "R02";
    else if(R == 0.4) jetparam = "R04";

    TString fileheader = "AuAu_200GeV_"+jetparam+"_"+cent+"_ptbin";
    TString datadir = "jet-trees-raw/";
    TString outdir = "../src/pre-processed-data/jet-trees-formatted/";

    for(int i =0; i<numptbins; i++){
        TString filename = datadir+fileheader+Form("%d",i)+".root";   
        cout << "Formatting: " << filename << endl;
        FormatRawJetTrees(filename, outdir, R , centbin);
    }

    return 0;
}
