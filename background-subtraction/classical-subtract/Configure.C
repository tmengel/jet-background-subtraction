#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>

#include <iostream>
#include <vector>
#include <string>


using namespace std;
using namespace ROOT;


// main function takes in
Int_t main(Int_t argc, Char_t** argv){

    if(argc < 3){
        std::cout << "Usage: " << argv[0] << " <matched.root> <outputfile.root>" << std::endl;
        return 1;
    }

    TString inputfile = argv[1];
    TString outputfile = argv[2];

    // open input file
    TFile *f = new TFile(inputfile, "READ");
    TTree *t = (TTree*)f->Get("tree");
    
    Float_t jet_pt_raw, jet_pt_truth, median_pt_over_npart;
    Int_t jet_nparts_pythia;

    // set branch addresses
    t->SetBranchAddress("jet_pt_raw", &jet_pt_raw);
    t->SetBranchAddress("jet_pt_truth", &jet_pt_truth);
    t->SetBranchAddress("jet_nparts_pythia", &jet_nparts_pythia);
    
    Int_t min_nparts_pythia = 100;
    Int_t max_nparts_pythia = 0;
    Float_t min_jet_pt_reco = 100;
    Float_t max_jet_pt_reco = 0;

    // get min and max values for nparts and jet pt
    Int_t nentries = t->GetEntries();
    for ( Int_t ientry = 0; ientry < nentries; ientry++){
        t->GetEntry(ientry);
        if (jet_nparts_pythia < min_nparts_pythia) min_nparts_pythia = jet_nparts_pythia;
        if (jet_nparts_pythia > max_nparts_pythia) max_nparts_pythia = jet_nparts_pythia;
        if (jet_pt_raw < min_jet_pt_reco) min_jet_pt_reco = jet_pt_raw;
        if (jet_pt_raw > max_jet_pt_reco) max_jet_pt_reco = jet_pt_raw;
    }

    // round min and max values to nearest 10
    min_nparts_pythia = int(min_nparts_pythia/10)*10;
    max_nparts_pythia = int(max_nparts_pythia/10)*10 + 10;
    min_jet_pt_reco = int(min_jet_pt_reco/10)*10;
    max_jet_pt_reco = int(max_jet_pt_reco/10)*10.0 + 10.0;

    // create jet pt bins
    Int_t n_jet_pt_bins = int((max_jet_pt_reco - min_jet_pt_reco)/10);
    Float_t jet_pt_bins[n_jet_pt_bins+1];
    for (Int_t i = 0; i < n_jet_pt_bins+1; i++){
        jet_pt_bins[i] = min_jet_pt_reco + i*10;
    }
    std::vector<Float_t> jet_pt_bins_vec(jet_pt_bins, jet_pt_bins + sizeof(jet_pt_bins)/sizeof(Float_t));
    // get integer npart bins
    Int_t n_npart_bins = int(max_nparts_pythia - min_nparts_pythia);


    // create output file
    TFile *fout = new TFile(outputfile, "RECREATE");
    // create TH2Ds for jet npart vs. jet pt
    TH2D *nparts_pythia_vs_jet_pt_truth = new TH2D("nparts_pythia_vs_jet_pt_truth", "nparts_pythia_vs_jet_pt_truth", n_jet_pt_bins, min_jet_pt_reco, max_jet_pt_reco, n_npart_bins, min_nparts_pythia, max_nparts_pythia);
    TH2D *nparts_pythia_vs_jet_pt_raw = new TH2D("nparts_pythia_vs_jet_pt_raw", "nparts_pythia_vs_jet_pt_raw", n_jet_pt_bins, min_jet_pt_reco, max_jet_pt_reco, n_npart_bins, min_nparts_pythia, max_nparts_pythia);

    // fill TH2Ds
    for ( Int_t ientry = 0; ientry < nentries; ientry++){
        t->GetEntry(ientry);
        nparts_pythia_vs_jet_pt_truth->Fill(jet_pt_truth, jet_nparts_pythia);
        nparts_pythia_vs_jet_pt_raw->Fill(jet_pt_raw, jet_nparts_pythia);
    }

    // create tprofiles
    TProfile *nparts_pythia_vs_jet_pt_truth_profile = nparts_pythia_vs_jet_pt_truth->ProfileX();
    TProfile *nparts_pythia_vs_jet_pt_raw_profile = nparts_pythia_vs_jet_pt_raw->ProfileX();

    // convert tprofiles to th1ds
    TH1D *nparts_pythia_vs_jet_pt_truth_hist = new TH1D("nparts_pythia_vs_jet_pt_truth_hist", "nparts_pythia_vs_jet_pt_truth_hist", n_jet_pt_bins, min_jet_pt_reco, max_jet_pt_reco);
    TH1D *nparts_pythia_vs_jet_pt_raw_hist = new TH1D("nparts_pythia_vs_jet_pt_raw_hist", "nparts_pythia_vs_jet_pt_raw_hist", n_jet_pt_bins, min_jet_pt_reco, max_jet_pt_reco);

    for ( Int_t ibin = 0; ibin <  nparts_pythia_vs_jet_pt_truth_profile->GetNbinsX(); ibin++){
        nparts_pythia_vs_jet_pt_truth_hist->SetBinContent(ibin+1, nparts_pythia_vs_jet_pt_truth_profile->GetBinContent(ibin+1));
        nparts_pythia_vs_jet_pt_truth_hist->SetBinError(ibin+1, nparts_pythia_vs_jet_pt_truth_profile->GetBinError(ibin+1));
    }
    for ( Int_t ibin = 0; ibin <  nparts_pythia_vs_jet_pt_raw_profile->GetNbinsX(); ibin++){
        nparts_pythia_vs_jet_pt_raw_hist->SetBinContent(ibin+1, nparts_pythia_vs_jet_pt_raw_profile->GetBinContent(ibin+1));
        nparts_pythia_vs_jet_pt_raw_hist->SetBinError(ibin+1, nparts_pythia_vs_jet_pt_raw_profile->GetBinError(ibin+1));
    }

    // save histograms
    fout->cd();
    nparts_pythia_vs_jet_pt_truth->Write();
    nparts_pythia_vs_jet_pt_raw->Write();
    nparts_pythia_vs_jet_pt_truth_profile->Write();
    nparts_pythia_vs_jet_pt_raw_profile->Write();
    nparts_pythia_vs_jet_pt_truth_hist->Write();
    nparts_pythia_vs_jet_pt_raw_hist->Write();

    fout->Close();

    // close input file
    f->Close();

    return 0;
   
}


