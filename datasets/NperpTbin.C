#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TH1D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"

#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace ROOT;

Float_t GetAverageJetPt(TString inputfile){

    ROOT::EnableImplicitMT();
    TString treeName = "tree";
    RDataFrame df(treeName.Data(), inputfile.Data());

    // filter events
    auto df_filter = df.Filter("jet_pt_raw > 10 && jet_pt_raw < 15");
    auto df_avg = df_filter.Mean("jet_pt_raw");
    Float_t avg_jet_pt = df_avg.GetValue();

    ROOT::DisableImplicitMT();
    return avg_jet_pt;
}

Int_t PlotNperpTbin(TString inputfile,TString outputfile){

    // Get average jet pt
    Float_t avg_jet_pt = GetAverageJetPt(inputfile);
    cout << "Average jet pt: " << avg_jet_pt << endl;
    // open input file
    TFile fin(inputfile,"READ");
    if(!fin.IsOpen()){
        cout << "Error: input file " << inputfile << " not found." << endl;
        return -1;
    }
    TString partTree = "tree";
    TTreeReader finTree(partTree, &fin);

    TTreeReaderValue<Int_t> event_cent_bin(finTree, "event_cent_bin");
    TTreeReaderValue<Int_t> event_pt_hard_bin(finTree, "event_pt_hard_bin");
    TTreeReaderValue<Float_t> event_weight(finTree, "event_weight");
    TTreeReaderValue<Float_t> jet_pt_raw(finTree, "jet_pt_raw");
    TTreeReaderValue<Int_t> jet_nparts(finTree, "jet_nparts");
    TTreeReaderValue<Int_t> jet_nparts_pythia(finTree, "jet_nparts_pythia");
    TTreeReaderArray<Float_t> jet_track_pt(finTree, "jet_track_pt");
    TTreeReaderArray<Float_t> jet_pythia_tracks_pt(finTree, "jet_pythia_tracks_pt");

    // Float_t track_pt_bins[17] =  {1, 1.3318, 1.7758,2.3677,3.1569,4.2092,5.6123,7.4831,9.9775,13.3033,17.7377,23.6503,31.5337,42.0449,56.0599,74.7466,99.66};
    // {0.207612457, 0.408304498, 0.598615917, 0.799307958, 0.996539792, 1.200692042, 1.397923875, 1.619377163, 1.816608997, 2.027681661, 2.200692042, 2.415224913, 2.629757785, 2.823529412, 3.006920415, 3.221453287, 3.429065744};
    //  {1, 1.3318, 1.7758,2.3677,3.1569,4.2092,5.6123,7.4831,9.9775,13.3033,17.7377,23.6503,31.5337,42.0449,56.0599,74.7466,99.66};

    // Float_t track_pt_bins[7] = {0.090717953, 0.135335283, 0.201896518, 0.301194212, 0.449328964, 0.670320046, 1.0};
    Float_t track_pt_bins[7] = {0.067205513, 0.100258844, 0.149568619, 0.22313016, 0.332871084, 0.496585304, 0.740818221};
    Int_t n_track_pt_bins = 6;

    for (Int_t i = 0; i < n_track_pt_bins+1; i++){
        track_pt_bins[i] = track_pt_bins[i]*avg_jet_pt;
    }
    
    // Create Histograms
    TH1D* h_nparts_vs_pt_track = new TH1D("h_nparts_vs_pt_track", "h_nparts_vs_pt_track", n_track_pt_bins, track_pt_bins);
    TH1D* h_nparts_vs_pt_track_pythia = new TH1D("h_nparts_vs_pt_track_pythia", "h_nparts_vs_pt_track_pythia", n_track_pt_bins, track_pt_bins);
    TH1D* h_jetpT = new TH1D("h_jetpT", "h_jetpT", 100, 0, 100);

    // jet counter
    Int_t n_jets = 0;
    Float_t pt_high = 15;
    Float_t pt_low = 10;

    Int_t njets_below_mean = 0;
    Int_t njets_above_mean = 0;
    // loop over events
    while(finTree.Next()){
        // loop over jets
        if(*jet_pt_raw < pt_low || *jet_pt_raw > pt_high) continue;

        h_jetpT->Fill(*jet_pt_raw);

        if(*jet_pt_raw < avg_jet_pt) njets_below_mean++;
        
        if(*jet_pt_raw > avg_jet_pt) njets_above_mean++;

        Int_t npythiatracks = *jet_nparts_pythia;
        Int_t ntracks = *jet_nparts;
        for(Int_t itrack = 0; itrack < ntracks; itrack++){
            Float_t track_pt = jet_track_pt[itrack];
            h_nparts_vs_pt_track->Fill(track_pt);
        }
        for (Int_t itrack = 0; itrack < npythiatracks; itrack++){
            Float_t track_pt = jet_pythia_tracks_pt[itrack];
            h_nparts_vs_pt_track_pythia->Fill(track_pt);
        }
        n_jets++;

    }

    cout << "Total number of jets: " << n_jets << endl;
    cout << "Number of jets below mean: " << njets_below_mean << endl;
    cout << "Number of jets above mean: " << njets_above_mean << endl;
    cout << "Percentage of jets below mean: " << (Float_t)njets_below_mean/(Float_t)n_jets << endl;
    cout << "Percentage of jets above mean: " << (Float_t)njets_above_mean/(Float_t)n_jets << endl;

    // normalize histograms by number of jets
    for(Int_t ix = 1; ix<h_nparts_vs_pt_track->GetNbinsX()+1; ix++){
        // Float_t bin_width = h_nparts_vs_pt_track->GetBinWidth(ix);
        Float_t bin_content = h_nparts_vs_pt_track->GetBinContent(ix);
        Float_t bin_error = h_nparts_vs_pt_track->GetBinError(ix);
        h_nparts_vs_pt_track->SetBinContent(ix, bin_content/n_jets);
        h_nparts_vs_pt_track->SetBinError(ix, bin_error/n_jets);
    }
    for(Int_t ix = 1; ix<h_nparts_vs_pt_track_pythia->GetNbinsX()+1; ix++){
        // Float_t bin_width = h_nparts_vs_pt_track_pythia->GetBinWidth(ix);
        Float_t bin_content = h_nparts_vs_pt_track_pythia->GetBinContent(ix);
        Float_t bin_error = h_nparts_vs_pt_track_pythia->GetBinError(ix);
        h_nparts_vs_pt_track_pythia->SetBinContent(ix, bin_content/n_jets);
        h_nparts_vs_pt_track_pythia->SetBinError(ix, bin_error/n_jets);
    }

    Float_t total_nparts = 0;
    Float_t total_nparts_pythia = 0;
    for(Int_t ix = 1; ix<h_nparts_vs_pt_track->GetNbinsX()+1; ix++){
        Float_t bin_content = h_nparts_vs_pt_track->GetBinContent(ix);
        Float_t bin_content_pythia = h_nparts_vs_pt_track_pythia->GetBinContent(ix);
        total_nparts += bin_content;
        total_nparts_pythia += bin_content_pythia;
    }
    cout << "Total number of particles: " << total_nparts << endl;
    cout << "Total number of particles (pythia): " << total_nparts_pythia << endl;

    // print percentage of particles in first bin
    cout << "Percentage of particles in first bin: " << 100*(h_nparts_vs_pt_track->GetBinContent(1)/total_nparts) << endl;
    cout << "Percentage of particles in first bin (pythia): " << 100*(h_nparts_vs_pt_track_pythia->GetBinContent(1)/total_nparts_pythia) << endl;
    cout << "Number of particles in first bin: " << h_nparts_vs_pt_track->GetBinContent(1) << endl;
    cout << "Number of particles in first bin (pythia): " << h_nparts_vs_pt_track_pythia->GetBinContent(1) << endl;
 

    // h_nparts_vs_pt_track->Scale(1.0/n_jets);
    // h_nparts_vs_pt_track_pythia->Scale(1.0/n_jets);

    // create output file
    TFile *fout = new TFile(outputfile,"RECREATE");
    fout->cd();
    h_nparts_vs_pt_track->Write();
    h_nparts_vs_pt_track_pythia->Write();
    h_jetpT->Write();
    fout->Close();

    // close input file
    fin.Close();

    return 0;
}

Int_t main(Int_t argc, Char_t** argv){


    if(argc < 3){
        cout << "Usage: ./NperpTbin <input file> <output file>" << endl;
        return -1;
    }
    TString input_file = argv[1];
    TString output_file = argv[2];
    return PlotNperpTbin(input_file, output_file);
}