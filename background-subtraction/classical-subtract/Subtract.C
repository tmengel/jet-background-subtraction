#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TSystem.h>

#include <iostream>
#include <vector>
#include <string>


using namespace std;
using namespace ROOT;


// main function takes in
Int_t main(Int_t argc, Char_t** argv){

    if(argc < 3){
        std::cout << "Usage: " << argv[0] << " <configure_file.root> <subtract_file.root>" << std::endl;
        return 1;
    }

    TString configurefile = argv[1];
    TString infile = argv[2];

    // open configure file
    TFile *f = new TFile(configurefile, "READ");
    TH1D* h_nparts_vs_jet_pt_reco = (TH1D*)f->Get("nparts_pythia_vs_jet_pt_raw_hist");

    // get jet pt bins
    Int_t n_jet_pt_bins = h_nparts_vs_jet_pt_reco->GetNbinsX();
    Float_t jet_pt_bins[n_jet_pt_bins+1];
    for (Int_t i = 0; i < n_jet_pt_bins+1; i++){
        jet_pt_bins[i] = h_nparts_vs_jet_pt_reco->GetXaxis()->GetBinLowEdge(i+1);
    }
    // get upper edge of last bin
    jet_pt_bins[n_jet_pt_bins] = h_nparts_vs_jet_pt_reco->GetXaxis()->GetBinUpEdge(n_jet_pt_bins);

    // get content of bins
    Float_t nparts_pythia[n_jet_pt_bins];
    for (Int_t i = 0; i < n_jet_pt_bins; i++){
        nparts_pythia[i] = h_nparts_vs_jet_pt_reco->GetBinContent(i+1);
    }

    // close configure file
    f->Close();
    
    // open subtract file
    TFile *fin = new TFile(infile, "READ");
    TTree *t = (TTree*)fin->Get("tree");

    // get list of branchs and their types
    std::vector<TString> int_branch_names;
    std::vector<TString> float_branch_names;
    Int_t nInts = 0;
    Int_t nFloats = 0;

    // nessicary branches
    std::vector<TString> used_variables = {"jet_nparts", "jet_pt_raw", "median_pt_over_npart", "jet_area", "median_pt_over_area"};
    TObjArray *branches = t->GetListOfBranches();
    for (Int_t i = 0; i < branches->GetEntries(); i++){
        TBranch *branch = (TBranch*)branches->At(i);
        TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);

        TString branch_name = branch->GetName();
        TString leaf_type = leaf->GetTypeName();

        Int_t save = 1;
        for (Int_t j = 0; j < used_variables.size(); j++){
            if (branch_name == used_variables[j]){
                save = 0;
            }
        }

        if(save){
            if (leaf_type == "Int_t"){
                int_branch_names.push_back(branch_name);
                nInts++;
            }
            else if (leaf_type == "Float_t"){
                float_branch_names.push_back(branch_name);
                nFloats++;
            }
        }
    }


    // set branch addresses
    Int_t jet_nparts;
    Float_t jet_pt_raw, median_pt_over_npart, jet_area, median_pt_over_area;
    Float_t jet_pt_multiplicity_corrected, jet_pt_area_corrected;
    Int_t passive_ints[nInts];
    Float_t passive_floats[nFloats];

    t->SetBranchAddress("jet_nparts", &jet_nparts);
    t->SetBranchAddress("jet_pt_raw", &jet_pt_raw);
    t->SetBranchAddress("median_pt_over_npart", &median_pt_over_npart);
    t->SetBranchAddress("jet_area", &jet_area);
    t->SetBranchAddress("median_pt_over_area", &median_pt_over_area);
    for (Int_t i = 0; i < nInts; i++){
        t->SetBranchAddress(int_branch_names[i], &passive_ints[i]);
    }
    for (Int_t i = 0; i < nFloats; i++){
        t->SetBranchAddress(float_branch_names[i], &passive_floats[i]);
    }

    // create tmp file
    TString tmp_file_name = infile;
    tmp_file_name.ReplaceAll(".root", "_tmp.root");
    TFile *tmp_file = new TFile(tmp_file_name, "RECREATE");
    TTree *tmp_tree = new TTree("tree", "tree");
    tmp_tree->Branch("jet_nparts", &jet_nparts, "jet_nparts/I");
    tmp_tree->Branch("jet_pt_raw", &jet_pt_raw, "jet_pt_raw/F");
    tmp_tree->Branch("median_pt_over_npart", &median_pt_over_npart, "median_pt_over_npart/F");
    tmp_tree->Branch("jet_pt_multiplicity_corrected", &jet_pt_multiplicity_corrected, "jet_pt_multiplicity_corrected/F");
    tmp_tree->Branch("jet_area", &jet_area, "jet_area/F");
    tmp_tree->Branch("median_pt_over_area", &median_pt_over_area, "median_pt_over_area/F");
    tmp_tree->Branch("jet_pt_area_corrected", &jet_pt_area_corrected, "jet_pt_area_corrected/F");

    for (Int_t i = 0; i < nInts; i++){
        tmp_tree->Branch(int_branch_names[i], &passive_ints[i], Form("%s/I", int_branch_names[i].Data()));
    }
    for (Int_t i = 0; i < nFloats; i++){
        tmp_tree->Branch(float_branch_names[i], &passive_floats[i], Form("%s/F", float_branch_names[i].Data()));
    }
       
    // loop over entries
    Int_t nentries = t->GetEntries();
    for (Int_t i = 0; i < nentries; i++){
        t->GetEntry(i);

        // get jet pt bin
        Int_t pt_bin = 0;
        for( Int_t bin = 0; bin < n_jet_pt_bins; bin++){
            if (jet_pt_raw >= jet_pt_bins[bin] && jet_pt_raw < jet_pt_bins[bin+1]){
                pt_bin = bin;
                break;
            }
        }
       
       // get nparts from pythia
        Int_t nparts = nparts_pythia[pt_bin];

        // calculate corrected jet pt
        jet_pt_multiplicity_corrected = jet_pt_raw - median_pt_over_npart*(jet_nparts - nparts);
        jet_pt_area_corrected = jet_pt_raw - median_pt_over_area*jet_area;
        // fill tree
        tmp_tree->Fill();
    }

    // write tmp tree to tmp file
    tmp_file->cd();
    tmp_tree->Write();
    tmp_file->Close();
    // close file
    fin->Close();

    // cp tmp file to infile
    TString cmd = Form("cp %s %s", tmp_file_name.Data(), infile.Data());
    gSystem->Exec(cmd);

    // delete tmp file
    cmd = Form("rm %s", tmp_file_name.Data());
    gSystem->Exec(cmd);
    return 0;
   
}


