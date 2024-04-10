#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TH1D.h"
#include "TH2D.h"
#include "TH2.h"
#include "TH1.h"

#include "TFile.h"
#include "TTree.h"

#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "THnSparse.h"

#include "TSystem.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldErrors.h"

#include "TString.h"

#include "TMath.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <iostream>
#include <vector>

#endif


using namespace std;
using namespace ROOT;

Float_t GetMax(TString inputfile, TString variable, TString counter_variable){
    TFile fin(inputfile, "READ");
    TTreeReader reader("tree", &fin);
    TTreeReaderValue<Int_t> in_counter(reader, counter_variable);
    TTreeReaderArray<Float_t> in_variable(reader, variable);
    Float_t max = 0.0;
    while(reader.Next()){
        Int_t counter = *in_counter;
        for(Int_t i = 0; i < counter; i++){
            Float_t value = in_variable[i];
            if(value > max) max = value;
        }
    }
    fin.Close();
    return max;

}

Float_t GetMin(TString inputfile, TString variable, TString counter_variable){
    TFile fin(inputfile, "READ");
    TTreeReader reader("tree", &fin);
    TTreeReaderValue<Int_t> in_counter(reader, counter_variable);
    TTreeReaderArray<Float_t> in_variable(reader, variable);
    Float_t min = 1000000.0;
    while(reader.Next()){
        Int_t counter = *in_counter;
        for(Int_t i = 0; i < counter; i++){
            Float_t value = in_variable[i];
            if(value < min) min = value;
        }
    }
    fin.Close();
    return min;

}

TString NormalizeDiffFlow(TString inputfile){
    // cout << "Normalizing diff flow" << endl;
    // open input file
    TFile fin(inputfile, "READ");
    // get tree
    TTreeReader reader("tree", &fin);
    TTreeReaderValue<Int_t> in_event_rfp_M(reader, "event_rfp_M");
    TTreeReaderValue<Int_t> in_n_reco_bins(reader, "n_reco_bins");
    TTreeReaderArray<Float_t> in_reco_correlation_2nd(reader, "reco_correlation_2nd");
    TTreeReaderArray<Float_t> in_reco_correlation_3rd(reader, "reco_correlation_3rd");
    TTreeReaderArray<Float_t> in_reco_correlation_4th(reader, "reco_correlation_4th");
    TTreeReaderArray<Float_t> in_reco_correlation_weight(reader, "reco_correlation_weight");
    TTreeReaderArray<Float_t> in_fake_correlation_2nd(reader, "fake_correlation_2nd");
    TTreeReaderArray<Float_t> in_fake_correlation_3rd(reader, "fake_correlation_3rd");
    TTreeReaderArray<Float_t> in_fake_correlation_4th(reader, "fake_correlation_4th");
    TTreeReaderArray<Float_t> in_fake_correlation_weight(reader, "fake_correlation_weight");
    TTreeReaderArray<Float_t> in_matched_correlation_2nd(reader, "matched_correlation_2nd");
    TTreeReaderArray<Float_t> in_matched_correlation_3rd(reader, "matched_correlation_3rd");
    TTreeReaderArray<Float_t> in_matched_correlation_4th(reader, "matched_correlation_4th");
    TTreeReaderArray<Float_t> in_matched_correlation_weight(reader, "matched_correlation_weight");
    TTreeReaderValue<Int_t> in_n_truth_bins(reader, "n_truth_bins");
    TTreeReaderArray<Float_t> in_truth_correlation_2nd(reader, "truth_correlation_2nd");
    TTreeReaderArray<Float_t> in_truth_correlation_3rd(reader, "truth_correlation_3rd");
    TTreeReaderArray<Float_t> in_truth_correlation_4th(reader, "truth_correlation_4th");
    TTreeReaderArray<Float_t> in_truth_correlation_weight(reader, "truth_correlation_weight");
    TTreeReaderArray<Float_t> in_missed_correlation_2nd(reader, "missed_correlation_2nd");
    TTreeReaderArray<Float_t> in_missed_correlation_3rd(reader, "missed_correlation_3rd");
    TTreeReaderArray<Float_t> in_missed_correlation_4th(reader, "missed_correlation_4th");
    TTreeReaderArray<Float_t> in_missed_correlation_weight(reader, "missed_correlation_weight");

    Int_t n_truth_bins = 9;
    Float_t truth_bins[10] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70, 80, 90, 100};
    Int_t n_reco_bins = 12;
    Float_t reco_bins[13] = {-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50, 60, 70, 80, 90, 100, 110};


    Int_t event_rfp_M;
    // output variables
    Float_t reco_correlation_2nd[12] = {0.0};
    Float_t reco_correlation_3rd[12] = {0.0};
    Float_t reco_correlation_4th[12] = {0.0};
    Float_t truth_correlation_2nd[12] = {0.0};
    Float_t truth_correlation_3rd[12] = {0.0};
    Float_t truth_correlation_4th[12] = {0.0};
    Float_t fake_correlation_2nd[12] = {0.0};
    Float_t fake_correlation_3rd[12] = {0.0};
    Float_t fake_correlation_4th[12] = {0.0};
    Float_t missed_correlation_2nd[12] = {0.0};
    Float_t missed_correlation_3rd[12] = {0.0};
    Float_t missed_correlation_4th[12] = {0.0};
    Float_t matched_correlation_2nd[12] = {0.0};
    Float_t matched_correlation_3rd[12] = {0.0};
    Float_t matched_correlation_4th[12] = {0.0};
    Float_t reco_correlation_weight[12] = {0.0};
    Float_t truth_correlation_weight[12] = {0.0};
    Float_t fake_correlation_weight[12] = {0.0};
    Float_t missed_correlation_weight[12] = {0.0};
    Float_t matched_correlation_weight[12] = {0.0};

    // outfile
    TString output_file = inputfile;
    output_file.ReplaceAll(".root", "_normalized.root");
    TFile *fout = new TFile(output_file, "RECREATE");
    TTree *tree = new TTree("tree", "tree");
    tree->Branch("event_rfp_M", &event_rfp_M, "event_rfp_M/I");
    tree->Branch("n_reco_bins", &n_reco_bins, "n_reco_bins/I");
    tree->Branch("reco_correlation_2nd", reco_correlation_2nd, "reco_correlation_2nd[n_reco_bins]/F");
    tree->Branch("reco_correlation_3rd", reco_correlation_3rd, "reco_correlation_3rd[n_reco_bins]/F");
    tree->Branch("reco_correlation_4th", reco_correlation_4th, "reco_correlation_4th[n_reco_bins]/F");
    tree->Branch("reco_correlation_weight", reco_correlation_weight, "reco_correlation_weight[n_reco_bins]/F");
    tree->Branch("fake_correlation_2nd", fake_correlation_2nd, "fake_correlation_2nd[n_reco_bins]/F");
    tree->Branch("fake_correlation_3rd", fake_correlation_3rd, "fake_correlation_3rd[n_reco_bins]/F");
    tree->Branch("fake_correlation_4th", fake_correlation_4th, "fake_correlation_4th[n_reco_bins]/F");
    tree->Branch("fake_correlation_weight", fake_correlation_weight, "fake_correlation_weight[n_reco_bins]/F");
    tree->Branch("matched_correlation_2nd", matched_correlation_2nd, "matched_correlation_2nd[n_reco_bins]/F");
    tree->Branch("matched_correlation_3rd", matched_correlation_3rd, "matched_correlation_3rd[n_reco_bins]/F");
    tree->Branch("matched_correlation_4th", matched_correlation_4th, "matched_correlation_4th[n_reco_bins]/F");
    tree->Branch("matched_correlation_weight", matched_correlation_weight, "matched_correlation_weight[n_reco_bins]/F");
    tree->Branch("n_truth_bins", &n_truth_bins, "n_truth_bins/I");
    tree->Branch("truth_correlation_2nd", truth_correlation_2nd, "truth_correlation_2nd[n_reco_bins]/F");
    tree->Branch("truth_correlation_3rd", truth_correlation_3rd, "truth_correlation_3rd[n_reco_bins]/F");
    tree->Branch("truth_correlation_4th", truth_correlation_4th, "truth_correlation_4th[n_reco_bins]/F");
    tree->Branch("truth_correlation_weight", truth_correlation_weight, "truth_correlation_weight[n_reco_bins]/F");
    tree->Branch("missed_correlation_2nd", missed_correlation_2nd, "missed_correlation_2nd[n_reco_bins]/F");
    tree->Branch("missed_correlation_3rd", missed_correlation_3rd, "missed_correlation_3rd[n_reco_bins]/F");
    tree->Branch("missed_correlation_4th", missed_correlation_4th, "missed_correlation_4th[n_reco_bins]/F");
    tree->Branch("missed_correlation_weight", missed_correlation_weight, "missed_correlation_weight[n_reco_bins]/F");

    cout << "Normalizing diff flow" << endl;
    while(reader.Next()){
       event_rfp_M = *in_event_rfp_M;
        n_reco_bins = *in_n_reco_bins;
        n_truth_bins = *in_n_truth_bins;

         // zero out arrays
        for(Int_t i = 0; i < n_reco_bins; i++){
            reco_correlation_2nd[i] = 0.0;
            reco_correlation_3rd[i] = 0.0;
            reco_correlation_4th[i] = 0.0;
            reco_correlation_weight[i] = 0.0;
            fake_correlation_2nd[i] = 0.0;
            fake_correlation_3rd[i] = 0.0;
            fake_correlation_4th[i] = 0.0;
            fake_correlation_weight[i] = 0.0;
            matched_correlation_2nd[i] = 0.0;
            matched_correlation_3rd[i] = 0.0;
            matched_correlation_4th[i] = 0.0;
            matched_correlation_weight[i] = 0.0;
        // }
        // for(Int_t i = 0; i < n_truth_bins; i++){
            truth_correlation_2nd[i] = 0.0;
            truth_correlation_3rd[i] = 0.0;
            truth_correlation_4th[i] = 0.0;
            truth_correlation_weight[i] = 0.0;
            missed_correlation_2nd[i] = 0.0;
            missed_correlation_3rd[i] = 0.0;
            missed_correlation_4th[i] = 0.0;
            missed_correlation_weight[i] = 0.0;
        }


        // loop over reco jets
        for (Int_t ijet=0; ijet< n_reco_bins; ijet++){
            Float_t jet_correlation_2nd = in_reco_correlation_2nd[ijet];
            Float_t jet_correlation_3rd = in_reco_correlation_3rd[ijet];
            Float_t jet_correlation_4th = in_reco_correlation_4th[ijet];
            Float_t jet_correlation_weight = in_reco_correlation_weight[ijet];
            if(jet_correlation_weight == 0.0){
                reco_correlation_2nd[ijet] = 0.0;
                reco_correlation_3rd[ijet] = 0.0;
                reco_correlation_4th[ijet] = 0.0;

            }
            else{
                reco_correlation_2nd[ijet] = jet_correlation_2nd/jet_correlation_weight;
                reco_correlation_3rd[ijet] = jet_correlation_3rd/jet_correlation_weight;
                reco_correlation_4th[ijet] = jet_correlation_4th/jet_correlation_weight;
            }
            Float_t total_weight = jet_correlation_weight*event_rfp_M;
            reco_correlation_weight[ijet] = total_weight;
        }

        // loop over fake jets
        for (Int_t ijet=0; ijet< n_reco_bins; ijet++){
            Float_t jet_correlation_2nd = in_fake_correlation_2nd[ijet];
            Float_t jet_correlation_3rd = in_fake_correlation_3rd[ijet];
            Float_t jet_correlation_4th = in_fake_correlation_4th[ijet];
            Float_t jet_correlation_weight = in_fake_correlation_weight[ijet];
            if(jet_correlation_weight == 0.0){
                fake_correlation_2nd[ijet] = 0.0;
                fake_correlation_3rd[ijet] = 0.0;
                fake_correlation_4th[ijet] = 0.0;

            }
            else{
                fake_correlation_2nd[ijet] = jet_correlation_2nd/jet_correlation_weight;
                fake_correlation_3rd[ijet] = jet_correlation_3rd/jet_correlation_weight;
                fake_correlation_4th[ijet] = jet_correlation_4th/jet_correlation_weight;
            }
            Float_t total_weight = jet_correlation_weight*event_rfp_M;
            fake_correlation_weight[ijet] = total_weight;
        }

        // loop over matched jets
        for (Int_t ijet=0; ijet< n_reco_bins; ijet++){
            Float_t jet_correlation_2nd = in_matched_correlation_2nd[ijet];
            Float_t jet_correlation_3rd = in_matched_correlation_3rd[ijet];
            Float_t jet_correlation_4th = in_matched_correlation_4th[ijet];
            Float_t jet_correlation_weight = in_matched_correlation_weight[ijet];
            if(jet_correlation_weight == 0.0){
                matched_correlation_2nd[ijet] = 0.0;
                matched_correlation_3rd[ijet] = 0.0;
                matched_correlation_4th[ijet] = 0.0;

            }
            else{
                matched_correlation_2nd[ijet] = jet_correlation_2nd/jet_correlation_weight;
                matched_correlation_3rd[ijet] = jet_correlation_3rd/jet_correlation_weight;
                matched_correlation_4th[ijet] = jet_correlation_4th/jet_correlation_weight;
            }
            Float_t total_weight = jet_correlation_weight*event_rfp_M;
            matched_correlation_weight[ijet] = total_weight;
        }

        // loop over truth jets
        for (Int_t ijet=0; ijet< n_reco_bins; ijet++){
            Float_t jet_correlation_2nd = in_truth_correlation_2nd[ijet];
            Float_t jet_correlation_3rd = in_truth_correlation_3rd[ijet];
            Float_t jet_correlation_4th = in_truth_correlation_4th[ijet];
            Float_t jet_correlation_weight = in_truth_correlation_weight[ijet];
            if(jet_correlation_weight == 0.0){
                truth_correlation_2nd[ijet] = 0.0;
                truth_correlation_3rd[ijet] = 0.0;
                truth_correlation_4th[ijet] = 0.0;

            }
            else{
                truth_correlation_2nd[ijet] = jet_correlation_2nd/jet_correlation_weight;
                truth_correlation_3rd[ijet] = jet_correlation_3rd/jet_correlation_weight;
                truth_correlation_4th[ijet] = jet_correlation_4th/jet_correlation_weight;
            }
            Float_t total_weight = jet_correlation_weight*event_rfp_M;
            truth_correlation_weight[ijet] = total_weight;
        }

        // loop over missed jets
        for (Int_t ijet=0; ijet< n_reco_bins; ijet++){
            Float_t jet_correlation_2nd = in_missed_correlation_2nd[ijet];
            Float_t jet_correlation_3rd = in_missed_correlation_3rd[ijet];
            Float_t jet_correlation_4th = in_missed_correlation_4th[ijet];
            Float_t jet_correlation_weight = in_missed_correlation_weight[ijet];
           if(jet_correlation_weight == 0.0){
                missed_correlation_2nd[ijet] = 0.0;
                missed_correlation_3rd[ijet] = 0.0;
                missed_correlation_4th[ijet] = 0.0;

            }
            else{
                missed_correlation_2nd[ijet] = jet_correlation_2nd/jet_correlation_weight;
                missed_correlation_3rd[ijet] = jet_correlation_3rd/jet_correlation_weight;
                missed_correlation_4th[ijet] = jet_correlation_4th/jet_correlation_weight;
            }
            Float_t total_weight = jet_correlation_weight*event_rfp_M;
            missed_correlation_weight[ijet] = total_weight;
        }


        tree->Fill();

    }

    cout << "All Done" << endl;
    fout->cd();
    fout->Write();
    fout->Close();

    fin.Close();

    return output_file;

}

Int_t Unfoldvn(TString diff_input, TString ref_input, TString output_file){
    Int_t n_truth_bins = 9;
    Float_t truth_bins[10] = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70, 80, 90, 100};
    std::vector<Double_t> truth_bins_vec = {10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70, 80, 90, 100};
    Int_t n_reco_bins = 12;
    Float_t reco_bins[13] = {-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50, 60, 70, 80, 90, 100, 110};
    std::vector<Double_t> reco_bins_vec = {-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50, 60, 70, 80, 90, 100, 110};

    TString norm_diff_input = NormalizeDiffFlow(diff_input);
    // norm_diff_input.ReplaceAll(".root", "_normalized.root");
    // Get max and min values
    Float_t max_reco_2nd = GetMax(norm_diff_input, "reco_correlation_2nd", "n_reco_bins");
    Float_t max_reco_3rd = GetMax(norm_diff_input, "reco_correlation_3rd", "n_reco_bins");
    Float_t max_reco_4th = GetMax(norm_diff_input, "reco_correlation_4th", "n_reco_bins");

    Float_t max_fake_2nd = GetMax(norm_diff_input, "fake_correlation_2nd", "n_reco_bins");
    Float_t max_fake_3rd = GetMax(norm_diff_input, "fake_correlation_3rd", "n_reco_bins");
    Float_t max_fake_4th = GetMax(norm_diff_input, "fake_correlation_4th", "n_reco_bins");

    Float_t max_matched_2nd = GetMax(norm_diff_input, "matched_correlation_2nd", "n_reco_bins");
    Float_t max_matched_3rd = GetMax(norm_diff_input, "matched_correlation_3rd", "n_reco_bins");
    Float_t max_matched_4th = GetMax(norm_diff_input, "matched_correlation_4th", "n_reco_bins");

    Float_t max_truth_2nd = GetMax(norm_diff_input, "truth_correlation_2nd", "n_truth_bins");
    Float_t max_truth_3rd = GetMax(norm_diff_input, "truth_correlation_3rd", "n_truth_bins");
    Float_t max_truth_4th = GetMax(norm_diff_input, "truth_correlation_4th", "n_truth_bins");

    Float_t max_missed_2nd = GetMax(norm_diff_input, "missed_correlation_2nd", "n_truth_bins");
    Float_t max_missed_3rd = GetMax(norm_diff_input, "missed_correlation_3rd", "n_truth_bins");
    Float_t max_missed_4th = GetMax(norm_diff_input, "missed_correlation_4th", "n_truth_bins");

    Float_t min_reco_2nd = GetMin(norm_diff_input, "reco_correlation_2nd", "n_reco_bins");
    Float_t min_reco_3rd = GetMin(norm_diff_input, "reco_correlation_3rd", "n_reco_bins");
    Float_t min_reco_4th = GetMin(norm_diff_input, "reco_correlation_4th", "n_reco_bins");

    Float_t min_fake_2nd = GetMin(norm_diff_input, "fake_correlation_2nd", "n_reco_bins");
    Float_t min_fake_3rd = GetMin(norm_diff_input, "fake_correlation_3rd", "n_reco_bins");
    Float_t min_fake_4th = GetMin(norm_diff_input, "fake_correlation_4th", "n_reco_bins");

    Float_t min_matched_2nd = GetMin(norm_diff_input, "matched_correlation_2nd", "n_reco_bins");
    Float_t min_matched_3rd = GetMin(norm_diff_input, "matched_correlation_3rd", "n_reco_bins");
    Float_t min_matched_4th = GetMin(norm_diff_input, "matched_correlation_4th", "n_reco_bins");

    Float_t min_truth_2nd = GetMin(norm_diff_input, "truth_correlation_2nd", "n_truth_bins");
    Float_t min_truth_3rd = GetMin(norm_diff_input, "truth_correlation_3rd", "n_truth_bins");
    Float_t min_truth_4th = GetMin(norm_diff_input, "truth_correlation_4th", "n_truth_bins");

    Float_t min_missed_2nd = GetMin(norm_diff_input, "missed_correlation_2nd", "n_truth_bins");
    Float_t min_missed_3rd = GetMin(norm_diff_input, "missed_correlation_3rd", "n_truth_bins");
    Float_t min_missed_4th = GetMin(norm_diff_input, "missed_correlation_4th", "n_truth_bins");

   
    std::vector<Float_t> max_2nd = {max_reco_2nd, max_fake_2nd, max_matched_2nd, max_truth_2nd, max_missed_2nd};
    std::vector<Float_t> max_3rd = {max_reco_3rd, max_fake_3rd, max_matched_3rd, max_truth_3rd, max_missed_3rd};
    std::vector<Float_t> max_4th = {max_reco_4th, max_fake_4th, max_matched_4th, max_truth_4th, max_missed_4th};

    std::vector<Float_t> min_2nd = {min_reco_2nd, min_fake_2nd, min_matched_2nd, min_truth_2nd, min_missed_2nd};
    std::vector<Float_t> min_3rd = {min_reco_3rd, min_fake_3rd, min_matched_3rd, min_truth_3rd, min_missed_3rd};
    std::vector<Float_t> min_4th = {min_reco_4th, min_fake_4th, min_matched_4th, min_truth_4th, min_missed_4th};


    Double_t max_2nd_all = -999.0;
    Double_t max_3rd_all = -999.0;
    Double_t max_4th_all = -999.0;
    Double_t min_2nd_all = 999.0;
    Double_t min_3rd_all = 999.0;
    Double_t min_4th_all = 999.0;
    for(Int_t i = 0; i < 5; i++){
        if(max_2nd[i] > max_2nd_all) max_2nd_all = max_2nd[i];
        if(max_3rd[i] > max_3rd_all) max_3rd_all = max_3rd[i];
        if(max_4th[i] > max_4th_all) max_4th_all = max_4th[i];
        if(min_2nd[i] < min_2nd_all) min_2nd_all = min_2nd[i];
        if(min_3rd[i] < min_3rd_all) min_3rd_all = min_3rd[i];
        if(min_4th[i] < min_4th_all) min_4th_all = min_4th[i];
    }


    TFile fin(norm_diff_input, "READ");
    if(!fin.IsOpen()){ cout << "Could not open input file" << endl; exit(-1); }
    // get tree
    TTreeReader reader("tree", &fin);
    TTreeReaderValue<Int_t> in_event_rfp_M(reader, "event_rfp_M");
    TTreeReaderValue<Int_t> in_n_reco_bins(reader, "n_reco_bins");
    TTreeReaderArray<Float_t> in_reco_correlation_2nd(reader, "reco_correlation_2nd");
    TTreeReaderArray<Float_t> in_reco_correlation_3rd(reader, "reco_correlation_3rd");
    TTreeReaderArray<Float_t> in_reco_correlation_4th(reader, "reco_correlation_4th");
    TTreeReaderArray<Float_t> in_reco_correlation_weight(reader, "reco_correlation_weight");
    TTreeReaderArray<Float_t> in_fake_correlation_2nd(reader, "fake_correlation_2nd");
    TTreeReaderArray<Float_t> in_fake_correlation_3rd(reader, "fake_correlation_3rd");
    TTreeReaderArray<Float_t> in_fake_correlation_4th(reader, "fake_correlation_4th");
    TTreeReaderArray<Float_t> in_fake_correlation_weight(reader, "fake_correlation_weight");
    TTreeReaderArray<Float_t> in_matched_correlation_2nd(reader, "matched_correlation_2nd");
    TTreeReaderArray<Float_t> in_matched_correlation_3rd(reader, "matched_correlation_3rd");
    TTreeReaderArray<Float_t> in_matched_correlation_4th(reader, "matched_correlation_4th");
    TTreeReaderArray<Float_t> in_matched_correlation_weight(reader, "matched_correlation_weight");
    TTreeReaderValue<Int_t> in_n_truth_bins(reader, "n_truth_bins");
    TTreeReaderArray<Float_t> in_truth_correlation_2nd(reader, "truth_correlation_2nd");
    TTreeReaderArray<Float_t> in_truth_correlation_3rd(reader, "truth_correlation_3rd");
    TTreeReaderArray<Float_t> in_truth_correlation_4th(reader, "truth_correlation_4th");
    TTreeReaderArray<Float_t> in_truth_correlation_weight(reader, "truth_correlation_weight");
    TTreeReaderArray<Float_t> in_missed_correlation_2nd(reader, "missed_correlation_2nd");
    TTreeReaderArray<Float_t> in_missed_correlation_3rd(reader, "missed_correlation_3rd");
    TTreeReaderArray<Float_t> in_missed_correlation_4th(reader, "missed_correlation_4th");
    TTreeReaderArray<Float_t> in_missed_correlation_weight(reader, "missed_correlation_weight");

    Int_t event_rfp_M;
    // output variables
    Float_t reco_correlation_2nd[12] = {0.0};
    Float_t reco_correlation_3rd[12] = {0.0};
    Float_t reco_correlation_4th[12] = {0.0};
    Float_t truth_correlation_2nd[12] = {0.0};
    Float_t truth_correlation_3rd[12] = {0.0};
    Float_t truth_correlation_4th[12] = {0.0};
    Float_t fake_correlation_2nd[12] = {0.0};
    Float_t fake_correlation_3rd[12] = {0.0};
    Float_t fake_correlation_4th[12] = {0.0};
    Float_t missed_correlation_2nd[12] = {0.0};
    Float_t missed_correlation_3rd[12] = {0.0};
    Float_t missed_correlation_4th[12] = {0.0};
    Float_t matched_correlation_2nd[12] = {0.0};
    Float_t matched_correlation_3rd[12] = {0.0};
    Float_t matched_correlation_4th[12] = {0.0};
    Float_t reco_correlation_weight[12] = {0.0};
    Float_t truth_correlation_weight[12] = {0.0};
    Float_t fake_correlation_weight[12] = {0.0};
    Float_t missed_correlation_weight[12] = {0.0};
    Float_t matched_correlation_weight[12] = {0.0};


    Int_t bins[4] = {n_reco_bins, 100, n_truth_bins, 100};
    Double_t xmin2[4] = {reco_bins[0], min_2nd_all, reco_bins[0], min_2nd_all};
    Double_t xmax2[4] = {reco_bins[n_reco_bins], max_2nd_all, reco_bins[n_truth_bins], max_2nd_all};
    Double_t xmin3[4] = {reco_bins[0], min_3rd_all, reco_bins[0], min_3rd_all};
    Double_t xmax3[4] = {reco_bins[n_reco_bins], max_3rd_all, reco_bins[n_truth_bins], max_3rd_all};
    Double_t xmin4[4] = {reco_bins[0], min_4th_all, reco_bins[0], min_4th_all};
    Double_t xmax4[4] = {reco_bins[n_reco_bins], max_4th_all, reco_bins[n_truth_bins], max_4th_all};

    THnSparseD * hn_response_2nd = new THnSparseD("hn_response_2nd", "hn_response_2nd", 4, bins, xmin2, xmax2);
    THnSparseD * hn_measured_2nd = new THnSparseD("hn_measured_2nd", "hn_measured_2nd", 4, bins, xmin2, xmax2);
    THnSparseD * hn_truth_2nd = new THnSparseD("hn_truth_2nd", "hn_truth_2nd", 4, bins, xmin2, xmax2);
    THnSparseD * hn_response_3rd = new THnSparseD("hn_response_3rd", "hn_response_3rd", 4, bins, xmin3, xmax3);
    THnSparseD * hn_measured_3rd = new THnSparseD("hn_measured_3rd", "hn_measured_3rd", 4, bins, xmin3, xmax3);
    THnSparseD * hn_truth_3rd = new THnSparseD("hn_truth_3rd", "hn_truth_3rd", 4, bins, xmin3, xmax3);
    THnSparseD * hn_response_4th = new THnSparseD("hn_response_4th", "hn_response_4th", 4, bins, xmin4, xmax4);
    THnSparseD * hn_measured_4th = new THnSparseD("hn_measured_4th", "hn_measured_4th", 4, bins, xmin4, xmax4);
    THnSparseD * hn_truth_4th = new THnSparseD("hn_truth_4th", "hn_truth_4th", 4, bins, xmin4, xmax4);

    hn_response_2nd->Sumw2();
    hn_measured_2nd->Sumw2();
    hn_truth_2nd->Sumw2();
    hn_response_3rd->Sumw2();
    hn_measured_3rd->Sumw2();
    hn_truth_3rd->Sumw2();
    hn_response_4th->Sumw2();
    hn_measured_4th->Sumw2();
    hn_truth_4th->Sumw2();


    // TH2D * h2_measured_2nd = new TH2D("h2_measured_2nd", "h2_measured_2nd", n_reco_bins, &reco_bins_vec[0], 100, min_2nd_all, max_2nd_all);
    // TH2D * h2_truth_2nd = new TH2D("h2_truth_2nd", "h2_truth_2nd", n_reco_bins, &reco_bins_vec[0], 100, min_2nd_all, max_2nd_all);
    TH1D * h_reco_n_poi = new TH1D("h_reco_n_poi", "h_reco_n_poi", n_reco_bins, &reco_bins_vec[0]);
    TH1D * h_fake_n_poi = new TH1D("h_fake_n_poi", "h_fake_n_poi", n_reco_bins, &reco_bins_vec[0]);
    TH1D * h_matched_n_poi = new TH1D("h_matched_n_poi", "h_matched_n_poi", n_reco_bins, &reco_bins_vec[0]);
    TH1D * h_truth_n_poi = new TH1D("h_truth_n_poi", "h_truth_n_poi", n_reco_bins, &reco_bins_vec[0]);
    TH1D * h_missed_n_poi = new TH1D("h_missed_n_poi", "h_missed_n_poi", n_reco_bins, &reco_bins_vec[0]);

    cout << "Filling Sparse" << endl;
    while(reader.Next()){
        event_rfp_M = *in_event_rfp_M;
        n_reco_bins = *in_n_reco_bins;
        //  read arrays
        for(Int_t i = 0; i < n_reco_bins; i++){
            Float_t pt = (reco_bins[i]+reco_bins[i+1])/2.0;
            reco_correlation_2nd[i] = in_reco_correlation_2nd[i];
            reco_correlation_3rd[i] = in_reco_correlation_3rd[i];
            reco_correlation_4th[i] = in_reco_correlation_4th[i];
            reco_correlation_weight[i] = in_reco_correlation_weight[i];
            fake_correlation_2nd[i] = in_fake_correlation_2nd[i];
            fake_correlation_3rd[i] = in_fake_correlation_3rd[i];
            fake_correlation_4th[i] = in_fake_correlation_4th[i];
            fake_correlation_weight[i] = in_fake_correlation_weight[i];
            matched_correlation_2nd[i] = in_matched_correlation_2nd[i];
            matched_correlation_3rd[i] = in_matched_correlation_3rd[i];
            matched_correlation_4th[i] = in_matched_correlation_4th[i];
            matched_correlation_weight[i] = in_matched_correlation_weight[i];
            truth_correlation_2nd[i] = in_truth_correlation_2nd[i];
            truth_correlation_3rd[i] = in_truth_correlation_3rd[i];
            truth_correlation_4th[i] = in_truth_correlation_4th[i];
            truth_correlation_weight[i] = in_truth_correlation_weight[i];
            missed_correlation_2nd[i] = in_missed_correlation_2nd[i];
            missed_correlation_3rd[i] = in_missed_correlation_3rd[i];
            missed_correlation_4th[i] = in_missed_correlation_4th[i];
            missed_correlation_weight[i] = in_missed_correlation_weight[i];

            h_reco_n_poi->Fill(pt, reco_correlation_weight[i]);
            h_fake_n_poi->Fill(pt, fake_correlation_weight[i]);
            h_matched_n_poi->Fill(pt, matched_correlation_weight[i]);
            h_truth_n_poi->Fill(pt, truth_correlation_weight[i]);
            h_missed_n_poi->Fill(pt, missed_correlation_weight[i]);

            Double_t in_array_fill2[4] = {pt, matched_correlation_2nd[i], pt, truth_correlation_2nd[i]};
            Double_t in_array_fake2[4] = {pt, fake_correlation_2nd[i], 0.0, 0.0};
            Double_t in_array_missed2[4] = {0.0, 0.0, pt, missed_correlation_2nd[i]};
            // 3rd
            Double_t in_array_fill3[4] = {pt, matched_correlation_3rd[i], pt, truth_correlation_3rd[i]};
            Double_t in_array_fake3[4] = {pt, fake_correlation_3rd[i], 0.0, 0.0};
            Double_t in_array_missed3[4] = {0.0, 0.0, pt, missed_correlation_3rd[i]};
            // 4th
            Double_t in_array_fill4[4] = {pt, matched_correlation_4th[i], pt, truth_correlation_4th[i]};
            Double_t in_array_fake4[4] = {pt, fake_correlation_4th[i], 0.0, 0.0};
            Double_t in_array_missed4[4] = {0.0, 0.0, pt, missed_correlation_4th[i]};

            hn_response_2nd->Fill(in_array_fill2, matched_correlation_weight[i]);
            hn_response_2nd->Fill(in_array_fake2, fake_correlation_weight[i]);
            hn_response_2nd->Fill(in_array_missed2, missed_correlation_weight[i]);

            hn_response_3rd->Fill(in_array_fill3, matched_correlation_weight[i]);
            hn_response_3rd->Fill(in_array_fake3, fake_correlation_weight[i]);
            hn_response_3rd->Fill(in_array_missed3, missed_correlation_weight[i]);

            hn_response_4th->Fill(in_array_fill4, matched_correlation_weight[i]);
            hn_response_4th->Fill(in_array_fake4, fake_correlation_weight[i]);
            hn_response_4th->Fill(in_array_missed4, missed_correlation_weight[i]);

            
            Double_t in_array_meas2[4] = {pt, reco_correlation_2nd[i], 0, 0};
            hn_measured_2nd->Fill(in_array_meas2, reco_correlation_weight[i]);

            Double_t in_array_truth2[4] = {0, 0, pt, truth_correlation_2nd[i]};
            // Double_t in_array_missed2[4] = {0, 0, pt, missed_correlation_2nd[i]};
            hn_truth_2nd->Fill(in_array_truth2, truth_correlation_weight[i]);
            hn_truth_2nd->Fill(in_array_missed2, missed_correlation_weight[i]);

            Double_t in_array_meas3[4] = {pt, reco_correlation_3rd[i], 0, 0};
            hn_measured_3rd->Fill(in_array_meas3, reco_correlation_weight[i]);

            Double_t in_array_truth3[4] = {0, 0, pt, truth_correlation_3rd[i]};
            // Double_t in_array_missed3[4] = {0, 0, pt, missed_correlation_3rd[i]};
            hn_truth_3rd->Fill(in_array_truth3, truth_correlation_weight[i]);
            hn_truth_3rd->Fill(in_array_missed3, missed_correlation_weight[i]);

            Double_t in_array_meas4[4] = {pt, reco_correlation_4th[i], 0, 0};
            hn_measured_4th->Fill(in_array_meas4, reco_correlation_weight[i]);

            Double_t in_array_truth4[4] = {0, 0, pt, truth_correlation_4th[i]};
            // Double_t in_array_missed4[4] = {0, 0, pt, missed_correlation_4th[i]};
            hn_truth_4th->Fill(in_array_truth4, truth_correlation_weight[i]);
            hn_truth_4th->Fill(in_array_missed4, missed_correlation_weight[i]);

        }
    
    }
   
    cout << "Opening Reference File" << endl;
    TFile * fout = new TFile(output_file, "RECREATE");
    fout->cd();

    cout << "Creating Projections" << endl;
    Int_t nprojections = n_reco_bins;
    TH1 * hMeas[nprojections];
    TH1 * hTruth[nprojections];
    TH2 * hResponse[nprojections];
    TH1D * hUnfolded[nprojections];

    TH1 * hMeas_2nd[nprojections];
    TH1 * hTruth_2nd[nprojections];
    TH2 * hResponse_2nd[nprojections];
    TH1D * hUnfolded_2nd[nprojections];

    TH1 * hMeas_3rd[nprojections];
    TH1 * hTruth_3rd[nprojections];
    TH2 * hResponse_3rd[nprojections];
    TH1D * hUnfolded_3rd[nprojections];

    TH1 * hMeas_4th[nprojections];
    TH1 * hTruth_4th[nprojections];
    TH2 * hResponse_4th[nprojections];
    TH1D * hUnfolded_4th[nprojections];




    for(Int_t i = 1; i < nprojections; i++){
        hn_measured_2nd->GetAxis(0)->SetRange(i, i+1);
        hn_truth_2nd->GetAxis(2)->SetRange(i, i+1);
        hn_response_2nd->GetAxis(0)->SetRange(i, i+1);
        hn_response_2nd->GetAxis(2)->SetRange(i, i+1);

        hn_measured_3rd->GetAxis(0)->SetRange(i, i+1);
        hn_truth_3rd->GetAxis(2)->SetRange(i, i+1);
        hn_response_3rd->GetAxis(0)->SetRange(i, i+1);
        hn_response_3rd->GetAxis(2)->SetRange(i, i+1);

        hn_measured_4th->GetAxis(0)->SetRange(i, i+1);
        hn_truth_4th->GetAxis(2)->SetRange(i, i+1);
        hn_response_4th->GetAxis(0)->SetRange(i, i+1);
        hn_response_4th->GetAxis(2)->SetRange(i, i+1);



        TH1 * mproj2 = (TH1*)hn_measured_2nd->Projection(1);
        TH1 * tproj2 = (TH1*)hn_truth_2nd->Projection(3);
        TH2 * rproj2 = (TH2*)hn_response_2nd->Projection(3,1);
        TH1 * mproj3 = (TH1*)hn_measured_3rd->Projection(1);
        TH1 * tproj3 = (TH1*)hn_truth_3rd->Projection(3);
        TH2 * rproj3 = (TH2*)hn_response_3rd->Projection(3,1);
        TH1 * mproj4 = (TH1*)hn_measured_4th->Projection(1);
        TH1 * tproj4 = (TH1*)hn_truth_4th->Projection(3);
        TH2 * rproj4 = (TH2*)hn_response_4th->Projection(3,1);
        mproj2->SetName(Form("hMeas_2nd_%i", i));
        tproj2->SetName(Form("hTruth_2nd_%i", i));
        rproj2->SetName(Form("hResponse_2nd_%i", i));
        mproj3->SetName(Form("hMeas_3rd_%i", i));
        tproj3->SetName(Form("hTruth_3rd_%i", i));
        rproj3->SetName(Form("hResponse_3rd_%i", i));
        mproj4->SetName(Form("hMeas_4th_%i", i));
        tproj4->SetName(Form("hTruth_4th_%i", i));
        rproj4->SetName(Form("hResponse_4th_%i", i));
        hMeas_2nd[i] = (TH1*)mproj2->Clone(Form("hMeas_2nd_%i", i));
        hTruth_2nd[i] = (TH1*)tproj2->Clone(Form("hTruth_2nd_%i", i));
        hResponse_2nd[i] = (TH2*)rproj2->Clone(Form("hResponse_2nd_%i", i));
        hMeas_3rd[i] = (TH1*)mproj3->Clone(Form("hMeas_3rd_%i", i));
        hTruth_3rd[i] = (TH1*)tproj3->Clone(Form("hTruth_3rd_%i", i));
        hResponse_3rd[i] = (TH2*)rproj3->Clone(Form("hResponse_3rd_%i", i));
        hMeas_4th[i] = (TH1*)mproj4->Clone(Form("hMeas_4th_%i", i));
        hTruth_4th[i] = (TH1*)tproj4->Clone(Form("hTruth_4th_%i", i));
        hResponse_4th[i] = (TH2*)rproj4->Clone(Form("hResponse_4th_%i", i));
        delete mproj2;
        delete tproj2;
        delete rproj2;
        delete mproj3;
        delete tproj3;
        delete rproj3;
        delete mproj4;
        delete tproj4;
        delete rproj4;
       
        hMeas_2nd[i]->Write();
        hTruth_2nd[i]->Write();
        hResponse_2nd[i]->Write();
        hMeas_3rd[i]->Write();
        hTruth_3rd[i]->Write();
        hResponse_3rd[i]->Write();
        hMeas_4th[i]->Write();
        hTruth_4th[i]->Write();
        hResponse_4th[i]->Write();

    }

    cout << "Creating RooUnfoldResponse" << endl;
    RooUnfoldResponse responses2[nprojections];
    RooUnfoldResponse responses3[nprojections];
    RooUnfoldResponse responses4[nprojections];
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    for(Int_t i = 1; i < nprojections; i++){
        responses2[i] = RooUnfoldResponse(hMeas_2nd[i], hTruth_2nd[i], hResponse_2nd[i], Form("responses2_%i", i), Form("responses2_%i", i));
        RooUnfoldBayes unfold2(&responses2[i], hMeas_2nd[i], 7);
        hUnfolded_2nd[i] = (TH1D*)unfold2.Hunfold(errorTreatment);
        hUnfolded_2nd[i]->SetName(Form("hUnfolded_2nd_%i", i));
        hUnfolded_2nd[i]->SetTitle(Form("hUnfolded_2nd_%i", i));
        hUnfolded_2nd[i]->Write();

        responses3[i] = RooUnfoldResponse(hMeas_3rd[i], hTruth_3rd[i], hResponse_3rd[i], Form("responses3_%i", i), Form("responses3_%i", i));
        RooUnfoldBayes unfold3(&responses3[i], hMeas_3rd[i], 7);
        hUnfolded_3rd[i] = (TH1D*)unfold3.Hunfold(errorTreatment);
        hUnfolded_3rd[i]->SetName(Form("hUnfolded_3rd_%i", i));
        hUnfolded_3rd[i]->SetTitle(Form("hUnfolded_3rd_%i", i));

        responses4[i] = RooUnfoldResponse(hMeas_4th[i], hTruth_4th[i], hResponse_4th[i], Form("responses4_%i", i), Form("responses4_%i", i));
        RooUnfoldBayes unfold4(&responses4[i], hMeas_4th[i], 7);
        hUnfolded_4th[i] = (TH1D*)unfold4.Hunfold(errorTreatment);
        hUnfolded_4th[i]->SetName(Form("hUnfolded_4th_%i", i));
        hUnfolded_4th[i]->SetTitle(Form("hUnfolded_4th_%i", i));
        hUnfolded_4th[i]->Write();

    }

    TH2D * h2_measured_2nd = new TH2D("h2_measured_2nd", "h2_measured_2nd", n_reco_bins, &reco_bins_vec[0], 100, min_2nd_all, max_2nd_all);
    TH2D * h2_truth_2nd = new TH2D("h2_truth_2nd", "h2_truth_2nd", n_reco_bins, &reco_bins_vec[0], 100, min_2nd_all, max_2nd_all);
    TH2D * h2_unfolded_2nd = new TH2D("h2_unfolded_2nd", "h2_unfolded_2nd", n_reco_bins, &reco_bins_vec[0], 100, min_2nd_all, max_2nd_all);

    TH2D * h2_measured_3rd = new TH2D("h2_measured_3rd", "h2_measured_3rd", n_reco_bins, &reco_bins_vec[0], 100, min_3rd_all, max_3rd_all);
    TH2D * h2_truth_3rd = new TH2D("h2_truth_3rd", "h2_truth_3rd", n_reco_bins, &reco_bins_vec[0], 100, min_3rd_all, max_3rd_all);
    TH2D * h2_unfolded_3rd = new TH2D("h2_unfolded_3rd", "h2_unfolded_3rd", n_reco_bins, &reco_bins_vec[0], 100, min_3rd_all, max_3rd_all);

    TH2D * h2_measured_4th = new TH2D("h2_measured_4th", "h2_measured_4th", n_reco_bins, &reco_bins_vec[0], 100, min_4th_all, max_4th_all);
    TH2D * h2_truth_4th = new TH2D("h2_truth_4th", "h2_truth_4th", n_reco_bins, &reco_bins_vec[0], 100, min_4th_all, max_4th_all);
    TH2D * h2_unfolded_4th = new TH2D("h2_unfolded_4th", "h2_unfolded_4th", n_reco_bins, &reco_bins_vec[0], 100, min_4th_all, max_4th_all);



    for (Int_t i = 1; i < nprojections; i++){
        for (Int_t j = 1; j < 101; j++){
            h2_measured_2nd->SetBinContent(i, j, hMeas_2nd[i]->GetBinContent(j));
            h2_truth_2nd->SetBinContent(i, j, hTruth_2nd[i]->GetBinContent(j));
            h2_unfolded_2nd->SetBinContent(i, j, hUnfolded_2nd[i]->GetBinContent(j));
            h2_measured_2nd->SetBinError(i, j, hMeas_2nd[i]->GetBinError(j));
            h2_truth_2nd->SetBinError(i, j, hTruth_2nd[i]->GetBinError(j));
            h2_unfolded_2nd->SetBinError(i, j, hUnfolded_2nd[i]->GetBinError(j));
            h2_measured_3rd->SetBinContent(i, j, hMeas_3rd[i]->GetBinContent(j));
            h2_truth_3rd->SetBinContent(i, j, hTruth_3rd[i]->GetBinContent(j));
            h2_unfolded_3rd->SetBinContent(i, j, hUnfolded_3rd[i]->GetBinContent(j));
            h2_measured_3rd->SetBinError(i, j, hMeas_3rd[i]->GetBinError(j));
            h2_truth_3rd->SetBinError(i, j, hTruth_3rd[i]->GetBinError(j));
            h2_unfolded_3rd->SetBinError(i, j, hUnfolded_3rd[i]->GetBinError(j));
            h2_measured_4th->SetBinContent(i, j, hMeas_4th[i]->GetBinContent(j));
            h2_truth_4th->SetBinContent(i, j, hTruth_4th[i]->GetBinContent(j));
            h2_unfolded_4th->SetBinContent(i, j, hUnfolded_4th[i]->GetBinContent(j));
            h2_measured_4th->SetBinError(i, j, hMeas_4th[i]->GetBinError(j));
            h2_truth_4th->SetBinError(i, j, hTruth_4th[i]->GetBinError(j));
            h2_unfolded_4th->SetBinError(i, j, hUnfolded_4th[i]->GetBinError(j));


        }
    }

    h2_measured_2nd->Write();
    h2_truth_2nd->Write();
    h2_unfolded_2nd->Write();
    h2_measured_3rd->Write();
    h2_truth_3rd->Write();
    h2_unfolded_3rd->Write();
    h2_measured_4th->Write();
    h2_truth_4th->Write();
    h2_unfolded_4th->Write();

    TProfile * p2_measured_2nd = h2_measured_2nd->ProfileX("p2_measured_2nd");
    TProfile * p2_truth_2nd = h2_truth_2nd->ProfileX("p2_truth_2nd");
    TProfile * p2_unfolded_2nd = h2_unfolded_2nd->ProfileX("p2_unfolded_2nd");
    TProfile * p2_measured_3rd = h2_measured_3rd->ProfileX("p2_measured_3rd");
    TProfile * p2_truth_3rd = h2_truth_3rd->ProfileX("p2_truth_3rd");
    TProfile * p2_unfolded_3rd = h2_unfolded_3rd->ProfileX("p2_unfolded_3rd");
    TProfile * p2_measured_4th = h2_measured_4th->ProfileX("p2_measured_4th");
    TProfile * p2_truth_4th = h2_truth_4th->ProfileX("p2_truth_4th");
    TProfile * p2_unfolded_4th = h2_unfolded_4th->ProfileX("p2_unfolded_4th");



    p2_measured_2nd->Write();
    p2_truth_2nd->Write();
    p2_unfolded_2nd->Write();
    p2_measured_3rd->Write();
    p2_truth_3rd->Write();
    p2_unfolded_3rd->Write();
    p2_measured_4th->Write();
    p2_truth_4th->Write();
    p2_unfolded_4th->Write();


    TFile ref_flow_file(ref_input, "READ");
    TTreeReader ref_flow_reader("tree", &ref_flow_file);
    TTreeReaderValue<Float_t> in_avgM(ref_flow_reader, "avgM");
    TTreeReaderValue<Float_t> in_v2(ref_flow_reader, "v2");
    TTreeReaderValue<Float_t> in_v3(ref_flow_reader, "v3");
    TTreeReaderValue<Float_t> in_v4(ref_flow_reader, "v4");
    TTreeReaderValue<Float_t> in_v2_stat(ref_flow_reader, "v2_stat");
    TTreeReaderValue<Float_t> in_v3_stat(ref_flow_reader, "v3_stat");
    TTreeReaderValue<Float_t> in_v4_stat(ref_flow_reader, "v4_stat");
    TTreeReaderValue<Float_t> in_v2_sys(ref_flow_reader, "v2_sys");
    TTreeReaderValue<Float_t> in_v3_sys(ref_flow_reader, "v3_sys");
    TTreeReaderValue<Float_t> in_v4_sys(ref_flow_reader, "v4_sys");
    TTreeReaderValue<Float_t> in_v2_substat(ref_flow_reader, "v2_substat");
    TTreeReaderValue<Float_t> in_v3_substat(ref_flow_reader, "v3_substat");
    TTreeReaderValue<Float_t> in_v4_substat(ref_flow_reader, "v4_substat");
    ref_flow_reader.Next();
    Float_t avgM = *in_avgM;
    Float_t v2 = *in_v2;
    Float_t v3 = *in_v3;
    Float_t v4 = *in_v4;
    Float_t v2_stat = *in_v2_stat;
    Float_t v3_stat = *in_v3_stat;
    Float_t v4_stat = *in_v4_stat;
    Float_t v2_sys = *in_v2_sys;
    Float_t v3_sys = *in_v3_sys;
    Float_t v4_sys = *in_v4_sys;
    Float_t v2_substat = *in_v2_substat;
    Float_t v3_substat = *in_v3_substat;
    Float_t v4_substat = *in_v4_substat;
    ref_flow_file.Close();

    fout->cd();
    TH1D * measured_v2 = (TH1D*)p2_measured_2nd->Clone("measured_v2");
    TH1D * truth_v2 = (TH1D*)p2_truth_2nd->Clone("truth_v2");
    TH1D * unfolded_v2 = (TH1D*)p2_unfolded_2nd->Clone("unfolded_v2");
    measured_v2->SetName("measured_v2");
    measured_v2->SetTitle("measured_v2");
    truth_v2->SetName("truth_v2");
    truth_v2->SetTitle("truth_v2");
    unfolded_v2->SetName("unfolded_v2");
    unfolded_v2->SetTitle("unfolded_v2");
    measured_v2->Scale(1.0/(v2));
    truth_v2->Scale(1.0/(v2));
    unfolded_v2->Scale(1.0/(v2));
    measured_v2->Write();
    truth_v2->Write();
    unfolded_v2->Write();

    TH1D * measured_v3 = (TH1D*)p2_measured_3rd->Clone("measured_v3");
    TH1D * truth_v3 = (TH1D*)p2_truth_3rd->Clone("truth_v3");
    TH1D * unfolded_v3 = (TH1D*)p2_unfolded_3rd->Clone("unfolded_v3");
    measured_v3->SetName("measured_v3");
    measured_v3->SetTitle("measured_v3");
    truth_v3->SetName("truth_v3");
    truth_v3->SetTitle("truth_v3");
    unfolded_v3->SetName("unfolded_v3");
    unfolded_v3->SetTitle("unfolded_v3");
    measured_v3->Scale(1.0/(v3));
    truth_v3->Scale(1.0/(v3));
    unfolded_v3->Scale(1.0/(v3));
    measured_v3->Write();
    truth_v3->Write();
    unfolded_v3->Write();

    TH1D * measured_v4 = (TH1D*)p2_measured_4th->Clone("measured_v4");
    TH1D * truth_v4 = (TH1D*)p2_truth_4th->Clone("truth_v4");
    TH1D * unfolded_v4 = (TH1D*)p2_unfolded_4th->Clone("unfolded_v4");
    measured_v4->SetName("measured_v4");
    measured_v4->SetTitle("measured_v4");
    truth_v4->SetName("truth_v4");
    truth_v4->SetTitle("truth_v4");
    unfolded_v4->SetName("unfolded_v4");
    unfolded_v4->SetTitle("unfolded_v4");
    measured_v4->Scale(1.0/(v4));
    truth_v4->Scale(1.0/(v4));
    unfolded_v4->Scale(1.0/(v4));
    measured_v4->Write();
    truth_v4->Write();
    unfolded_v4->Write();


    
    cout << "Creating RooUnfoldResponse" << endl;
    fout->Close();
    fin.Close();
    return 0;

}
#ifndef __CINT__
Int_t main(Int_t argc, Char_t** argv){
    cout << "Starting " << argv[0] << endl;
    // get arguments
    // check if correct number of arguments
    if(argc != 4){ cout << "Incorrect number of arguments" << endl; exit(-1);}
   // get arguments
    TString diff_input = argv[1];
    TString ref_input = argv[2];
    TString output_file = argv[3];
    // print arguments
    cout << "diff_input: " << diff_input << endl;
    cout << "ref_input: " << ref_input << endl;
    cout << "output_file: " << output_file << endl;

    return Unfoldvn(diff_input, ref_input, output_file);

}
#endif