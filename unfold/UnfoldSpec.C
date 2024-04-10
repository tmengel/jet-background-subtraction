#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;

#include "TH1D.h"
#include "TH2D.h"
#include "TH2.h"
#include "TH1.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TList.h"

#include "TFile.h"
#include "TTree.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldErrors.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"


#include "TString.h"
#include "TMath.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include <vector>


#endif

using namespace std;
using namespace ROOT;

// reco jet pt bins

// unfolding iteratons
const Int_t n_iterations_min = 7;
const Int_t n_iterations_max = 9;
Int_t n_iterations = n_iterations_max-n_iterations_min+1;

Double_t GetChi2(TH1D *Unfolded, TH1D *Truth)
{
  Double_t chi2 = 0;
  for(Int_t i = 1; i <= Unfolded->GetNbinsX(); i++)
  {
    chi2 += TMath::Power((Truth->GetBinContent(i)-Unfolded->GetBinContent(i))/Unfolded->GetBinError(i),2.);
  }
  return chi2;
}

Double_t GetMax(TString Variable, TString FileName){
    cout << "Getting max of " << Variable.Data() << " from " << FileName.Data() << endl;
    ROOT::EnableImplicitMT();
    TString treeName = "tree";
    RDataFrame df(treeName.Data(), FileName.Data());
    Double_t max = df.Max(Variable.Data()).GetValue();
    ROOT::DisableImplicitMT();
    return max;
}

Double_t GetMin(TString Variable, TString FileName){
    cout << "Getting min of " << Variable.Data() << " from " << FileName.Data() << endl;
    ROOT::EnableImplicitMT();
    TString treeName = "tree";
    RDataFrame df(treeName.Data(), FileName.Data());
    Double_t min = df.Min(Variable.Data()).GetValue();
    ROOT::DisableImplicitMT();
    return min;
}

void UnfoldSpec(TString inputfile, TString responsefile, TString outputfile){
    
    cout << "Unfolding " << inputfile <<  " with response " << responsefile << endl;
    cout << "Output file " << outputfile << endl;

    // min max pt 
    Double_t min_measured_pt = GetMin("accepted_matched_jet_pt", inputfile);
    Double_t max_measured_pt = GetMax("accepted_matched_jet_pt", inputfile);
    Double_t min_truth_pt = GetMin("accepted_truth_jet_pt", inputfile);
    Double_t max_truth_pt = GetMax("accepted_truth_jet_pt", inputfile);
   

    // round down min and round up max to nearest 10
    Double_t pt_bin_width = 5.0;
    Double_t reco_pt_min = TMath::Floor(min_measured_pt/pt_bin_width)*pt_bin_width;
    Double_t reco_pt_max = TMath::Ceil(max_measured_pt/pt_bin_width)*pt_bin_width;
    Double_t truth_pt_min = TMath::Floor(min_truth_pt/pt_bin_width)*pt_bin_width;
    Double_t truth_pt_max = TMath::Ceil(max_truth_pt/pt_bin_width)*pt_bin_width;
    reco_pt_min = 5.0;
    truth_pt_min = 10.0;
    reco_pt_max = 90.0;
    truth_pt_max = 90.0;
    Int_t n_reco_pt_bins = Int_t((reco_pt_max-reco_pt_min)/pt_bin_width);
    Int_t n_truth_pt_bins = Int_t((truth_pt_max-truth_pt_min)/pt_bin_width);
   
    const Int_t n_pt_bin_edges_reco = n_reco_pt_bins+1;
    const Int_t n_pt_bin_edges_truth = n_truth_pt_bins+1;
    std::vector<Double_t> pt_bin_edges_reco_vec;
    std::vector<Double_t> pt_bin_edges_truth_vec;
    for(Int_t i = 0; i < n_pt_bin_edges_reco; i++){
        pt_bin_edges_reco_vec.push_back(reco_pt_min+i*pt_bin_width);
    }
    for(Int_t i = 0; i < n_pt_bin_edges_truth; i++){
        pt_bin_edges_truth_vec.push_back(truth_pt_min+i*pt_bin_width);
    }

   
    TFile fin(inputfile.Data(), "READ");
    TFile fresponse(responsefile.Data(), "READ");
    if(!fin.IsOpen() || !fresponse.IsOpen()){
        cout << "Could not open input file or response file" << endl; exit(-1);
    }
    // set up tree reader
    TTreeReader reader("tree", &fin); // get measured and truth from input file
    TTreeReaderValue<Int_t> in_event_rfp_M(reader, "event_rfp_M");
    TTreeReaderValue<Double_t> in_accepted_truth_jet_pt(reader, "accepted_truth_jet_pt");
    TTreeReaderValue<Double_t> in_accepted_truth_correlation(reader, "accepted_truth_correlation");
    TTreeReaderValue<Double_t> in_accepted_matched_jet_pt(reader, "accepted_matched_jet_pt");
    TTreeReaderValue<Double_t> in_accepted_matched_correlation(reader, "accepted_matched_correlation");
    TTreeReaderValue<Double_t> in_missed_truth_jet_pt(reader, "missed_truth_jet_pt");
    TTreeReaderValue<Double_t> in_missed_truth_correlation(reader, "missed_truth_correlation");
    TTreeReaderValue<Int_t> in_njets(reader, "njets");
    TTreeReaderArray<Double_t> in_fake_jet_pt(reader, "fake_jet_pt");
    TTreeReaderArray<Double_t> in_fake_jet_correlation(reader, "fake_jet_correlation");

    TTreeReader reader_response("tree", &fresponse); // get accepted truth, accepted matched, and missed truth and fake from response file
    TTreeReaderValue<Int_t> in_response_event_rfp_M(reader_response, "event_rfp_M");
    TTreeReaderValue<Double_t> in_response_accepted_truth_jet_pt(reader_response, "accepted_truth_jet_pt");
    TTreeReaderValue<Double_t> in_response_accepted_truth_correlation(reader_response, "accepted_truth_correlation");
    TTreeReaderValue<Double_t> in_response_accepted_matched_jet_pt(reader_response, "accepted_matched_jet_pt");
    TTreeReaderValue<Double_t> in_response_accepted_matched_correlation(reader_response, "accepted_matched_correlation");
    TTreeReaderValue<Double_t> in_response_missed_truth_jet_pt(reader_response, "missed_truth_jet_pt");
    TTreeReaderValue<Double_t> in_response_missed_truth_correlation(reader_response, "missed_truth_correlation");
    TTreeReaderValue<Int_t> in_response_njets(reader_response, "njets");
    TTreeReaderArray<Double_t> in_response_fake_jet_pt(reader_response, "fake_jet_pt");
    TTreeReaderArray<Double_t> in_response_fake_jet_correlation(reader_response, "fake_jet_correlation");


    // set up output file
    TFile * fout = new TFile(outputfile.Data(), "RECREATE");
    // create TH2Ds for correlations
    TH1D * hMeasuredJetPt = new TH1D("hMeasuredJetPt", "Measured Jet Pt", n_reco_pt_bins, &pt_bin_edges_reco_vec[0]);
    TH1D * hMissedJetPt = new TH1D("hMissedJetPt", "Missed Jet Pt", n_truth_pt_bins, &pt_bin_edges_truth_vec[0]);
    TH1D * hTruthJetPt = new TH1D("hTruthJetPt", "Truth Jet Pt", n_truth_pt_bins, &pt_bin_edges_truth_vec[0]);
    TH1D * hFakeJetPt = new TH1D("hFakeJetPt", "Fake Jet Pt", n_reco_pt_bins, &pt_bin_edges_reco_vec[0]);
    TH1D * hEvents = new TH1D("hEvents", "Events", 1, 0.0, 1.0);
    // debug response histos
    TH1D * hResponseMeasuredJetPt = new TH1D("hResponseMeasuredJetPt", "Response Measured Jet Pt", n_reco_pt_bins, &pt_bin_edges_reco_vec[0]);
    TH1D * hResponseTruthJetPt = new TH1D("hResponseTruthJetPt", "Response Truth Jet Pt", n_truth_pt_bins, &pt_bin_edges_truth_vec[0]);

    TH1D * hUnfoldedJetPt[n_iterations];

    // set up roounfold
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse response_momentum_object = RooUnfoldResponse(n_reco_pt_bins, reco_pt_min, reco_pt_max, n_truth_pt_bins, truth_pt_min, truth_pt_max);

    cout << "Filling measured and truth histos" << endl;
    // loop through events in input file and fill truth and measured histos
    Int_t event_counter = 0;
       while(reader.Next()){
        // fill measured histos
        Double_t x = *in_accepted_matched_jet_pt;
        Double_t w = 1.0;
        Double_t xt = *in_accepted_truth_jet_pt;
        Double_t xm = *in_missed_truth_jet_pt;

        if(xt > 0.0){
            hMeasuredJetPt->Fill(x);
            hTruthJetPt->Fill(xt);

        }
        if(xm > 0.0){
            hTruthJetPt->Fill(xm);
            hMissedJetPt->Fill(xm);
        }

        for (Int_t i = 0; i < *in_njets; i++){
            Double_t fx = in_fake_jet_pt[i];
            hMeasuredJetPt->Fill(fx);
            hFakeJetPt->Fill(fx);
        }
        event_counter++;

    }
    Double_t nevents = 1.0*event_counter;
    hEvents->SetBinContent(1, nevents);

    cout << "Filling response object" << endl;
    // loop through events in response file and fill response histos
    event_counter = 0;
    while(reader_response.Next()){
       
        Double_t xt = *in_response_accepted_truth_jet_pt;
        Double_t x = *in_response_accepted_matched_jet_pt;
        Double_t xm = *in_response_missed_truth_jet_pt;
        if(xt > 0.0){
            response_momentum_object.Fill(x, xt);
            hResponseMeasuredJetPt->Fill(x);
            hResponseTruthJetPt->Fill(xt);
        }
        if(xm > 0.0){
            response_momentum_object.Miss(xm);
            hResponseTruthJetPt->Fill(xm);
        }

        for (Int_t i = 0; i < *in_response_njets; i++){
            Double_t fx = in_response_fake_jet_pt[i];
            response_momentum_object.Fake(fx); 
            hResponseMeasuredJetPt->Fill(fx);
        }
    }
   
  
    cout << "Unfolding momentum" << endl;
    for(Int_t i = 0; i < n_iterations; i++){
        cout << "Iteration " << i+n_iterations_min << endl;
        RooUnfoldBayes unfold_object(&response_momentum_object, hMeasuredJetPt, i+n_iterations_min);
        hUnfoldedJetPt[i] = (TH1D*) unfold_object.Hunfold(errorTreatment);
        hUnfoldedJetPt[i]->SetName(Form("hUnfoldedJetPt_iter%i", i+n_iterations_min));
        hUnfoldedJetPt[i]->SetTitle(Form("Unfolded Jet Pt Iteration %i", i+n_iterations_min));
        hUnfoldedJetPt[i]->Write();
    }

    cout << "Done unfolding" << endl;

    fout->cd();
    hEvents->Write();
    hResponseMeasuredJetPt->Write();
    hResponseTruthJetPt->Write();   
    hMeasuredJetPt->Write();
    hTruthJetPt->Write();
    hFakeJetPt->Write();
    hMissedJetPt->Write();
    // clean up
    fin.Close();
    fresponse.Close();

}

#ifndef __CINT__
Int_t main(Int_t argc, Char_t** argv){

    cout << "Starting " << argv[0] << endl;
    // get arguments
    // check if correct number of arguments
    if(argc != 4){
        cout << "Usage: " << argv[0] << " <inputfile> <response> <output> " << endl;
        cout << "Incorrect number of arguments" << endl; exit(-1);
    }
   // get arguments
    TString inputfile = argv[1];
    TString responsefile = argv[2];
    TString outputfile = argv[3];
    UnfoldSpec(inputfile, responsefile, outputfile);
    return 0;
}
#endif