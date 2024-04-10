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
const Int_t n_iterations_min = 4;
const Int_t n_iterations_max = 8;
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

void Unfold2D(TString inputfile, TString responsefile, TString outputfile){
    
    cout << "Unfolding " << inputfile <<  " with response " << responsefile << endl;
    cout << "Output file " << outputfile << endl;

    // get min and max of correlations
    Double_t min_measured_correlation = GetMin("accepted_matched_correlation", inputfile);
    Double_t max_measured_correlation = GetMax("accepted_matched_correlation", inputfile);
    Double_t min_truth_correlation = GetMin("accepted_truth_correlation", inputfile);
    Double_t max_truth_correlation = GetMax("accepted_truth_correlation", inputfile);
    Double_t min_correlation = TMath::Min(min_measured_correlation, min_truth_correlation);
    Double_t max_correlation = TMath::Max(max_measured_correlation, max_truth_correlation);
    // min_measured_correlation= min_correlation;
    // max_measured_correlation = max_correlation;
    // min_truth_correlation = min_correlation;
    // max_truth_correlation = max_correlation;


    // min max pt 
    Double_t min_measured_pt = GetMin("accepted_matched_jet_pt", inputfile);
    Double_t max_measured_pt = GetMax("accepted_matched_jet_pt", inputfile);
    Double_t min_truth_pt = GetMin("accepted_truth_jet_pt", inputfile);
    Double_t max_truth_pt = GetMax("accepted_truth_jet_pt", inputfile);

    // min max jet_vn_mc
    Double_t min_jet_vn_mc = GetMin("jet_vn_mc", inputfile);
    Double_t max_jet_vn_mc = GetMax("jet_vn_mc", inputfile);
    min_jet_vn_mc -= 0.1*min_jet_vn_mc;
    max_jet_vn_mc += 0.1*max_jet_vn_mc;
   
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
    cout << "Measured pt min: " << reco_pt_min << endl;
    cout << "Measured pt max: " << reco_pt_max << endl;
    cout << "Truth pt min: " << truth_pt_min << endl;
    cout << "Truth pt max: " << truth_pt_max << endl;
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
    TTreeReaderValue<Double_t> in_jet_vn_mc(reader, "jet_vn_mc");
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
    TH2D * hMeasuredCorrelation = new TH2D("hMeasuredCorrelation", "Measured Correlation", n_reco_pt_bins, &pt_bin_edges_reco_vec[0], 100, min_measured_correlation, max_measured_correlation);
    TH2D * hTruthCorrelation = new TH2D("hTruthCorrelation", "Truth Correlation", n_truth_pt_bins, &pt_bin_edges_truth_vec[0], 100, min_truth_correlation, max_truth_correlation);
    TH2D * hFakeCorrelation = new TH2D("hFakeCorrelation", "Fake Correlation", n_reco_pt_bins, &pt_bin_edges_reco_vec[0], 100, min_measured_correlation, max_measured_correlation);

    // debug response histos
    TH2D * hResponseMeasuredCorrelation = new TH2D("hResponseMeasuredCorrelation", "Response Measured Correlation", n_reco_pt_bins, &pt_bin_edges_reco_vec[0], 100, min_measured_correlation, max_measured_correlation);
    TH2D * hResponseTruthCorrelation = new TH2D("hResponseTruthCorrelation", "Response Truth Correlation", n_truth_pt_bins, &pt_bin_edges_truth_vec[0], 100, min_truth_correlation, max_truth_correlation);

    TH2D * hJetVnMC_truth = new TH2D("hJetVnMC_truth", "Jet Vn MC Truth", n_truth_pt_bins, &pt_bin_edges_truth_vec[0], 100, min_jet_vn_mc, max_jet_vn_mc);
    TH2D * hJetVnMC = new TH2D("hJetVnMC", "Jet Vn MC", n_reco_pt_bins, &pt_bin_edges_reco_vec[0], 100, min_jet_vn_mc, max_jet_vn_mc);
    TH2D * hUnfolded[n_iterations];

    // set up roounfold
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse *response_object = new RooUnfoldResponse("response_object", "");
    response_object->Setup(hMeasuredCorrelation, hTruthCorrelation);

    cout << "Filling measured and truth histos" << endl;
    // loop through events in input file and fill truth and measured histos
    while(reader.Next()){
        // fill measured histos
        Double_t x = *in_accepted_matched_jet_pt;
        Double_t y = *in_accepted_matched_correlation;
        Double_t w = 1.0;
        Double_t xt = *in_accepted_truth_jet_pt;
        Double_t yt = *in_accepted_truth_correlation;
        Double_t xm = *in_missed_truth_jet_pt;
        Double_t ym = *in_missed_truth_correlation;
        if(xt > 0.0){
            hMeasuredCorrelation->Fill(x, y);
            hTruthCorrelation->Fill(xt, yt);
        }
        if(xm > 0.0){
            hTruthCorrelation->Fill(xm, ym);
        }

        for (Int_t i = 0; i < *in_njets; i++){
            Double_t fx = in_fake_jet_pt[i];
            Double_t fy = in_fake_jet_correlation[i];
            hFakeCorrelation->Fill(fx, fy);
            hMeasuredCorrelation->Fill(fx, fy);
        }

        hJetVnMC_truth->Fill(*in_accepted_truth_jet_pt, *in_jet_vn_mc);
        hJetVnMC->Fill(*in_accepted_matched_jet_pt, *in_jet_vn_mc);

    }

    cout << "Filling response object" << endl;
    // loop through events in response file and fill response histos
    Int_t event_counter = 0;
    while(reader_response.Next()){
       
        Double_t xt = *in_response_accepted_truth_jet_pt;
        Double_t x = *in_response_accepted_matched_jet_pt;
        Double_t y = *in_response_accepted_matched_correlation;
        Double_t yt = *in_response_accepted_truth_correlation;
        Double_t w = 1.0;
        Double_t xm = *in_response_missed_truth_jet_pt;
        Double_t ym = *in_response_missed_truth_correlation;
        if(xt > 0.0){
            response_object->Fill(x, y, xt, yt, w);
            hResponseMeasuredCorrelation->Fill(x, y);
            hResponseTruthCorrelation->Fill(xt, yt);
        }
        if(xm > 0.0){
            response_object->Miss(xm, ym, w);
            hResponseTruthCorrelation->Fill(xm, ym);
        }

        for (Int_t i = 0; i < *in_response_njets; i++){
            Double_t fx = in_response_fake_jet_pt[i];
            Double_t fy = in_response_fake_jet_correlation[i];
            response_object->Fake(fx, fy, w);
            hResponseMeasuredCorrelation->Fill(fx, fy);
        }
    }
   
    cout << "Unfolding" << endl;
    // unfold
    for(Int_t i = 0; i < n_iterations; i++){
        cout << "Iteration " << i+n_iterations_min << endl;
        RooUnfoldBayes unfold_object(response_object, hMeasuredCorrelation, i+n_iterations_min);
        hUnfolded[i] = (TH2D*) unfold_object.Hunfold(errorTreatment);
        hUnfolded[i]->SetName(Form("hUnfolded_iter%i", i+n_iterations_min));
        hUnfolded[i]->SetTitle(Form("Unfolded Iteration %i", i+n_iterations_min));
        hUnfolded[i]->Write();
    }
    cout << "Done unfolding" << endl;
    // get profiles of all iterations
    for (Int_t i = 0; i < n_iterations; i++){
        TProfile * hUnfoldedProfile = hUnfolded[i]->ProfileX(Form("hProfile_iter%i", i+n_iterations_min));
        TH1D * hUnfolded1D = hUnfoldedProfile->ProjectionX(Form("hUnfolded1D_iter%i", i+n_iterations_min));
        hUnfolded1D->Write();
    }


    // TProfile * hUnfoldedProfile = hUnfolded[n_iterations-2]->ProfileX("hProfile");
    // TH1D * hUnfolded1D = hUnfoldedProfile->ProjectionX("hUnfolded1D");
    // hUnfolded1D->Write();


    // TProfile * hUnfoldedProfile_nminus = hUnfolded[n_iterations-3]->ProfileX("hProfile_nminus");
    // TH1D * hUnfolded1D_nminus = hUnfoldedProfile_nminus->ProjectionX("hUnfolded1D_nminus");
    // hUnfolded1D_nminus->Write();

    // TProfile * hUnfoldedProfile_nplus = hUnfolded[n_iterations-1]->ProfileX("hProfile_nplus");
    // TH1D * hUnfolded1D_nplus = hUnfoldedProfile_nplus->ProjectionX("hUnfolded1D_nplus");
    // hUnfolded1D_nplus->Write();

    // get measured 
    TProfile * hMeasuredProfile = hMeasuredCorrelation->ProfileX("hMeasuredProfile");
    TH1D * hMeasured1D = hMeasuredProfile->ProjectionX("hMeasured1D");
    hMeasured1D->Write();

    // get truth
    TProfile * hTruthProfile = hTruthCorrelation->ProfileX("hTruthProfile");
    TH1D * hTruth1D = hTruthProfile->ProjectionX("hTruth1D");
    hTruth1D->Write();

    // get fake
    TProfile * hFakeProfile = hFakeCorrelation->ProfileX("hFakeProfile");
    TH1D * hFake1D = hFakeProfile->ProjectionX("hFake1D");
    hFake1D->Write();

    // get vn 
    TProfile * hJetVnMC_truthProfile = hJetVnMC_truth->ProfileX("hJetVnMC_truthProfile");
    TH1D * hJetVnMC_truth1D = hJetVnMC_truthProfile->ProjectionX("hJetVnMC_truth1D");
    hJetVnMC_truth1D->Write();

    TProfile * hJetVnMCProfile = hJetVnMC->ProfileX("hJetVnMCProfile");
    TH1D * hJetVnMC1D = hJetVnMCProfile->ProjectionX("hJetVnMC1D");
    hJetVnMC1D->Write();

    hJetVnMC_truth->Write();
    hJetVnMC->Write();


    // write histos
    
    fout->cd();
    hResponseMeasuredCorrelation->Write();
    hResponseTruthCorrelation->Write();
    hMeasuredCorrelation->Write();
    hTruthCorrelation->Write();
    hFakeCorrelation->Write();
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
    Unfold2D(inputfile, responsefile, outputfile);
    return 0;
}
#endif