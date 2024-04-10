#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)
#include <iostream>
using std::cout;
using std::endl;


#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TSystem.h"
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldErrors.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDF/RInterface.hxx"

#include <iostream>
#include <string>
#include <vector>

#endif

using namespace std;
using namespace ROOT;

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

Float_t GetChi2(TH1D *Unfolded, TH1D *Truth)
{
  Float_t chi2 = 0;
  for(Int_t i = 1; i <= Unfolded->GetNbinsX(); i++)
  {
    chi2 += TMath::Power((Truth->GetBinContent(i)-Unfolded->GetBinContent(i))/Unfolded->GetBinError(i),2.);
  }
  return chi2;
}

TMatrixD CalculatePearsonCoefficients(TMatrixD covmat)
{
  Int_t nrows = covmat.GetNrows();
  Int_t ncols = covmat.GetNcols();

  TMatrixD pearsonCoefs(nrows,ncols);

  Float_t pearson = 0.;

  for(Int_t row = 0; row<nrows; row++)
  {
    for(Int_t col = 0; col<ncols; col++)
    {
      if(covmat(row,row)!=0. && covmat(col,col)!=0.)
      {
        pearson = covmat(row,col)/TMath::Sqrt(covmat(row,row)*covmat(col,col));
      }
      pearsonCoefs(row,col) = pearson;
    }
  }

  return pearsonCoefs;
}

Int_t Accepted(Float_t pT, Float_t pTmin, Float_t pTmax){
  if(pT>pTmin && pT<pTmax) return 1;
  else return 0;
}

Float_t GetSigmaPt(TString TargetName, TString TruthName, TString FileName){

    cout << "Getting sigma of " << TargetName.Data() << " from " << FileName.Data() << endl;
  
    ROOT::EnableImplicitMT();
    TString treeName = "tree";
    RDataFrame df(treeName.Data(), FileName.Data());
    Float_t sigma = df.Define("deltaPt", Form("%s-%s", TargetName.Data(), TruthName.Data())).Histo1D({"hDeltaPt","hDeltaPt", 100, -49.0, 49.0}, "deltaPt")->GetRMS();
    ROOT::DisableImplicitMT();
    return sigma;
}

Float_t GetMaxPt(TString Variable, TString FileName){

    cout << "Getting max of " << Variable.Data() << " from " << FileName.Data() << endl;
    ROOT::EnableImplicitMT();
    TString treeName = "tree";
    RDataFrame df(treeName.Data(), FileName.Data());
    Float_t maxPt = df.Max(Variable.Data()).GetValue();
    ROOT::DisableImplicitMT();
    return maxPt;
}

Float_t GetMinPt(TString Variable, TString FileName){

    cout << "Getting min of " << Variable.Data() << " from " << FileName.Data() << endl;
    ROOT::EnableImplicitMT();
    TString treeName = "tree";
    RDataFrame df(treeName.Data(), FileName.Data());
    Float_t minPt = df.Min(Variable.Data()).GetValue();
    ROOT::DisableImplicitMT();
    return minPt;
}

Float_t GetBinWidth(TString Variable, TString FileName){

  cout << "Getting bin width of " << Variable.Data() << " from " << FileName.Data() << endl;
  ROOT::EnableImplicitMT();
  TString treeName = "tree";
  RDataFrame df(treeName.Data(), FileName.Data());
  Float_t std = df.StdDev(Variable.Data()).GetValue();
  Int_t nSamples = df.Count().GetValue();
  Float_t binWidth = (3.49*std)/TMath::Power(nSamples, 0.333)*20;

  ROOT::DisableImplicitMT();
  return binWidth;

}

Int_t Unfold(const TString matched_file, const TString unmatched_file , const TString missed_file, const TString fake_file, const TString truth_file, const TString output_file)
{


  Int_t nIter_low = 2;
  Int_t nIter_high = 9;
  Int_t totalIter = nIter_high - nIter_low + 1;

  std::vector<TString> models = {"Area", "Mult", "DNN", "SNN"};
  std::vector<TString> target_variables = {"jet_pt_area_corrected", "jet_pt_multiplicity_corrected", "jet_pt_dnn_corrected", "jet_pt_snn_corrected"};
  TString truth_variable = "jet_pt_truth";
  TString truth_variable_pp = "jet_pt";

  Float_t min_jet_pT_measurement[4];
  Float_t max_jet_pT_measurement[4];
  Int_t n_measurement_bins[4];
  Float_t measurement_bin_width[4];

  Float_t sigma_delta_pT[4];

  cout << "Truth binning" << endl;
  Float_t min_jet_pT_truth = 10.0;
  Float_t max_jet_pT_truth = GetMaxPt(truth_variable, matched_file);
  // if (max_jet_pT_truth > 300) max_jet_pT_truth = 300;

  Float_t truth_bin_width = 5.0;// GetBinWidth(truth_variable, matched_file);
  Int_t n_truth_bins = Int_t((max_jet_pT_truth-min_jet_pT_truth)/truth_bin_width) + 1;
  max_jet_pT_truth = (n_truth_bins*truth_bin_width) + min_jet_pT_truth;
  Float_t measured_binwidth_constant =5.0; // GetBinWidth(target_variables[0], matched_file);

  for(Int_t imodel=0; imodel<models.size(); imodel++){
      cout << "Model = " << models[imodel] << endl;
      // measurement binning
      sigma_delta_pT[imodel] = GetSigmaPt(target_variables[imodel], truth_variable, matched_file);
      min_jet_pT_measurement[imodel] =  GetMinPt(target_variables[imodel], matched_file);
      //5.0*sigma_delta_pT[imodel];
      //GetMinPt(target_variables[imodel], matched_file);
      max_jet_pT_measurement[imodel] = GetMaxPt(target_variables[imodel], matched_file); // max_jet_pT_truth;//90.0; //GetMaxPt(target_variables[imodel], matched_file);
      // if (max_jet_pT_measurement[imodel] > 300) max_jet_pT_measurement[imodel] = 300;

      measurement_bin_width[imodel] = measured_binwidth_constant; //GetBinWidth(target_variables[imodel], matched_file);
      n_measurement_bins[imodel] = Int_t((max_jet_pT_measurement[imodel]-min_jet_pT_measurement[imodel])/measurement_bin_width[imodel]) + 1;
      max_jet_pT_measurement[imodel] = (n_measurement_bins[imodel]*measurement_bin_width[imodel]) + min_jet_pT_measurement[imodel];

     
  }

  //cout for debugging
  for(Int_t imodel=0; imodel<models.size(); imodel++){
      cout << "Model: " << models[imodel] << endl;
      cout << "Truth binning: " << endl;
      cout << "min_jet_pT_truth = " << min_jet_pT_truth << endl;
      cout << "max_jet_pT_truth = " << max_jet_pT_truth << endl;
      cout << "n_truth_bins = " << n_truth_bins << endl;
      cout << "truth_bin_width = " << truth_bin_width << endl;
      cout << "Measurement binning: " << endl;
      cout << "min_jet_pT_measurement = " << min_jet_pT_measurement[imodel] << endl;
      cout << "max_jet_pT_measurement = " << max_jet_pT_measurement[imodel] << endl;
      cout << "n_measurement_bins = " << n_measurement_bins[imodel] << endl;
      cout << "measurement_bin_width = " << measurement_bin_width[imodel] << endl;
      cout << "sigma_delta_pT = " << sigma_delta_pT[imodel] << endl;

  }

  // ----------------------------------------- //
  // Histogram declarations for each variable //
  
  // Unmatched histograms //
  TH1D* measurement_area_reco = new TH1D("measurement_area_reco", "measurement_area_reco", n_measurement_bins[0], min_jet_pT_measurement[0], max_jet_pT_measurement[0]);
  TH1D* measurement_mult_reco = new TH1D("measurement_mult_reco", "measurement_mult_reco", n_measurement_bins[1], min_jet_pT_measurement[1], max_jet_pT_measurement[1]);
  TH1D* measurement_dnn_reco = new TH1D("measurement_dnn_reco", "measurement_dnn_reco", n_measurement_bins[2], min_jet_pT_measurement[2], max_jet_pT_measurement[2]);
  TH1D* measurement_snn_reco = new TH1D("measurement_snn_reco", "measurement_snn_reco", n_measurement_bins[3], min_jet_pT_measurement[3], max_jet_pT_measurement[3]);
  
  // Truth histograms //
  TH1D* h_truth = new TH1D("h_truth", "h_truth", n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  TH1D* h_truth_measurement_area_bins = new TH1D("h_truth_measurement_area_bins", "h_truth_measurement_area_bins", n_measurement_bins[0], min_jet_pT_measurement[0], max_jet_pT_measurement[0]);
  TH1D* h_truth_measurement_mult_bins = new TH1D("h_truth_measurement_mult_bins", "h_truth_measurement_mult_bins", n_measurement_bins[1], min_jet_pT_measurement[1], max_jet_pT_measurement[1]);
  TH1D* h_truth_measurement_dnn_bins = new TH1D("h_truth_measurement_dnn_bins", "h_truth_measurement_dnn_bins", n_measurement_bins[2], min_jet_pT_measurement[2], max_jet_pT_measurement[2]);
  TH1D* h_truth_measurement_snn_bins = new TH1D("h_truth_measurement_snn_bins", "h_truth_measurement_snn_bins", n_measurement_bins[3], min_jet_pT_measurement[3], max_jet_pT_measurement[3]);

  TH1D* area_reco_uncut = new TH1D("area_reco_uncut", "area_reco_uncut", n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  TH1D* mult_reco_uncut = new TH1D("mult_reco_uncut", "mult_reco_uncut", n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  TH1D* dnn_reco_uncut = new TH1D("dnn_reco_uncut", "dnn_reco_uncut", n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  TH1D* snn_reco_uncut = new TH1D("snn_reco_uncut", "snn_reco_uncut", n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);

  // ratio of measured to truth
  TH1D* area_measurement_over_truth;
  TH1D* mult_measurement_over_truth;
  TH1D* dnn_measurement_over_truth;
  TH1D* snn_measurement_over_truth;

  TH1D* area_measurement_over_truth_normalized;
  TH1D* mult_measurement_over_truth_normalized;
  TH1D* dnn_measurement_over_truth_normalized;
  TH1D* snn_measurement_over_truth_normalized;

  // Unfolded histograms
  TH1D *unfolded_area[totalIter];
  TH1D *unfolded_mult[totalIter];
  TH1D *unfolded_dnn[totalIter];
  TH1D *unfolded_snn[totalIter]; 

  // ratio of unfolded to truth 
  TH1D* area_unfolded_over_truth[totalIter];
  TH1D* mult_unfolded_over_truth[totalIter];
  TH1D* dnn_unfolded_over_truth[totalIter];
  TH1D* snn_unfolded_over_truth[totalIter];

  // Pearson Correlation Coefficients
  TH2D* area_pearson[totalIter];
  TH2D* mult_pearson[totalIter];
  TH2D* dnn_pearson[totalIter];
  TH2D* snn_pearson[totalIter];

  TH2D* area_covariance[totalIter];
  TH2D* mult_covariance[totalIter];
  TH2D* dnn_covariance[totalIter];
  TH2D* snn_covariance[totalIter];

  //------------------------------------------//
  // 1. Get RooUnfoldResponse for each model  //
  RooUnfoldResponse response_area = RooUnfoldResponse(n_measurement_bins[0], min_jet_pT_measurement[0], max_jet_pT_measurement[0], n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  RooUnfoldResponse response_mult = RooUnfoldResponse(n_measurement_bins[1], min_jet_pT_measurement[1], max_jet_pT_measurement[1], n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  RooUnfoldResponse response_dnn = RooUnfoldResponse(n_measurement_bins[2], min_jet_pT_measurement[2], max_jet_pT_measurement[2], n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  RooUnfoldResponse response_snn = RooUnfoldResponse(n_measurement_bins[3], min_jet_pT_measurement[3], max_jet_pT_measurement[3], n_truth_bins, min_jet_pT_truth, max_jet_pT_truth);
  RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance; // error method


  /////////////////////////////////////////////////////////////////////
  // 2. Get RooUnfoldResponse for each model and fill with histograms //
  ///////////////////////////
  cout << "Reading trees..." << endl;
  // matched file //
  
  TFile *matched= new TFile(matched_file, "READ");
  cout << "Matched file: " << matched_file << endl;
  TTree *matchedtree = (TTree*)matched->Get("tree");
  Float_t jet_pt_truth, jet_pt_area, jet_pt_mult, jet_pt_dnn, jet_pt_snn;
  Float_t jet_track_pt_0_truth;
  matchedtree->SetBranchAddress(truth_variable, &jet_pt_truth);
  matchedtree->SetBranchAddress(target_variables[0], &jet_pt_area);
  matchedtree->SetBranchAddress(target_variables[1], &jet_pt_mult);
  matchedtree->SetBranchAddress(target_variables[2], &jet_pt_dnn);
  matchedtree->SetBranchAddress(target_variables[3], &jet_pt_snn);
  matchedtree->SetBranchAddress("jet_track_pt_0", &jet_track_pt_0_truth);


  for(Int_t i=0; i<matchedtree->GetEntries(); i++){

      matchedtree->GetEntry(i);

      
      // fill response
      // if(jet_track_pt_0_truth > 7.0){
        response_area.Fill(jet_pt_area, jet_pt_truth);
        response_mult.Fill(jet_pt_mult, jet_pt_truth);
        response_dnn.Fill(jet_pt_dnn, jet_pt_truth);
        response_snn.Fill(jet_pt_snn, jet_pt_truth);

        h_truth->Fill(jet_pt_truth);
        h_truth_measurement_area_bins->Fill(jet_pt_truth);
        h_truth_measurement_mult_bins->Fill(jet_pt_truth);
        h_truth_measurement_dnn_bins->Fill(jet_pt_truth);
        h_truth_measurement_snn_bins->Fill(jet_pt_truth);
      // }

      // if(jet_track_pt_0_truth < 5.0){
      //   //fill missed response
      //   response_area.Miss(jet_pt_truth);
      //   response_mult.Miss(jet_pt_truth);
      //   response_dnn.Miss(jet_pt_truth);
      //   response_snn.Miss(jet_pt_truth);
      // }

  }


  ///////////////////////////
  // missed file //
  TFile *missed = new TFile(missed_file, "READ");
  cout << "Missed file: " << missed_file << endl;
  TTree *missedtree = (TTree*)missed->Get("tree");
  Float_t jet_pt_truth_missed;
  missedtree->SetBranchAddress("jet_pt_truth", &jet_pt_truth_missed);
  
  for (Int_t i = 0; i < missedtree->GetEntries(); i++)
  {
      missedtree->GetEntry(i);
      // fill missed response
      response_area.Miss(jet_pt_truth_missed);
      response_mult.Miss(jet_pt_truth_missed);
      response_dnn.Miss(jet_pt_truth_missed);
      response_snn.Miss(jet_pt_truth_missed);
  }

  ///////////////////////////
  // fake file //
  TFile *fake = new TFile(fake_file, "READ");
  cout << "Fake file: " << fake_file << endl;
  TTree *faketree = (TTree*)fake->Get("tree");
  Float_t jet_pt_area_fake, jet_pt_mult_fake, jet_pt_dnn_fake, jet_pt_snn_fake;
  Float_t jet_track_pt_0_fake;
  faketree->SetBranchAddress(target_variables[0], &jet_pt_area_fake);
  faketree->SetBranchAddress(target_variables[1], &jet_pt_mult_fake);
  faketree->SetBranchAddress(target_variables[2], &jet_pt_dnn_fake);
  faketree->SetBranchAddress(target_variables[3], &jet_pt_snn_fake);
  faketree->SetBranchAddress("jet_track_pt_0", &jet_track_pt_0_fake);

  for (Int_t i = 0; i < faketree->GetEntries(); i++)
  {
      faketree->GetEntry(i);
      // if(jet_track_pt_0_fake < 7.0) continue; // cut on leading track pT
      // fill fake response
      response_area.Fake(jet_pt_area_fake);
      response_mult.Fake(jet_pt_mult_fake);
      response_dnn.Fake(jet_pt_dnn_fake);
      response_snn.Fake(jet_pt_snn_fake);

  }

  ///////////////////////////
  // unmatched file //
  TFile *unmatched= new TFile(unmatched_file, "READ");
  cout << "Unmatched file: " << unmatched_file << endl;
  TTree *unmatchedtree = (TTree*)unmatched->Get("tree");
  Float_t jet_pt_area_unmatched, jet_pt_mult_unmatched, jet_pt_dnn_unmatched, jet_pt_snn_unmatched;
  Float_t jet_track_pt_0_unmatched;
  unmatchedtree->SetBranchAddress(target_variables[0], &jet_pt_area_unmatched);
  unmatchedtree->SetBranchAddress(target_variables[1], &jet_pt_mult_unmatched);
  unmatchedtree->SetBranchAddress(target_variables[2], &jet_pt_dnn_unmatched);
  unmatchedtree->SetBranchAddress(target_variables[3], &jet_pt_snn_unmatched);
  unmatchedtree->SetBranchAddress("jet_track_pt_0", &jet_track_pt_0_unmatched);

  for ( Int_t i = 0; i < unmatchedtree->GetEntries(); i++)
  {
      unmatchedtree->GetEntry(i);
      // fill measurement histograms
      //if(jet_track_pt_0_unmatched < 7.0) continue; // cut on leading track pT
      measurement_area_reco->Fill(jet_pt_area_unmatched);
      measurement_mult_reco->Fill(jet_pt_mult_unmatched);
      measurement_dnn_reco->Fill(jet_pt_dnn_unmatched);
      measurement_snn_reco->Fill(jet_pt_snn_unmatched);

      area_reco_uncut->Fill(jet_pt_area_unmatched);
      mult_reco_uncut->Fill(jet_pt_mult_unmatched); 
      dnn_reco_uncut->Fill(jet_pt_dnn_unmatched);
      snn_reco_uncut->Fill(jet_pt_snn_unmatched);
  }

  ///////////////////////////
  // truth file //
  TFile *truth = new TFile(truth_file, "READ");
  cout << "Truth file: " << truth_file << endl;
  TTree *truthtree = (TTree*)truth->Get("tree");
  Float_t jet_pt_true;
  truthtree->SetBranchAddress("jet_pt", &jet_pt_true);

  // for (Int_t i = 0; i < truthtree->GetEntries(); i++)
  // {
  //     truthtree->GetEntry(i);
  //     // fill truth histograms
  //     h_truth->Fill(jet_pt_true);
  //     h_truth_measurement_area_bins->Fill(jet_pt_true);
  //     h_truth_measurement_mult_bins->Fill(jet_pt_true);
  //     h_truth_measurement_dnn_bins->Fill(jet_pt_true);
  //     h_truth_measurement_snn_bins->Fill(jet_pt_true);

  // }


  // -------------------------//
  // Unfolding
  // Chi Square declaration
  Double_t chi2_area[totalIter];
  Double_t chi2_mult[totalIter];
  Double_t chi2_dnn[totalIter];
  Double_t chi2_snn[totalIter];

  for(Int_t i = nIter_low; i<nIter_high+1; i++){

    cout << "Iteration: " << i << endl;
    // Unfold Histograms
    RooUnfoldBayes unfoldArea(&response_area, measurement_area_reco, i);
    RooUnfoldBayes unfoldMult(&response_mult, measurement_mult_reco, i);
    RooUnfoldBayes unfoldDNN(&response_dnn, measurement_dnn_reco, i);
    RooUnfoldBayes unfoldSNN(&response_snn, measurement_snn_reco, i);

    // Get Unfolded Histograms
    unfolded_area[i-nIter_low] = (TH1D*)unfoldArea.Hunfold(errorTreatment);
    unfolded_mult[i-nIter_low] = (TH1D*)unfoldMult.Hunfold(errorTreatment);
    unfolded_dnn[i-nIter_low] = (TH1D*)unfoldDNN.Hunfold(errorTreatment);
    unfolded_snn[i-nIter_low] = (TH1D*)unfoldSNN.Hunfold(errorTreatment);
    unfolded_area[i-nIter_low]->SetName(Form("unfolded_area_Iteration%d",i));
    unfolded_mult[i-nIter_low]->SetName(Form("unfolded_mult_Iteration%d",i));
    unfolded_dnn[i-nIter_low]->SetName(Form("unfolded_dnn_Iteration%d",i));
    unfolded_snn[i-nIter_low]->SetName(Form("unfolded_snn_Iteration%d",i));

    //Chi2
    chi2_area[i-nIter_low] = GetChi2(unfolded_area[i-nIter_low], h_truth);
    chi2_mult[i-nIter_low] = GetChi2(unfolded_mult[i-nIter_low], h_truth);
    chi2_dnn[i-nIter_low] = GetChi2(unfolded_dnn[i-nIter_low], h_truth);
    chi2_snn[i-nIter_low] = GetChi2(unfolded_snn[i-nIter_low], h_truth);
    
    // Get Unfolded Ratio Histograms
    area_unfolded_over_truth[i-nIter_low] = (TH1D*)unfolded_area[i-nIter_low]->Clone(Form("area_unfolded_over_truth_Iteration%d",i));
    mult_unfolded_over_truth[i-nIter_low] = (TH1D*)unfolded_mult[i-nIter_low]->Clone(Form("mult_unfolded_over_truth_Iteration%d",i));
    dnn_unfolded_over_truth[i-nIter_low] = (TH1D*)unfolded_dnn[i-nIter_low]->Clone(Form("dnn_unfolded_over_truth_Iteration%d",i));
    snn_unfolded_over_truth[i-nIter_low] = (TH1D*)unfolded_snn[i-nIter_low]->Clone(Form("snn_unfolded_over_truth_Iteration%d",i));
    area_unfolded_over_truth[i-nIter_low]->Divide(h_truth);
    mult_unfolded_over_truth[i-nIter_low]->Divide(h_truth);
    dnn_unfolded_over_truth[i-nIter_low]->Divide(h_truth);
    snn_unfolded_over_truth[i-nIter_low]->Divide(h_truth);

    // Get Pearson Correlation Coefficients
    TMatrixD covmatMult = unfoldMult.Eunfold(errorTreatment);
    TMatrixD covmatArea = unfoldArea.Eunfold(errorTreatment);
    TMatrixD covmatDNN = unfoldDNN.Eunfold(errorTreatment);
    TMatrixD covmatSNN = unfoldSNN.Eunfold(errorTreatment);

    area_covariance[i-nIter_low] = new TH2D(covmatArea);
    mult_covariance[i-nIter_low] = new TH2D(covmatMult);
    dnn_covariance[i-nIter_low] = new TH2D(covmatDNN);
    snn_covariance[i-nIter_low] = new TH2D(covmatSNN);

    area_covariance[i-nIter_low]->SetName(Form("area_covariance_Iteration%d",i));
    mult_covariance[i-nIter_low]->SetName(Form("mult_covariance_Iteration%d",i));
    dnn_covariance[i-nIter_low]->SetName(Form("dnn_covariance_Iteration%d",i));
    snn_covariance[i-nIter_low]->SetName(Form("snn_covariance_Iteration%d",i));

    TMatrixD MpearsonArea = CalculatePearsonCoefficients(covmatArea);
    TMatrixD MpearsonMult = CalculatePearsonCoefficients(covmatMult);
    TMatrixD MpearsonSNN = CalculatePearsonCoefficients(covmatSNN);
    TMatrixD MpearsonDNN = CalculatePearsonCoefficients(covmatDNN);

    area_pearson[i-nIter_low] = new TH2D(MpearsonArea);
    mult_pearson[i-nIter_low] = new TH2D(MpearsonMult);
    dnn_pearson[i-nIter_low] = new TH2D(MpearsonDNN);
    snn_pearson[i-nIter_low] = new TH2D(MpearsonSNN);
    area_pearson[i-nIter_low]->SetName(Form("area_pearson_Iteration%d",i));
    mult_pearson[i-nIter_low]->SetName(Form("mult_pearson_Iteration%d",i));
    dnn_pearson[i-nIter_low]->SetName(Form("dnn_pearson_Iteration%d",i));
    snn_pearson[i-nIter_low]->SetName(Form("snn_pearson_Iteration%d",i));
  
}
 
  // Calculate ratios of measured to truth
  area_measurement_over_truth = (TH1D*)measurement_area_reco->Clone("area_measurement_over_truth");
  mult_measurement_over_truth = (TH1D*)measurement_mult_reco->Clone("mult_measurement_over_truth");
  dnn_measurement_over_truth = (TH1D*)measurement_dnn_reco->Clone("dnn_measurement_over_truth");
  snn_measurement_over_truth = (TH1D*)measurement_snn_reco->Clone("snn_measurement_over_truth");
  area_measurement_over_truth->Divide(h_truth_measurement_area_bins);
  mult_measurement_over_truth->Divide(h_truth_measurement_mult_bins);
  dnn_measurement_over_truth->Divide(h_truth_measurement_dnn_bins);
  snn_measurement_over_truth->Divide(h_truth_measurement_snn_bins);
 
  // area_measurement_over_truth = (TH1D*)area_reco_uncut->Clone("area_measurement_over_truth");
  // mult_measurement_over_truth = (TH1D*)mult_reco_uncut->Clone("mult_measurement_over_truth");
  // dnn_measurement_over_truth = (TH1D*)dnn_reco_uncut->Clone("dnn_measurement_over_truth");
  // snn_measurement_over_truth = (TH1D*)snn_reco_uncut->Clone("snn_measurement_over_truth");
  // area_measurement_over_truth->Divide(h_truth);
  // mult_measurement_over_truth->Divide(h_truth);
  // dnn_measurement_over_truth->Divide(h_truth);
  // snn_measurement_over_truth->Divide(h_truth);
  
  area_measurement_over_truth_normalized = (TH1D*)area_reco_uncut->Clone("area_measurement_over_truth_normalized");
  mult_measurement_over_truth_normalized = (TH1D*)mult_reco_uncut->Clone("mult_measurement_over_truth_normalized");
  dnn_measurement_over_truth_normalized = (TH1D*)dnn_reco_uncut->Clone("dnn_measurement_over_truth_normalized");
  snn_measurement_over_truth_normalized = (TH1D*)snn_reco_uncut->Clone("snn_measurement_over_truth_normalized");
  TH1D *h_truth_normalized = (TH1D*)h_truth->Clone("h_truth_normalized");
  h_truth_normalized->Scale(1./h_truth_normalized->Integral());
  area_measurement_over_truth_normalized->Scale(1./area_measurement_over_truth_normalized->Integral());
  mult_measurement_over_truth_normalized->Scale(1./mult_measurement_over_truth_normalized->Integral());
  dnn_measurement_over_truth_normalized->Scale(1./dnn_measurement_over_truth_normalized->Integral());
  snn_measurement_over_truth_normalized->Scale(1./snn_measurement_over_truth_normalized->Integral());
  area_measurement_over_truth_normalized->Divide(h_truth_normalized);
  mult_measurement_over_truth_normalized->Divide(h_truth_normalized);
  dnn_measurement_over_truth_normalized->Divide(h_truth_normalized);
  snn_measurement_over_truth_normalized->Divide(h_truth_normalized);






  // Write to file
  cout << "Writing to file" << endl;
  TFile *f = new TFile(output_file.Data(),"RECREATE");
  f->cd();
  TTree *t = new TTree("tree","tree");
  Double_t chimult, chiarea, chidnn, chisnn;
  t->Branch("chimult",&chimult);
  t->Branch("chiarea",&chiarea);
  t->Branch("chidnn",&chidnn);
  t->Branch("chisnn",&chisnn);
  // fill tree
  for(int i = nIter_low; i < nIter_high; i++){
    chimult = chi2_area[i-nIter_low];
    chiarea = chi2_mult[i-nIter_low];
    chidnn = chi2_dnn[i-nIter_low];
    chisnn = chi2_snn[i-nIter_low];
    t->Fill();
  }
  t->Write();

  // write histograms to file
  h_truth->Write();
  measurement_area_reco->Write();
  measurement_mult_reco->Write();
  measurement_dnn_reco->Write();
  measurement_snn_reco->Write();
  area_measurement_over_truth->Write();
  mult_measurement_over_truth->Write();
  dnn_measurement_over_truth->Write();
  snn_measurement_over_truth->Write();  
  area_measurement_over_truth_normalized->Write();
  mult_measurement_over_truth_normalized->Write();
  dnn_measurement_over_truth_normalized->Write();
  snn_measurement_over_truth_normalized->Write();

  for (int i = nIter_low; i < nIter_high; i++){
    area_pearson[i-nIter_low]->Write();
    mult_pearson[i-nIter_low]->Write();
    dnn_pearson[i-nIter_low]->Write();
    snn_pearson[i-nIter_low]->Write();
    unfolded_area[i-nIter_low]->Write();
    unfolded_mult[i-nIter_low]->Write();
    unfolded_dnn[i-nIter_low]->Write();
    unfolded_snn[i-nIter_low]->Write();
    area_unfolded_over_truth[i-nIter_low]->Write();
    mult_unfolded_over_truth[i-nIter_low]->Write();
    dnn_unfolded_over_truth[i-nIter_low]->Write();
    snn_unfolded_over_truth[i-nIter_low]->Write();
    area_covariance[i-nIter_low]->Write();
    mult_covariance[i-nIter_low]->Write();
    dnn_covariance[i-nIter_low]->Write();
    snn_covariance[i-nIter_low]->Write();
  }


  f->Close();
  cout << "Done" << endl;
  matched->Close();
  unmatched->Close();
  fake->Close();
  missed->Close();
  truth->Close();

  return 0;

}

#ifndef __CINT__
Int_t main (Int_t argc, char **argv){
  if(argc < 7){
    cout << "Usage: ./Unfold <matched file> <unmatched file> <missed file> <fake file> <truth file> <output file>" << endl;
    return 0;
  }
  TString matchedfile = argv[1];
  TString unmatchedfile = argv[2];
  TString missedfile = argv[3];
  TString fakefile = argv[4];
  TString truthfile = argv[5];
  TString outputfile = argv[6];

  cout << "Matched file: " << matchedfile << endl;
  cout << "Unmatched file: " << unmatchedfile << endl;
  cout << "Missed file: " << missedfile << endl;
  cout << "Fake file: " << fakefile << endl;
  cout << "Truth file: " << truthfile << endl;
  cout << "Output file: " << outputfile << endl;

  // return 0;
  return Unfold(matchedfile, unmatchedfile, missedfile, fakefile, truthfile, outputfile);
 
}
#endif
