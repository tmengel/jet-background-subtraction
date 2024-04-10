#if !(defined(__CINT__) || defined(__CLING__)) || defined(__ACLIC__)

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

#include <iostream>
#include <vector>

#endif

const Int_t n_pt_bins = 12;
const Double_t pt_bins[13] = {-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50, 60, 70, 80, 90, 100, 110};
std::vector<Double_t> pt_bins_vec = {-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 50, 60, 70, 80, 90, 100, 110};
  
Double_t GetChi2(TH1D *Unfolded, TH1D *Truth)
{
  Double_t chi2 = 0;
  for(Int_t i = 1; i <= Unfolded->GetNbinsX(); i++)
  {
    chi2 += TMath::Power((Truth->GetBinContent(i)-Unfolded->GetBinContent(i))/Unfolded->GetBinError(i),2.);
  }
  return chi2;
}

void UnfoldCorrelations(TString inputfile){
    
    cout << "Starting " << inputfile << endl;
    // get path to file
    TString inputdirectory = inputfile;
    inputdirectory.Remove(inputdirectory.Last('/'), inputdirectory.Length() - inputdirectory.Last('/'));
    TString outputfile = Form("%s/Unfolded.root", inputdirectory.Data());
    cout << "Input directory: " << inputdirectory << endl;
    cout << "Output file: " << outputfile << endl;
    TFile * fout = new TFile(outputfile.Data(), "RECREATE");
    TFile * fin = new TFile(inputfile.Data(), "READ");
    std::vector<TString> rotation_options = {"A", "B", "C", "cos","linear","log"};
    const Int_t n_rotation_options = rotation_options.size();
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;

    for(Int_t irotation = 0; irotation < rotation_options.size(); irotation++){
        cout << "Unfolding " << rotation_options[irotation] << endl;
        TH2D * h2_diff_measured = (TH2D*)fin->Get(Form("h2_diff_measured_%s", rotation_options[irotation].Data()));
        TH2D * h2_diff_truth = (TH2D*)fin->Get(Form("h2_diff_truth_%s", rotation_options[irotation].Data()));
        Double_t xmin = h2_diff_measured->GetXaxis()->GetXmin();
        Double_t xmax = h2_diff_measured->GetXaxis()->GetXmax();
        Double_t ymin = h2_diff_measured->GetYaxis()->GetXmin();
        Double_t ymax = h2_diff_measured->GetYaxis()->GetXmax();
        Int_t nxbins = h2_diff_measured->GetNbinsX();
        Int_t nybins = h2_diff_measured->GetNbinsY();
        TH2D * h2_diff_unfolded[n_rotation_options];
        for (Int_t i = 0; i < n_rotation_options; i++){
            h2_diff_unfolded[i] = new TH2D(Form("h2_diff_unfolded_%s_%s", rotation_options[irotation].Data(), rotation_options[i].Data()), "", nxbins, xmin, xmax, nybins, ymin, ymax);
        }
        
        for(Int_t ipt = 1; ipt < n_pt_bins; ipt++){
            cout << "\tUnfolding pt bin " << ipt << endl;
            TH1D * hMeas = (TH1D*)fin->Get(Form("hMeas_%d_%s", ipt, rotation_options[irotation].Data()));
            TH1D * hTruth = (TH1D*)fin->Get(Form("hTruth_%d_%s", ipt, rotation_options[irotation].Data()));
            for(Int_t jrotation = 0; jrotation < n_rotation_options; jrotation++){
                cout << "\t\tUnfolding rotation option " << rotation_options[jrotation] << endl;
                TH2D * hResponse = (TH2D*)fin->Get(Form("hResponse_%d_%s", ipt, rotation_options[jrotation].Data()));
                RooUnfoldResponse response(hMeas, hTruth, hResponse, "response", "response");
                RooUnfoldBayes unfold(&response, hMeas, 9);
                TH1D * hUnfolded = (TH1D*)unfold.Hunfold(errorTreatment);
                hUnfolded->SetName(Form("hUnfolded_%d_%s", ipt, rotation_options[jrotation].Data()));
                hUnfolded->SetTitle(Form("hUnfolded_%d_%s", ipt, rotation_options[jrotation].Data()));
                hUnfolded->SetDirectory(0);
                for(Int_t iy = 1; iy <= nybins; iy++){
                    h2_diff_unfolded[jrotation]->SetBinContent(ipt+1, iy, hUnfolded->GetBinContent(iy));
                    h2_diff_unfolded[jrotation]->SetBinError(ipt+1, iy, hUnfolded->GetBinError(iy));
                }
                delete hUnfolded;
            }
            delete hMeas;
            delete hTruth;
        }
        fout->cd();
        for(Int_t i = 0; i < n_rotation_options; i++){
            h2_diff_unfolded[i]->Write();
            h2_diff_measured->Write();
            h2_diff_truth->Write();
        }
        
        delete h2_diff_measured;
        delete h2_diff_truth;
        for(Int_t i = 0; i < n_rotation_options; i++){
            delete h2_diff_unfolded[i];
        }
        cout << "Finished " << rotation_options[irotation] << endl;
    }
    cout << "Finished " << inputfile << endl;
    fout->cd();
    fout->Write();
    fout->Close();
    fin->Close();
    cout << "Done" << endl;

}

#ifndef __CINT__
Int_t main(Int_t argc, Char_t** argv){

    cout << "Starting " << argv[0] << endl;
    // get arguments
    // check if correct number of arguments
    if(argc != 2){
        cout << "Usage: " << argv[0] << " <inputfile> " << endl;
        cout << "Incorrect number of arguments" << endl; exit(-1);
    }
   // get arguments
    TString inputfile = argv[1];
    UnfoldCorrelations(inputfile);
    return 0;


}
#endif