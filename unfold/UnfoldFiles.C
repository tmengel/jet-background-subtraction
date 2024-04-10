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

void UnfoldFiles(TString inputfile, TString responsefile, TString outputfile){
    
    cout << "Unfolding " << inputfile <<  " with response " << responsefile << endl;
    cout << "Output file " << outputfile << endl;



    TFile * fout = new TFile(outputfile.Data(), "RECREATE");
    TH1D * hMeasured[n_pt_bins];
    TH1D * hTruth[n_pt_bins];
    TH1D * hUnfolded[n_pt_bins];


    TFile * fin = new TFile(inputfile.Data(), "READ");
    TFile * fresponse = new TFile(responsefile.Data(), "READ");

    // set up error treatment
    RooUnfold::ErrorTreatment errorTreatment = RooUnfold::kCovariance;
    RooUnfoldResponse responses[n_pt_bins];

    for(Int_t ipt = 1; ipt < n_pt_bins; ipt++){
        cout << "\tUnfolding pt bin " << ipt << endl;
        TH1D * hMeas = (TH1D*)fin->Get(Form("hMeas_%d", ipt));
        TH1D * hTru = (TH1D*)fin->Get(Form("hTruth_%d", ipt));
        TH2D * hResponse = (TH2D*)fresponse->Get(Form("hResponse_%d", ipt));

        responses[ipt] = RooUnfoldResponse(hMeas, hTru, hResponse,Form("response_%d", ipt), Form("response_%d", ipt));
        RooUnfoldBayes unfold(&responses[ipt], hMeas, 4);
        TH1D * hUnfold = (TH1D*)unfold.Hunfold(errorTreatment);
        
        // save histos to array
        hMeasured[ipt] = (TH1D*)hMeas->Clone(Form("hMeasured_%d", ipt));
        hTruth[ipt] = (TH1D*)hTru->Clone(Form("hTruth_%d", ipt));
        hUnfolded[ipt] = (TH1D*)hUnfold->Clone(Form("hUnfolded_%d", ipt));
        hMeasured[ipt]->SetDirectory(0);
        hTruth[ipt]->SetDirectory(0);
        hUnfolded[ipt]->SetDirectory(0);

        // delete histos
        delete hMeas;
        delete hTru;
        delete hUnfold;
        delete hResponse;
    }

    // save histos to file
    fout->cd();
    for(Int_t ipt = 1; ipt < n_pt_bins; ipt++){
        hMeasured[ipt]->Write();
        hTruth[ipt]->Write();
        hUnfolded[ipt]->Write();
    }
    fout->Write();
    fout->Close();
    cout << "Done" << endl;
    fin->Close();
    fresponse->Close();
    cout << "Saved" << endl;

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
    UnfoldFiles(inputfile, responsefile, outputfile);
    return 0;


}
#endif