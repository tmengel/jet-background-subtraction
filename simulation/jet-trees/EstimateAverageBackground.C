#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"


#include <iostream>
#include <string>
#include <vector>

using namespace std;


Int_t main(Int_t argc, Char_t** argv){

    TString input = argv[1];
    TString output = argv[2];

    TFile tenngen_file(input.Data());
    if(!tenngen_file.IsOpen()){
        cout << "Error: could not open tenngen file" << endl;
        return 1;
    }

    TTreeReader tenngen_parts("tree", &tenngen_file);
    TTreeReaderValue<Int_t> tenngen_numparts(tenngen_parts, "nparts");
    TTreeReaderArray<Float_t> tenngen_px(tenngen_parts, "particle_px");
    TTreeReaderArray<Float_t> tenngen_py(tenngen_parts, "particle_py");
    TTreeReaderArray<Float_t> tenngen_pz(tenngen_parts, "particle_pz");
    TTreeReaderArray<Float_t> tenngen_E(tenngen_parts, "particle_e");

    Float_t avgpT;
    TFile *output_file = new TFile(output.Data(), "RECREATE");
    TTree *output_tree = new TTree("tree", "tree");
    output_tree->Branch("avgpT", &avgpT, "avgpT/F");


    while(tenngen_parts.Next()){
        Float_t pT_sum = 0;
        Int_t num_particles =0;
        for(Int_t i = 0; i < *tenngen_numparts; i++){
            Float_t pTtmp = TMath::Sqrt(TMath::Power(tenngen_px[i], 2) + TMath::Power(tenngen_py[i], 2));
            if(pTtmp > 0.15){
                pT_sum += pTtmp;
                num_particles++;
            }
        }
        avgpT = pT_sum / num_particles;
        output_tree->Fill();
    }

    output_file->Write();
    output_file->Close();
    tenngen_file.Close();

    return 0;


}
