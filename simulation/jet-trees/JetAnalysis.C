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

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/config.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace fastjet;

struct Args{
    Int_t collen;
    Int_t ptbin;
    Float_t jetparam;
    TString mode;
    Int_t nEvents;
    Int_t verbose;
};

struct Config{
    Int_t n_cent_bins;
    TString collsys;
    std::vector<TString> cent_dir_names;
    std::vector<Float_t> cent_event_start_fractions;
    std::vector<Float_t> cent_event_end_fractions;
    TString output_base_directory;
    TString pythia_base_directory;
    TString tenngen_base_directory;
    TString pythia_filename;
    TString pythia_jets_output_filename;
    std::vector<TString> output_filenames;
    std::vector<TString> tenngen_filenames;
    TString merged_output_filename;
};

struct JetFindingParameters{
    Int_t total_pythia_events;
    std::vector<Int_t> n_mixed_events;
    std::vector<Int_t> start_mixed_events;
    Int_t n_pythia_events;
    Float_t particle_max_eta;
    Float_t jet_pt_min;
    Float_t jetparam;
};

void PythiaJets(const TString input_file, const TString output_file, const struct JetFindingParameters jetparams){
    
    cout << "Starting PythiaJets()" << endl;
    cout << "Input file: " << input_file.Data() << endl;
    cout << "Output file: " << output_file.Data() << endl;
    cout << "jetparam: " << jetparams.jetparam << endl;
    cout << "nEvents: " << jetparams.n_pythia_events << endl;
    cout << "particle_max_eta: " << jetparams.particle_max_eta << endl;
    cout << "jet_pt_min: " << jetparams.jet_pt_min << endl;
    //////////////////////////////////////////////////////////
    // open input file
    TFile pythia_file(input_file.Data());
    if(!pythia_file.IsOpen()){ cout << "PYTHIA file not open" << endl; exit(1); }

    TString partTree = "tree";
    TString eventTree = "eventInfo";

    TTreeReader pythia_parts(partTree.Data(), &pythia_file);
    TTreeReader pythia_event(eventTree.Data(), &pythia_file);

    TTreeReaderValue<Int_t> pythia_nevent(pythia_event, "nevent");
    TTreeReaderValue<Int_t> pythia_ptbin(pythia_event, "ptbinID");
    TTreeReaderValue<Float_t> pythia_xsec(pythia_event, "xsec_over_eventweight");

    TTreeReaderValue<Float_t> pythia_weight(pythia_parts, "weight");
    TTreeReaderValue<Int_t> pythia_numparts(pythia_parts, "nparts");
    TTreeReaderArray<Float_t> pythia_px(pythia_parts, "particle_px");
    TTreeReaderArray<Float_t> pythia_py(pythia_parts, "particle_py");
    TTreeReaderArray<Float_t> pythia_pz(pythia_parts, "particle_pz");
    TTreeReaderArray<Float_t> pythia_E(pythia_parts, "particle_e");
    //////////////////////////////////////////////////////////

    const Int_t MAX_JETS = 100;
    Int_t event_pt_hard_bin;
    Int_t number_of_jets;
    Float_t event_weight;
    Float_t jet_pt[MAX_JETS], jet_eta[MAX_JETS], jet_phi[MAX_JETS];    
    ////////////////////////////////////////////////////////////////////////////////
    // Output file
    TFile *fout = new TFile(output_file.Data() ,"RECREATE");
    cout << "Output Tfile created: " <<  output_file.Data() << endl;
    TTree* TruthTree = new TTree("tree","tree");
    TruthTree->Branch("event_pt_hard_bin", &event_pt_hard_bin, "event_pt_hard_bin/I");
    TruthTree->Branch("event_weight", &event_weight, "event_weight/F");
    TruthTree->Branch("number_of_jets", &number_of_jets, "number_of_jets/I");
    TruthTree->Branch("jet_pt", jet_pt, "jet_pt[number_of_jets]/F");
    TruthTree->Branch("jet_eta", jet_eta, "jet_eta[number_of_jets]/F");
    TruthTree->Branch("jet_phi", jet_phi, "jet_phi[number_of_jets]/F");

    ///////////////////////////////////////////////////////
    // loop over events
    // Get pythia event info
    
    if(!pythia_event.Next()){ cout << "Reached end of pythia file" << endl; exit(1); }

    event_pt_hard_bin = *pythia_ptbin;
    Float_t pythia_xsec_over_eventweight = *pythia_xsec;
    // Initialize event counters
    Int_t event_number = 0;

    std::vector<fastjet::PseudoJet> particles;
    std::vector<fastjet::PseudoJet> jets;
    Float_t antikt_jet_abs_eta_max = jetparams.particle_max_eta - jetparams.jetparam;
    fastjet::GhostedAreaSpec ghost_area_spec(jetparams.particle_max_eta);
    fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparams.jetparam);
            
    ////////////////////////////////////////////////////////////////////////////////

    cout << "Starting pythia event loop" << endl;
    Int_t nJets_Total = 0;
    Int_t nEvents_processed = 0;
    for(Int_t iEvent = 0; iEvent < jetparams.n_pythia_events; iEvent++){

            if(!pythia_parts.Next()){ cout << "Reached end of pythia file" << endl; break;}
            
            // clear vectors
            particles.clear();
            jets.clear();
            nEvents_processed++;
            ////////////////////////////////////////////////////////////////////////////////
            // Fill the tree
            Int_t nParts_pythia = *pythia_numparts;
            // Get event weight
            event_weight = *pythia_weight;
            event_weight *= pythia_xsec_over_eventweight;
            // Get particle info
            Int_t nParts_accepeted = 0;

            for(Int_t ipart = 0; ipart < nParts_pythia; ipart++){   
                PseudoJet particle_temp(pythia_px[ipart], pythia_py[ipart], pythia_pz[ipart], pythia_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                nParts_accepeted++;
            }


            // find antikt jets
            fastjet::ClusterSequence antikt_cs(particles, antikt_jet_def);
            jets = antikt_cs.inclusive_jets(jetparams.jet_pt_min);
            Int_t number_of_jets_all = jets.size();

            // loop over jets
            number_of_jets = 0;
            for(Int_t ijet = 0; ijet < number_of_jets_all; ijet++){

                if(TMath::Abs(jets[ijet].eta()) > antikt_jet_abs_eta_max) continue;
                if(jets[ijet].pt() < jetparams.jet_pt_min) continue;
                if(ijet > MAX_JETS) continue;

                jet_pt[number_of_jets] = jets[ijet].pt();
                jet_eta[number_of_jets] = jets[ijet].eta();
                jet_phi[number_of_jets] = jets[ijet].phi();

                number_of_jets++;
                nJets_Total++;
            }

            if(number_of_jets == 0) continue;

            TruthTree->Fill();
            event_number++;

            for(Int_t ijet = 0; ijet < number_of_jets; ijet++){
                jet_pt[ijet] = 0;
                jet_eta[ijet] = 0;
                jet_phi[ijet] = 0;
            }
            
            jets.clear();
            particles.clear();
 
    }
    
    //end of event/cent loop
    ////////////////////////////////////////////////////////////////////////////////
    fout->Write();
    fout->Close();
    pythia_file.Close();

    cout << "Number of events: " << nEvents_processed << endl;
    cout << "Number of jets: " << nJets_Total << endl;
}

void Jets(const TString input_file_pythia, const TString input_file_tenngen, const TString output_file, const Int_t centbin, const struct JetFindingParameters jetparams){ 
    
    cout << "Starting Jets()" << endl;
    cout << "Input pythia file: " << input_file_pythia << endl;
    cout << "Input TennGen file: " << input_file_tenngen << endl;
    cout << "Output file: " << output_file << endl;
    cout << "Jet param: " << jetparams.jetparam << endl;
    cout << "Starting event: " << jetparams.start_mixed_events.at(centbin) << endl;
    cout << "Total events: " << jetparams.n_mixed_events.at(centbin) << endl;


    // open input file
    TFile pythia_file(input_file_pythia.Data());
    TFile tenngen_file(input_file_tenngen.Data());
    if(!pythia_file.IsOpen()){  cout << "PYTHIA file not open" << endl; exit(1); }
    if(!tenngen_file.IsOpen()){  cout << "TennGen file not open" << endl; exit(1); }
    
    TString partTree = "tree";
    TString eventTree = "eventInfo";
    TTreeReader pythia_parts(partTree.Data(), &pythia_file);
    TTreeReader pythia_event(eventTree.Data(), &pythia_file);
    TTreeReaderValue<Int_t> pythia_nevent(pythia_event, "nevent");
    TTreeReaderValue<Float_t> pythia_ptmax(pythia_event, "ptmin");
    TTreeReaderValue<Float_t> pythia_ptmin(pythia_event, "ptmax");
    TTreeReaderValue<Int_t> pythia_ptbin(pythia_event, "ptbinID");
    TTreeReaderValue<Float_t> pythia_xsec(pythia_event, "xsec_over_eventweight");
    TTreeReaderValue<Float_t> pythia_weight(pythia_parts, "weight");
    TTreeReaderValue<Int_t> pythia_numparts(pythia_parts, "nparts");
    TTreeReaderArray<Float_t> pythia_px(pythia_parts, "particle_px");
    TTreeReaderArray<Float_t> pythia_py(pythia_parts, "particle_py");
    TTreeReaderArray<Float_t> pythia_pz(pythia_parts, "particle_pz");
    TTreeReaderArray<Float_t> pythia_E(pythia_parts, "particle_e");
    //////////////////////////////////////////////////////////////
    TTreeReader tenngen_parts(partTree.Data(), &tenngen_file);
    TTreeReaderValue<Int_t> tenngen_numparts(tenngen_parts, "nparts");
    TTreeReaderArray<Float_t> tenngen_px(tenngen_parts, "particle_px");
    TTreeReaderArray<Float_t> tenngen_py(tenngen_parts, "particle_py");
    TTreeReaderArray<Float_t> tenngen_pz(tenngen_parts, "particle_pz");
    TTreeReaderArray<Float_t> tenngen_E(tenngen_parts, "particle_e");
    // check if the trees are empty
    if(!pythia_event.Next()){  cout << "PYTHIA tree is empty" << endl; exit(1); }
    if(!tenngen_parts.Next()){ cout << "TennGen tree is empty" << endl; exit(1); }
    //////////////////////// Create output file //////////////////////////

    const Int_t MAX_TRACKS= 250;
    Int_t event_cent_bin;
    Int_t event_pt_hard_bin;
    Float_t median_pt_over_npart, median_pt_over_area;
    Float_t event_weight;
    Float_t jet_pt_raw;
    Float_t jet_eta, jet_phi, jet_area;
    Float_t jet_pt_pythia;
    // Float_t jet_track_pt[MAX_TRACKS], jet_track_eta[MAX_TRACKS], jet_track_phi[MAX_TRACKS];
    Int_t jet_nparts_pythia, jet_nparts;
    Float_t jet_track_pt_0, jet_track_pt_1, jet_track_pt_2, jet_track_pt_3, jet_track_pt_4, jet_track_pt_5, jet_track_pt_6, jet_track_pt_7;
    Float_t jet_angularity;

    ////////////////////////////////////////////////////////////////////////////////
    // Output file
    TFile *fout = new TFile(output_file.Data(), "RECREATE");
    TTree *outTree = new TTree("tree", "tree");
    outTree->Branch("event_cent_bin", &event_cent_bin, "event_cent_bin/I");
    outTree->Branch("event_pt_hard_bin", &event_pt_hard_bin, "event_pt_hard_bin/I");
    outTree->Branch("median_pt_over_npart", &median_pt_over_npart, "median_pt_over_npart/F");
    outTree->Branch("median_pt_over_area", &median_pt_over_area, "median_pt_over_area/F");
    outTree->Branch("event_weight", &event_weight, "event_weight/F");
    outTree->Branch("jet_pt_raw", &jet_pt_raw, "jet_pt_raw/F");
    outTree->Branch("jet_pt_pythia", &jet_pt_pythia, "jet_pt_pythia/F");
    // outTree->Branch("jet_eta", &jet_eta, "jet_eta/F");
    // outTree->Branch("jet_phi", &jet_phi, "jet_phi/F");
    outTree->Branch("jet_area", &jet_area, "jet_area/F");
    outTree->Branch("jet_nparts", &jet_nparts, "jet_nparts/I");
    outTree->Branch("jet_nparts_pythia", &jet_nparts_pythia, "jet_nparts_pythia/I");
    outTree->Branch("jet_angularity", &jet_angularity, "jet_angularity/F");
    outTree->Branch("jet_track_pt_0", &jet_track_pt_0, "jet_track_pt_0/F");
    outTree->Branch("jet_track_pt_1", &jet_track_pt_1, "jet_track_pt_1/F");
    outTree->Branch("jet_track_pt_2", &jet_track_pt_2, "jet_track_pt_2/F");
    outTree->Branch("jet_track_pt_3", &jet_track_pt_3, "jet_track_pt_3/F");
    outTree->Branch("jet_track_pt_4", &jet_track_pt_4, "jet_track_pt_4/F");
    outTree->Branch("jet_track_pt_5", &jet_track_pt_5, "jet_track_pt_5/F");
    outTree->Branch("jet_track_pt_6", &jet_track_pt_6, "jet_track_pt_6/F");
    outTree->Branch("jet_track_pt_7", &jet_track_pt_7, "jet_track_pt_7/F");
    // outTree->Branch("jet_track_pt", jet_track_pt, "jet_track_pt[jet_nparts]/F");
    // outTree->Branch("jet_track_eta", jet_track_eta, "jet_track_eta[jet_nparts]/F");
    // outTree->Branch("jet_track_phi", jet_track_phi, "jet_track_phi[jet_nparts]/F");
    
    ////////////////////////////////////////////////////////////////////////////////
    // loop over events
    // Get pythia event info

    event_cent_bin = centbin;
    event_pt_hard_bin = *pythia_ptbin;
    Float_t pythia_xsec_over_eventweight = *pythia_xsec;
    Int_t pythia_nevents = *pythia_nevent;

    // Initialize event counters
    
    Int_t event_number = 0;

    std::vector<fastjet::PseudoJet> particles;
    std::vector<fastjet::PseudoJet> jets;
    std::vector<fastjet::PseudoJet> constituents;
    std::vector<Float_t> median_pt_over_npart_calc_vec;
    std::vector<Float_t> median_pt_over_area_calc_vec;
    Float_t antikt_jet_abs_eta_max = jetparams.particle_max_eta-jetparams.jetparam;
    fastjet::GhostedAreaSpec ghost_area_spec(jetparams.particle_max_eta);
    fastjet::JetDefinition kt_jet_def(fastjet::kt_algorithm, 0.4);

    // skip pythia events before start_event
    for (Int_t ievent = 0; ievent < jetparams.start_mixed_events.at(centbin); ievent++)
    { 
        if(!pythia_parts.Next()){ cout << "Reached end of pythia file at event " << ievent << endl; exit(1);}
    }
    Int_t nEvents_processed = 0;
    Int_t nJets_total = 0;
    for(Int_t ievent = 0; ievent < jetparams.n_mixed_events.at(centbin); ievent++){
            if(!tenngen_parts.Next()){ cout << "Reached end of tenngen file" << endl; break;}
            if(!pythia_parts.Next()){ cout << "Reached end of pythia file" << endl; break;}
          
            //clear vectors
            particles.clear();
            jets.clear();
            constituents.clear();
            median_pt_over_npart_calc_vec.clear();
            median_pt_over_area_calc_vec.clear();
            ////////////////////////////////////////////////////////////////////////////////
            // Fill the tree
            Int_t nParts_pythia = *pythia_numparts;
            Int_t nParts_tenngen = *tenngen_numparts;
            Int_t nParts_total=0;
            // Get event weight
            event_weight = *pythia_weight;
            event_weight *= pythia_xsec_over_eventweight;
            // Get particle info
         
            for(Int_t ipart = 0; ipart < nParts_pythia; ipart++){   

                PseudoJet particle_temp(pythia_px[ipart], pythia_py[ipart], pythia_pz[ipart], pythia_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(1);
                nParts_total++;

            }
            for(Int_t ipart = 0; ipart < nParts_tenngen; ipart++){
                
                PseudoJet particle_temp(tenngen_px[ipart], tenngen_py[ipart], tenngen_pz[ipart], tenngen_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(0);
                nParts_total++;

            }
            
            // find kt jets
            fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);  
            fastjet::ClusterSequenceArea kt_cs(particles, kt_jet_def, area_def);
            fastjet::Selector jet_selector_kt = (!SelectorNHardest(2));
            jets.clear();
            jets = jet_selector_kt(kt_cs.inclusive_jets());
            // find median kt jet pt
            Int_t nktjets = jets.size();
            median_pt_over_area_calc_vec.clear();
            median_pt_over_npart_calc_vec.clear();
            for(Int_t ijet = 0; ijet < nktjets; ijet++){

                if(!jets[ijet].has_area()) continue;
                Float_t kt_jet_area = jets[ijet].area();

                if(jets[ijet].constituents().size() == 0) continue;
                Int_t kt_jet_nparts = jets[ijet].constituents().size();

                if(jets[ijet].pt() == 0) continue;
                Float_t kt_jet_pt = jets[ijet].pt();
                
                Float_t pt_over_area = 1.0*(kt_jet_pt/kt_jet_area);
                Float_t pt_over_npart = 1.0*(kt_jet_pt/kt_jet_nparts);
                median_pt_over_area_calc_vec.push_back(pt_over_area);
                median_pt_over_npart_calc_vec.push_back(pt_over_npart);

            }
            // sort vectors
            std::sort(median_pt_over_npart_calc_vec.begin(), median_pt_over_npart_calc_vec.end());
            std::sort(median_pt_over_area_calc_vec.begin(), median_pt_over_area_calc_vec.end());
            median_pt_over_npart = 0;
            median_pt_over_area = 0;
            // find median
            if(median_pt_over_npart_calc_vec.size() > 0) median_pt_over_npart = median_pt_over_npart_calc_vec[median_pt_over_npart_calc_vec.size()/2];
            if(median_pt_over_area_calc_vec.size() > 0) median_pt_over_area = median_pt_over_area_calc_vec[median_pt_over_area_calc_vec.size()/2];
            //clear vectors
            jets.clear();
            median_pt_over_npart_calc_vec.clear();
            median_pt_over_area_calc_vec.clear();
            ////////////////////////////////////////////////////////////////////////////////
            // find antikt jets
            fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparams.jetparam);
            fastjet::ClusterSequenceArea antikt_cs(particles, antikt_jet_def, area_def);
            jets = antikt_cs.inclusive_jets(jetparams.jet_pt_min);

            // loop over antikt jets
            Int_t n_antikt_jets = jets.size();
            for(Int_t ijet = 0; ijet < n_antikt_jets; ijet++){

                if(!jets[ijet].has_area()) continue;
                if(jets[ijet].constituents().size() == 0) continue;
                if(jets[ijet].pt() < jetparams.jet_pt_min) continue;
                if(TMath::Abs(jets[ijet].eta()) > antikt_jet_abs_eta_max) continue;

                constituents.clear();
                constituents = sorted_by_pt(jets[ijet].constituents());
                Int_t n_constituents = constituents.size();

                jet_pt_pythia = 0;
                jet_pt_raw = jets[ijet].pt();
                jet_area = jets[ijet].area();
                jet_nparts = 0;
                jet_nparts_pythia = 0;
                jet_eta = jets[ijet].eta();
                jet_phi = jets[ijet].phi();
                Float_t jet_track_pt[8] = {0,0,0,0,0,0,0,0};
                Float_t temp_angularity = 0;
                for(Int_t ipart =0; ipart< n_constituents; ipart++){
                    if(constituents[ipart].user_index() == 1 || constituents[ipart].user_index() == 0){
                        Float_t deta = constituents[ipart].eta() - jet_eta;
                        Float_t dphi = constituents[ipart].phi() - jet_phi;
                        if(dphi < 0) dphi = -1.0*dphi;
                        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                        Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
                        temp_angularity += (dR/jetparams.jetparam)*(constituents[ipart].pt()/jet_pt_raw);
                        if(jet_nparts < 8){
                            jet_track_pt[jet_nparts] = constituents[ipart].pt();
                        }
                        jet_nparts++;
                    }
                    if(constituents[ipart].user_index() == 1){
                        jet_pt_pythia+= constituents[ipart].pt();
                        jet_nparts_pythia++;
                    }
                }

                jet_angularity = temp_angularity;
                jet_track_pt_0 = jet_track_pt[0];
                jet_track_pt_1 = jet_track_pt[1];
                jet_track_pt_2 = jet_track_pt[2];
                jet_track_pt_3 = jet_track_pt[3];
                jet_track_pt_4 = jet_track_pt[4];
                jet_track_pt_5 = jet_track_pt[5];
                jet_track_pt_6 = jet_track_pt[6];
                jet_track_pt_7 = jet_track_pt[7];

                if(jet_nparts == 0) continue;
                // if(jet_nparts_pythia ==0 ) continue;
                // if(jet_pt_pythia < jet) continue;
                
                outTree->Fill();
                nJets_total++;
                constituents.clear();
            }
            

            nEvents_processed++;
            jets.clear();
            particles.clear();
             
    }//end of event/cent loop
    fout->Write();
    fout->Close();
    tenngen_file.Close();
    pythia_file.Close();

    cout << "nEvents_processed: " << nEvents_processed << endl;
    cout << "nJets_total: " << nJets_total << endl;
}

void JetMatch(const TString input_file_pythia, const TString input_file_tenngen, const TString output_file, const Int_t centbin, const struct JetFindingParameters jetparams){
    
    cout << "Starting JetMatch()" << endl;
    cout << "Input pythia file: " << input_file_pythia << endl;
    cout << "Input TennGen file: " << input_file_tenngen << endl;
    cout << "Output file: " << output_file << endl;
    cout << "Jet param: " << jetparams.jetparam << endl;
    cout << "Starting event: " << jetparams.start_mixed_events.at(centbin) << endl;
    cout << "Total events: " << jetparams.n_mixed_events.at(centbin) << endl;
    // open input file
    TFile pythia_file(input_file_pythia.Data());
        TFile tenngen_file(input_file_tenngen.Data());
    if(!pythia_file.IsOpen()){  cout << "PYTHIA file not open" << endl; exit(1); }
    if(!tenngen_file.IsOpen()){  cout << "TennGen file not open" << endl; exit(1); }
    
    TString partTree = "tree";
    TString eventTree = "eventInfo";
    TTreeReader pythia_parts(partTree.Data(), &pythia_file);
    TTreeReader pythia_event(eventTree.Data(), &pythia_file);
    TTreeReaderValue<Int_t> pythia_nevent(pythia_event, "nevent");
    TTreeReaderValue<Float_t> pythia_ptmax(pythia_event, "ptmin");
    TTreeReaderValue<Float_t> pythia_ptmin(pythia_event, "ptmax");
    TTreeReaderValue<Int_t> pythia_ptbin(pythia_event, "ptbinID");
    TTreeReaderValue<Float_t> pythia_xsec(pythia_event, "xsec_over_eventweight");
    TTreeReaderValue<Float_t> pythia_weight(pythia_parts, "weight");
    TTreeReaderValue<Int_t> pythia_numparts(pythia_parts, "nparts");
    TTreeReaderArray<Float_t> pythia_px(pythia_parts, "particle_px");
    TTreeReaderArray<Float_t> pythia_py(pythia_parts, "particle_py");
    TTreeReaderArray<Float_t> pythia_pz(pythia_parts, "particle_pz");
    TTreeReaderArray<Float_t> pythia_E(pythia_parts, "particle_e");
    //////////////////////////////////////////////////////////////
    TTreeReader tenngen_parts(partTree.Data(), &tenngen_file);
    TTreeReaderValue<Int_t> tenngen_numparts(tenngen_parts, "nparts");
    TTreeReaderArray<Float_t> tenngen_px(tenngen_parts, "particle_px");
    TTreeReaderArray<Float_t> tenngen_py(tenngen_parts, "particle_py");
    TTreeReaderArray<Float_t> tenngen_pz(tenngen_parts, "particle_pz");
    TTreeReaderArray<Float_t> tenngen_E(tenngen_parts, "particle_e");
    // check if the trees are empty
    if(!pythia_event.Next()){  cout << "PYTHIA tree is empty" << endl; exit(1); }
    if(!tenngen_parts.Next()){ cout << "TennGen tree is empty" << endl; exit(1); }
    //////////////////////// Create output file //////////////////////////
    const Int_t MAX_TRACKS= 250;
    Int_t event_id_number;
    Int_t event_cent_bin;
    Int_t event_pt_hard_bin;
    Float_t median_pt_over_npart, median_pt_over_area;
    Float_t event_weight;
    Float_t jet_pt_raw;
    Float_t jet_pt_truth;
    Float_t jet_eta, jet_phi, jet_area;
    Float_t jet_pt_pythia;
    Int_t jet_nparts, jet_nparts_pythia;
    // Float_t jet_track_pt[MAX_TRACKS], jet_track_eta[MAX_TRACKS], jet_track_phi[MAX_TRACKS];
    Float_t jet_track_pt_0, jet_track_pt_1, jet_track_pt_2, jet_track_pt_3, jet_track_pt_4, jet_track_pt_5, jet_track_pt_6, jet_track_pt_7;
    Float_t jet_angularity;
    ////////////////////////////////////////////////////////////////////////////////
    // Output file
    TFile *fout = new TFile(output_file.Data(), "RECREATE");
    TTree *outTree = new TTree("tree", "tree");
    outTree->Branch("event_cent_bin", &event_cent_bin, "event_cent_bin/I");
    outTree->Branch("event_pt_hard_bin", &event_pt_hard_bin, "event_pt_hard_bin/I");
    outTree->Branch("median_pt_over_npart", &median_pt_over_npart, "median_pt_over_npart/F");
    outTree->Branch("median_pt_over_area", &median_pt_over_area, "median_pt_over_area/F");
    outTree->Branch("event_weight", &event_weight, "event_weight/F");
    outTree->Branch("jet_pt_raw", &jet_pt_raw, "jet_pt_raw/F");
    outTree->Branch("jet_pt_truth", &jet_pt_truth, "jet_pt_truth/F");
    outTree->Branch("jet_pt_pythia", &jet_pt_pythia, "jet_pt_pythia/F");
    // outTree->Branch("jet_eta", &jet_eta, "jet_eta/F");
    // outTree->Branch("jet_phi", &jet_phi, "jet_phi/F");
    outTree->Branch("jet_area", &jet_area, "jet_area/F");
    outTree->Branch("jet_nparts", &jet_nparts, "jet_nparts/I");
    outTree->Branch("jet_nparts_pythia", &jet_nparts_pythia, "jet_nparts_pythia/I");
    outTree->Branch("jet_angularity", &jet_angularity, "jet_angularity/F");
    outTree->Branch("jet_track_pt_0", &jet_track_pt_0, "jet_track_pt_0/F");
    outTree->Branch("jet_track_pt_1", &jet_track_pt_1, "jet_track_pt_1/F");
    outTree->Branch("jet_track_pt_2", &jet_track_pt_2, "jet_track_pt_2/F");
    outTree->Branch("jet_track_pt_3", &jet_track_pt_3, "jet_track_pt_3/F");
    outTree->Branch("jet_track_pt_4", &jet_track_pt_4, "jet_track_pt_4/F");
    outTree->Branch("jet_track_pt_5", &jet_track_pt_5, "jet_track_pt_5/F");
    outTree->Branch("jet_track_pt_6", &jet_track_pt_6, "jet_track_pt_6/F");
    outTree->Branch("jet_track_pt_7", &jet_track_pt_7, "jet_track_pt_7/F");

    // outTree->Branch("jet_track_pt", jet_track_pt, "jet_track_pt[jet_nparts]/F");
    // outTree->Branch("jet_track_eta", jet_track_eta, "jet_track_eta[jet_nparts]/F");
    // outTree->Branch("jet_track_phi", jet_track_phi, "jet_track_phi[jet_nparts]/F");
    ////////////////////////////////////////////////////////////////////////////////
    // loop over events
    // Get pythia event info

    event_cent_bin = centbin;
    event_pt_hard_bin = *pythia_ptbin;
    Float_t pythia_xsec_over_eventweight = *pythia_xsec;
    Int_t pythia_nevents = *pythia_nevent;
    // Initialize event counters

    Float_t antikt_jet_abs_eta_max = jetparams.particle_max_eta - jetparams.jetparam;
    std::vector<fastjet::PseudoJet> particles;
    std::vector<fastjet::PseudoJet> particles_pythia;
    std::vector<fastjet::PseudoJet> jets;
    std::vector<fastjet::PseudoJet> jets_truth;
    std::vector<fastjet::PseudoJet> constituents;
    std::vector<Float_t> median_pt_over_npart_calc_vec;
    std::vector<Float_t> median_pt_over_area_calc_vec;

    fastjet::GhostedAreaSpec ghost_area_spec(jetparams.particle_max_eta);
    fastjet::JetDefinition kt_jet_def(fastjet::kt_algorithm, 0.4);

    ////////////////////////////////////////////////////////////////////////////////
    // skip pythia events before start_event
    for (Int_t ievent = 0; ievent < jetparams.start_mixed_events.at(centbin); ievent++)
    { 
        if(!pythia_parts.Next()){ cout << "Reached end of pythia file at event " << ievent << endl; exit(1);}
    }
    Int_t nEvents_processed = 0;
    Int_t nJets_total = 0;

    for(Int_t ievent = 0; ievent < jetparams.n_mixed_events.at(centbin); ievent++){
            if(!tenngen_parts.Next()){ cout << "Reached end of tenngen file" << endl; break;}
            if(!pythia_parts.Next()){ cout << "Reached end of pythia file" << endl; break;}
          
            //clear vectors
            particles.clear();
            jets.clear();
            constituents.clear();
            jets_truth.clear();
            particles_pythia.clear();
            median_pt_over_npart_calc_vec.clear();
            median_pt_over_area_calc_vec.clear();
            ////////////////////////////////////////////////////////////////////////////////
            // Fill some tree variables
            event_id_number = ievent;
            // Get event weight
            Float_t temp_event_weight = *pythia_weight;
            event_weight = temp_event_weight*pythia_xsec_over_eventweight;
            ////////////////////////////////////////////////////////////////////////////////
            // Get particle info
            Int_t nParts_pythia = *pythia_numparts;
            Int_t nParts_tenngen = *tenngen_numparts;
            Int_t nParts_total=0;
            for(Int_t ipart = 0; ipart < nParts_pythia; ipart++){   

                PseudoJet particle_temp(pythia_px[ipart], pythia_py[ipart], pythia_pz[ipart], pythia_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles_pythia.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(1);
                particles_pythia.at(nParts_total).set_user_index(1);
                nParts_total++;

            }
            for(Int_t ipart = 0; ipart < nParts_tenngen; ipart++){
                
                PseudoJet particle_temp(tenngen_px[ipart], tenngen_py[ipart], tenngen_pz[ipart], tenngen_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(0);
                nParts_total++;

            }

            ////////////////////////////////////////////////////////////////////////////////
            // find kt jets
            fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);  
            fastjet::ClusterSequenceArea kt_cs(particles, kt_jet_def, area_def);
            fastjet::Selector jet_selector_kt = (!SelectorNHardest(2));
            // clear vectors
            jets.clear();
            median_pt_over_area_calc_vec.clear();
            median_pt_over_npart_calc_vec.clear();
            // get kt jets
            jets = jet_selector_kt(kt_cs.inclusive_jets());
            // find median kt jet pt
            Int_t nktjets = jets.size();
            for(Int_t ijet = 0; ijet < nktjets; ijet++){

                if(!jets[ijet].has_area()) continue;
                Float_t kt_jet_area = jets[ijet].area();

                if(jets[ijet].constituents().size() == 0) continue;
                Int_t kt_jet_nparts = jets[ijet].constituents().size();

                if(jets[ijet].pt() == 0) continue;
                Float_t kt_jet_pt = jets[ijet].pt();
                
                Float_t pt_over_area = 1.0*(kt_jet_pt/kt_jet_area);
                Float_t pt_over_npart = 1.0*(kt_jet_pt/kt_jet_nparts);
                median_pt_over_area_calc_vec.push_back(pt_over_area);
                median_pt_over_npart_calc_vec.push_back(pt_over_npart);

            }
            // sort vectors
            std::sort(median_pt_over_npart_calc_vec.begin(), median_pt_over_npart_calc_vec.end());
            std::sort(median_pt_over_area_calc_vec.begin(), median_pt_over_area_calc_vec.end());
            median_pt_over_npart = 0;
            median_pt_over_area = 0;
            // find median
            if(median_pt_over_npart_calc_vec.size() > 0) median_pt_over_npart = median_pt_over_npart_calc_vec[median_pt_over_npart_calc_vec.size()/2];
            if(median_pt_over_area_calc_vec.size() > 0) median_pt_over_area = median_pt_over_area_calc_vec[median_pt_over_area_calc_vec.size()/2];
            //clear vectors
            jets.clear();
            median_pt_over_npart_calc_vec.clear();
            median_pt_over_area_calc_vec.clear();

            ////////////////////////////////////////////////////////////////////////////////
            // find antikt jets
            fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparams.jetparam);
            fastjet::ClusterSequenceArea antikt_cs(particles, antikt_jet_def, area_def);
            fastjet::ClusterSequenceArea antikt_cs_pythia(particles_pythia, antikt_jet_def, area_def);
            jets = antikt_cs.inclusive_jets(jetparams.jet_pt_min);
            jets_truth = antikt_cs_pythia.inclusive_jets(jetparams.jet_pt_min);
            // loop over truth antikt jets
            Int_t n_antikt_jets = jets.size();
            Int_t n_antikt_jets_truth = jets_truth.size();
            Int_t matchedJets = 0;
            Int_t n_matched_jets_no_bijections = 0;
            for(Int_t itruth =0; itruth < n_antikt_jets_truth; itruth++){

                // check if truth jet is in acceptance
                if(TMath::Abs(jets_truth[itruth].eta()) > antikt_jet_abs_eta_max) continue;
                if(jets_truth[itruth].pt() < jetparams.jet_pt_min) continue;
                if(jets_truth[itruth].constituents().size() == 0) continue;
                //if(!jets_truth[itruth].has_area()) continue;

                Float_t dR_minimum = 0.101;
                Int_t matched_reco_index = -1;
                Int_t matched_pythia_index = -1;
                Float_t dR_minimum_pythia = 0.101;
                
                // loop over reco antikt jets
                for(Int_t ireco =0; ireco<n_antikt_jets; ireco++){
                    // check if reco jet is in acceptance
                    if(TMath::Abs(jets[ireco].eta()) > antikt_jet_abs_eta_max) continue;
                    if(jets[ireco].pt() < jetparams.jet_pt_min) continue;
                    if(jets[ireco].constituents().size() == 0) continue;
                    //if(!jets[ireco].has_area()) continue;

                    Float_t deta = jets_truth[itruth].eta() - jets[ireco].eta();
                    Float_t dphi = jets_truth[itruth].phi() - jets[ireco].phi();
                    if(dphi < 0) dphi = -1.0*dphi;
                    if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                    Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
                    if(dR<dR_minimum){
                        dR_minimum = dR;
                        matched_reco_index = ireco;
                    }
                }

                if(matched_reco_index == -1) continue;
                n_matched_jets_no_bijections++;
                // cross check
                for(Int_t ipythia_check =0; ipythia_check < n_antikt_jets_truth; ipythia_check++){
                    //if(ipythia_check == itruth) continue;
                    // check if truth jet is in acceptance
                    if(TMath::Abs(jets_truth[ipythia_check].eta()) > antikt_jet_abs_eta_max) continue;
                    if(jets_truth[ipythia_check].pt() < jetparams.jet_pt_min) continue;
                    if(jets_truth[ipythia_check].constituents().size() == 0) continue;
                    //if(!jets_truth[ipythia_check].has_area()) continue;

                    Float_t deta_check = jets[matched_reco_index].eta() - jets_truth[ipythia_check].eta();
                    Float_t dphi_check = jets[matched_reco_index].phi() - jets_truth[ipythia_check].phi();
                    if(dphi_check < 0) dphi_check = -1.0*dphi_check;
                    if(dphi_check > TMath::Pi()) dphi_check = 2*TMath::Pi() - dphi_check;
                    Float_t dR_check = TMath::Sqrt(deta_check*deta_check + dphi_check*dphi_check);
                    if(dR_check<dR_minimum_pythia){
                        dR_minimum_pythia = dR_check;
                        matched_pythia_index = ipythia_check;
                    }
                }

                if(matched_pythia_index != itruth) continue;
                matchedJets++;

                constituents.clear();
                jet_pt_pythia = 0;
                jet_pt_raw = jets[matched_reco_index].pt();
                jet_eta = jets[matched_reco_index].eta();
                jet_phi = jets[matched_reco_index].phi();
                jet_area = jets[matched_reco_index].area();
                jet_pt_truth = jets_truth[matched_pythia_index].pt();
                jet_nparts = 0;
                jet_nparts_pythia = 0;
                Float_t jet_track_pt[8] = {0,0,0,0,0,0,0,0};
                Float_t temp_angularity = 0;
                constituents = sorted_by_pt(jets[matched_reco_index].constituents());
                Int_t n_constituents = constituents.size();
                Int_t active_tracks = 0;
                for(Int_t ipart =0; ipart< n_constituents; ipart++){
                    if(constituents[ipart].user_index() == 1 || constituents[ipart].user_index() == 0){
                        Float_t deta = constituents[ipart].eta() - jet_eta;
                        Float_t dphi = constituents[ipart].phi() - jet_phi;
                        if(dphi < 0) dphi = -1.0*dphi;
                        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                        Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
                        temp_angularity += (dR/jetparams.jetparam)*(constituents[ipart].pt()/jet_pt_raw);
                        if(jet_nparts < 8){
                            jet_track_pt[jet_nparts] = constituents[ipart].pt();
                        }
                        jet_nparts++;
                    }
                    if(constituents[ipart].user_index() == 1){
                        jet_pt_pythia+= constituents[ipart].pt();
                        jet_nparts_pythia++;
                    }
                }

                jet_angularity = temp_angularity;
                jet_track_pt_0 = jet_track_pt[0];
                jet_track_pt_1 = jet_track_pt[1];
                jet_track_pt_2 = jet_track_pt[2];
                jet_track_pt_3 = jet_track_pt[3];
                jet_track_pt_4 = jet_track_pt[4];
                jet_track_pt_5 = jet_track_pt[5];
                jet_track_pt_6 = jet_track_pt[6];
                jet_track_pt_7 = jet_track_pt[7];

                // fill tree
                outTree->Fill();

                // clear vectors
                constituents.clear();
                // reset variable
                nJets_total++;

                
            }
            
            // clear vectors
            jets.clear();
            jets_truth.clear();
            particles.clear();
            particles_pythia.clear();
            nEvents_processed++;
            if(ievent%10000==0) cout<<"Event "<<ievent <<" contained " << matchedJets << " matched jets" << endl;
             
    }//end of event/cent loop


    fout->Write();
    fout->Close();
    tenngen_file.Close();
    pythia_file.Close();
    cout << "Total number of events processed: " << nEvents_processed << endl;
    cout << "Total number of matched jets: " << nJets_total << endl;

}

void MissedPythiaJets(const TString input_file_pythia, const TString input_file_tenngen, const TString output_file, const Int_t centbin, const struct JetFindingParameters jetparams){
    
    cout << "Starting MissedPythiaJets()" << endl;
    cout << "Input pythia file: " << input_file_pythia << endl;
    cout << "Input TennGen file: " << input_file_tenngen << endl;
    cout << "Output file: " << output_file << endl;
    cout << "Jet param: " << jetparams.jetparam << endl;
    cout << "Starting event: " << jetparams.start_mixed_events.at(centbin) << endl;
    cout << "Total events: " << jetparams.n_mixed_events.at(centbin) << endl;
    // open input file
    TFile pythia_file(input_file_pythia.Data());
        TFile tenngen_file(input_file_tenngen.Data());
    if(!pythia_file.IsOpen()){  cout << "PYTHIA file not open" << endl; exit(1); }
    if(!tenngen_file.IsOpen()){  cout << "TennGen file not open" << endl; exit(1); }
    
    TString partTree = "tree";
    TString eventTree = "eventInfo";
    TTreeReader pythia_parts(partTree.Data(), &pythia_file);
    TTreeReader pythia_event(eventTree.Data(), &pythia_file);

    TTreeReaderValue<Int_t> pythia_ptbin(pythia_event, "ptbinID");
    TTreeReaderValue<Int_t> pythia_numparts(pythia_parts, "nparts");
    TTreeReaderArray<Float_t> pythia_px(pythia_parts, "particle_px");
    TTreeReaderArray<Float_t> pythia_py(pythia_parts, "particle_py");
    TTreeReaderArray<Float_t> pythia_pz(pythia_parts, "particle_pz");
    TTreeReaderArray<Float_t> pythia_E(pythia_parts, "particle_e");
    //////////////////////////////////////////////////////////////
    TTreeReader tenngen_parts(partTree.Data(), &tenngen_file);
    TTreeReaderValue<Int_t> tenngen_numparts(tenngen_parts, "nparts");
    TTreeReaderArray<Float_t> tenngen_px(tenngen_parts, "particle_px");
    TTreeReaderArray<Float_t> tenngen_py(tenngen_parts, "particle_py");
    TTreeReaderArray<Float_t> tenngen_pz(tenngen_parts, "particle_pz");
    TTreeReaderArray<Float_t> tenngen_E(tenngen_parts, "particle_e");
    // check if the trees are empty
    if(!pythia_event.Next()){  cout << "PYTHIA tree is empty" << endl; exit(1); }
    if(!tenngen_parts.Next()){ cout << "TennGen tree is empty" << endl; exit(1); }
    //////////////////////// Create output file //////////////////////////
    const Int_t MAX_TRACKS= 250;
    Int_t event_cent_bin;
    Int_t event_pt_hard_bin;
    Float_t jet_pt_truth;

    ////////////////////////////////////////////////////////////////////////////////
    // Output file
    TFile *fout = new TFile(output_file.Data(), "RECREATE");
    TTree *outTree = new TTree("tree", "tree");
    outTree->Branch("event_cent_bin", &event_cent_bin, "event_cent_bin/I");
    outTree->Branch("event_pt_hard_bin", &event_pt_hard_bin, "event_pt_hard_bin/I");
    outTree->Branch("jet_pt_truth", &jet_pt_truth, "jet_pt_truth/F");
    ////////////////////////////////////////////////////////////////////////////////
    // loop over events
    // Get pythia event info

    event_cent_bin = centbin;
    event_pt_hard_bin = *pythia_ptbin;
    // Initialize event counters

    Float_t antikt_jet_abs_eta_max = jetparams.particle_max_eta - jetparams.jetparam;
    std::vector<fastjet::PseudoJet> particles;
    std::vector<fastjet::PseudoJet> particles_pythia;
    std::vector<fastjet::PseudoJet> jets;
    std::vector<fastjet::PseudoJet> jets_truth;

    fastjet::GhostedAreaSpec ghost_area_spec(jetparams.particle_max_eta);
    fastjet::JetDefinition kt_jet_def(fastjet::kt_algorithm, 0.4);

    ////////////////////////////////////////////////////////////////////////////////
    // skip pythia events before start_event
    for (Int_t ievent = 0; ievent < jetparams.start_mixed_events.at(centbin); ievent++)
    { 
        if(!pythia_parts.Next()){ cout << "Reached end of pythia file at event " << ievent << endl; exit(1);}
    }
    Int_t nEvents_processed = 0;
    Int_t nJets_total = 0;

    for(Int_t ievent = 0; ievent < jetparams.n_mixed_events.at(centbin); ievent++){
            if(!tenngen_parts.Next()){ cout << "Reached end of tenngen file" << endl; break;}
            if(!pythia_parts.Next()){ cout << "Reached end of pythia file" << endl; break;}
          
            //clear vectors
            particles.clear();
            jets.clear();
            jets_truth.clear();
            particles_pythia.clear();
            ////////////////////////////////////////////////////////////////////////////////
            // Get particle info
            Int_t nParts_pythia = *pythia_numparts;
            Int_t nParts_tenngen = *tenngen_numparts;
            Int_t nParts_total=0;
            for(Int_t ipart = 0; ipart < nParts_pythia; ipart++){   

                PseudoJet particle_temp(pythia_px[ipart], pythia_py[ipart], pythia_pz[ipart], pythia_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles_pythia.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(1);
                particles_pythia.at(nParts_total).set_user_index(1);
                nParts_total++;

            }
            for(Int_t ipart = 0; ipart < nParts_tenngen; ipart++){
                
                PseudoJet particle_temp(tenngen_px[ipart], tenngen_py[ipart], tenngen_pz[ipart], tenngen_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(0);
                nParts_total++;

            }

            ////////////////////////////////////////////////////////////////////////////////
            // find antikt jets
            fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);  
            fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparams.jetparam);
            fastjet::ClusterSequenceArea antikt_cs(particles, antikt_jet_def, area_def);
            fastjet::ClusterSequenceArea antikt_cs_pythia(particles_pythia, antikt_jet_def, area_def);
            jets = antikt_cs.inclusive_jets(jetparams.jet_pt_min);
            jets_truth = antikt_cs_pythia.inclusive_jets(jetparams.jet_pt_min);
            // loop over truth antikt jets
            Int_t n_antikt_jets = jets.size();
            Int_t n_antikt_jets_truth = jets_truth.size();
            Int_t missedJets = 0;
            for(Int_t itruth =0; itruth < n_antikt_jets_truth; itruth++){

                Int_t jet_is_matched = 1;
                // check if truth jet is in acceptance
                if(TMath::Abs(jets_truth[itruth].eta()) > antikt_jet_abs_eta_max) continue;
                if(jets_truth[itruth].pt() < jetparams.jet_pt_min) continue;
                if(jets_truth[itruth].constituents().size() == 0) continue;
                //if(!jets_truth[itruth].has_area()) continue;

                Float_t dR_minimum = 0.101;
                Int_t matched_reco_index = -1;
                Int_t matched_pythia_index = -1;
                Float_t dR_minimum_pythia = 0.101;
                
                // loop over reco antikt jets
                for(Int_t ireco =0; ireco<n_antikt_jets; ireco++){
                    // check if reco jet is in acceptance
                    if(TMath::Abs(jets[ireco].eta()) > antikt_jet_abs_eta_max) continue;
                    if(jets[ireco].pt() < jetparams.jet_pt_min) continue;
                    if(jets[ireco].constituents().size() == 0) continue;
                    //if(!jets[ireco].has_area()) continue;

                    Float_t deta = jets_truth[itruth].eta() - jets[ireco].eta();
                    Float_t dphi = jets_truth[itruth].phi() - jets[ireco].phi();
                    if(dphi < 0) dphi = -1.0*dphi;
                    if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                    Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
                    if(dR<dR_minimum){
                        dR_minimum = dR;
                        matched_reco_index = ireco;
                    }
                }

                if (matched_reco_index == -1){
                    jet_is_matched = 0;
                }
                else if (matched_reco_index != -1){
                    // cross check
                    for(Int_t ipythia_check =0; ipythia_check < n_antikt_jets_truth; ipythia_check++){
                        // check if truth jet is in acceptance
                        if(TMath::Abs(jets_truth[ipythia_check].eta()) > antikt_jet_abs_eta_max) continue;
                        if(jets_truth[ipythia_check].pt() < jetparams.jet_pt_min) continue;
                        if(jets_truth[ipythia_check].constituents().size() == 0) continue;
                        //if(!jets_truth[ipythia_check].has_area()) continue;

                        Float_t deta_check = jets[matched_reco_index].eta() - jets_truth[ipythia_check].eta();
                        Float_t dphi_check = jets[matched_reco_index].phi() - jets_truth[ipythia_check].phi();
                        if(dphi_check < 0) dphi_check = -1.0*dphi_check;
                        if(dphi_check > TMath::Pi()) dphi_check = 2*TMath::Pi() - dphi_check;
                        Float_t dR_check = TMath::Sqrt(deta_check*deta_check + dphi_check*dphi_check);
                        if(dR_check<dR_minimum_pythia){
                            dR_minimum_pythia = dR_check;
                            matched_pythia_index = ipythia_check;
                        }
                    }

                    if (matched_pythia_index != itruth) jet_is_matched = 0;
                }
                // only fill tree if jet is not matched to another jet
                if(jet_is_matched) continue;

                missedJets++;
                nJets_total++;

                jet_pt_truth = jets_truth[itruth].pt();
                // fill tree
                outTree->Fill();                
            }
            
            // clear vectors
            jets.clear();
            jets_truth.clear();
            particles.clear();
            particles_pythia.clear();
            nEvents_processed++;
            if(ievent%10000==0) cout<<"Event "<<ievent <<" contained " << missedJets << " missed jets" << endl;
             
    }//end of event/cent loop


    fout->Write();
    fout->Close();
    tenngen_file.Close();
    pythia_file.Close();
    cout << "Total number of events processed: " << nEvents_processed << endl;
    cout << "Total number of missed jets: " << nJets_total << endl;

}

void UnmatchedJets(const TString input_file_pythia, const TString input_file_tenngen, const TString output_file, const Int_t centbin, const struct JetFindingParameters jetparams){
    
    cout << "Starting UnmatchedJets()" << endl;
    cout << "Input pythia file: " << input_file_pythia << endl;
    cout << "Input TennGen file: " << input_file_tenngen << endl;
    cout << "Output file: " << output_file << endl;
    cout << "Jet param: " << jetparams.jetparam << endl;
    cout << "Starting event: " << jetparams.start_mixed_events.at(centbin) << endl;
    cout << "Total events: " << jetparams.n_mixed_events.at(centbin) << endl;
    // open input file
    TFile pythia_file(input_file_pythia.Data());
        TFile tenngen_file(input_file_tenngen.Data());
    if(!pythia_file.IsOpen()){  cout << "PYTHIA file not open" << endl; exit(1); }
    if(!tenngen_file.IsOpen()){  cout << "TennGen file not open" << endl; exit(1); }
    
    TString partTree = "tree";
    TString eventTree = "eventInfo";
    TTreeReader pythia_parts(partTree.Data(), &pythia_file);
    TTreeReader pythia_event(eventTree.Data(), &pythia_file);
    TTreeReaderValue<Int_t> pythia_nevent(pythia_event, "nevent");
    TTreeReaderValue<Float_t> pythia_ptmax(pythia_event, "ptmin");
    TTreeReaderValue<Float_t> pythia_ptmin(pythia_event, "ptmax");
    TTreeReaderValue<Int_t> pythia_ptbin(pythia_event, "ptbinID");
    TTreeReaderValue<Float_t> pythia_xsec(pythia_event, "xsec_over_eventweight");
    TTreeReaderValue<Float_t> pythia_weight(pythia_parts, "weight");
    TTreeReaderValue<Int_t> pythia_numparts(pythia_parts, "nparts");
    TTreeReaderArray<Float_t> pythia_px(pythia_parts, "particle_px");
    TTreeReaderArray<Float_t> pythia_py(pythia_parts, "particle_py");
    TTreeReaderArray<Float_t> pythia_pz(pythia_parts, "particle_pz");
    TTreeReaderArray<Float_t> pythia_E(pythia_parts, "particle_e");
    //////////////////////////////////////////////////////////////
    TTreeReader tenngen_parts(partTree.Data(), &tenngen_file);
    TTreeReaderValue<Int_t> tenngen_numparts(tenngen_parts, "nparts");
    TTreeReaderArray<Float_t> tenngen_px(tenngen_parts, "particle_px");
    TTreeReaderArray<Float_t> tenngen_py(tenngen_parts, "particle_py");
    TTreeReaderArray<Float_t> tenngen_pz(tenngen_parts, "particle_pz");
    TTreeReaderArray<Float_t> tenngen_E(tenngen_parts, "particle_e");
    // check if the trees are empty
    if(!pythia_event.Next()){  cout << "PYTHIA tree is empty" << endl; exit(1); }
    if(!tenngen_parts.Next()){ cout << "TennGen tree is empty" << endl; exit(1); }
    //////////////////////// Create output file //////////////////////////
    const Int_t MAX_TRACKS= 250;
    Int_t event_id_number;
    Int_t event_cent_bin;
    Int_t event_pt_hard_bin;
    Float_t median_pt_over_npart, median_pt_over_area;
    Float_t event_weight;
    Float_t jet_pt_raw;
    Float_t jet_pt_truth;
    Float_t jet_eta, jet_phi, jet_area;
    Float_t jet_pt_pythia;
    Int_t jet_nparts, jet_nparts_pythia;
    // Float_t jet_track_pt[MAX_TRACKS], jet_track_eta[MAX_TRACKS], jet_track_phi[MAX_TRACKS];
    Float_t jet_track_pt_0, jet_track_pt_1, jet_track_pt_2, jet_track_pt_3, jet_track_pt_4, jet_track_pt_5, jet_track_pt_6, jet_track_pt_7;
    Float_t jet_angularity;
    ////////////////////////////////////////////////////////////////////////////////
    // Output file
    TFile *fout = new TFile(output_file.Data(), "RECREATE");
    TTree *outTree = new TTree("tree", "tree");
    outTree->Branch("event_cent_bin", &event_cent_bin, "event_cent_bin/I");
    outTree->Branch("event_pt_hard_bin", &event_pt_hard_bin, "event_pt_hard_bin/I");
    outTree->Branch("median_pt_over_npart", &median_pt_over_npart, "median_pt_over_npart/F");
    outTree->Branch("median_pt_over_area", &median_pt_over_area, "median_pt_over_area/F");
    outTree->Branch("event_weight", &event_weight, "event_weight/F");
    outTree->Branch("jet_pt_raw", &jet_pt_raw, "jet_pt_raw/F");
    outTree->Branch("jet_pt_truth", &jet_pt_truth, "jet_pt_truth/F");
    outTree->Branch("jet_pt_pythia", &jet_pt_pythia, "jet_pt_pythia/F");
    outTree->Branch("jet_area", &jet_area, "jet_area/F");
    outTree->Branch("jet_nparts", &jet_nparts, "jet_nparts/I");
    outTree->Branch("jet_nparts_pythia", &jet_nparts_pythia, "jet_nparts_pythia/I");
    outTree->Branch("jet_angularity", &jet_angularity, "jet_angularity/F");
    outTree->Branch("jet_track_pt_0", &jet_track_pt_0, "jet_track_pt_0/F");
    outTree->Branch("jet_track_pt_1", &jet_track_pt_1, "jet_track_pt_1/F");
    outTree->Branch("jet_track_pt_2", &jet_track_pt_2, "jet_track_pt_2/F");
    outTree->Branch("jet_track_pt_3", &jet_track_pt_3, "jet_track_pt_3/F");
    outTree->Branch("jet_track_pt_4", &jet_track_pt_4, "jet_track_pt_4/F");
    outTree->Branch("jet_track_pt_5", &jet_track_pt_5, "jet_track_pt_5/F");
    outTree->Branch("jet_track_pt_6", &jet_track_pt_6, "jet_track_pt_6/F");
    outTree->Branch("jet_track_pt_7", &jet_track_pt_7, "jet_track_pt_7/F");
  ////////////////////////////////////////////////////////////////////////////////
    // loop over events
    // Get pythia event info

    event_cent_bin = centbin;
    event_pt_hard_bin = *pythia_ptbin;
    Float_t pythia_xsec_over_eventweight = *pythia_xsec;
    Int_t pythia_nevents = *pythia_nevent;
    // Initialize event counters

    Float_t antikt_jet_abs_eta_max = jetparams.particle_max_eta - jetparams.jetparam;
    std::vector<fastjet::PseudoJet> particles;
    std::vector<fastjet::PseudoJet> particles_pythia;
    std::vector<fastjet::PseudoJet> jets;
    std::vector<fastjet::PseudoJet> jets_truth;
    std::vector<fastjet::PseudoJet> constituents;
    std::vector<Float_t> median_pt_over_npart_calc_vec;
    std::vector<Float_t> median_pt_over_area_calc_vec;

    fastjet::GhostedAreaSpec ghost_area_spec(jetparams.particle_max_eta);
    fastjet::JetDefinition kt_jet_def(fastjet::kt_algorithm, 0.4);

    ////////////////////////////////////////////////////////////////////////////////
    // skip pythia events before start_event
    for (Int_t ievent = 0; ievent < jetparams.start_mixed_events.at(centbin); ievent++)
    { 
        if(!pythia_parts.Next()){ cout << "Reached end of pythia file at event " << ievent << endl; exit(1);}
    }
    Int_t nEvents_processed = 0;
    Int_t nJets_total = 0;

    for(Int_t ievent = 0; ievent < jetparams.n_mixed_events.at(centbin); ievent++){
            if(!tenngen_parts.Next()){ cout << "Reached end of tenngen file" << endl; break;}
            if(!pythia_parts.Next()){ cout << "Reached end of pythia file" << endl; break;}
          
            //clear vectors
            particles.clear();
            jets.clear();
            constituents.clear();
            jets_truth.clear();
            particles_pythia.clear();
            median_pt_over_npart_calc_vec.clear();
            median_pt_over_area_calc_vec.clear();
            ////////////////////////////////////////////////////////////////////////////////
            // Fill some tree variables
            event_id_number = ievent;
            // Get event weight
            Float_t temp_event_weight = *pythia_weight;
            event_weight = temp_event_weight*pythia_xsec_over_eventweight;
            ////////////////////////////////////////////////////////////////////////////////
            // Get particle info
            Int_t nParts_pythia = *pythia_numparts;
            Int_t nParts_tenngen = *tenngen_numparts;
            Int_t nParts_total=0;
            for(Int_t ipart = 0; ipart < nParts_pythia; ipart++){   

                PseudoJet particle_temp(pythia_px[ipart], pythia_py[ipart], pythia_pz[ipart], pythia_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles_pythia.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(1);
                particles_pythia.at(nParts_total).set_user_index(1);
                nParts_total++;

            }
            for(Int_t ipart = 0; ipart < nParts_tenngen; ipart++){
                
                PseudoJet particle_temp(tenngen_px[ipart], tenngen_py[ipart], tenngen_pz[ipart], tenngen_E[ipart]);
                if(particle_temp.pt() < 0.15) continue;
                if(TMath::Abs(particle_temp.eta()) > jetparams.particle_max_eta) continue;
                if(particle_temp.px() == 0 && particle_temp.py() == 0 && particle_temp.pz() == 0 && particle_temp.E() == 0) continue;
                particles.push_back(particle_temp);
                particles.at(nParts_total).set_user_index(0);
                nParts_total++;

            }

            ////////////////////////////////////////////////////////////////////////////////
            // find kt jets
            fastjet::AreaDefinition area_def(fastjet::active_area, ghost_area_spec);  
            fastjet::ClusterSequenceArea kt_cs(particles, kt_jet_def, area_def);
            fastjet::Selector jet_selector_kt = (!SelectorNHardest(2));
            // clear vectors
            jets.clear();
            median_pt_over_area_calc_vec.clear();
            median_pt_over_npart_calc_vec.clear();
            // get kt jets
            jets = jet_selector_kt(kt_cs.inclusive_jets());
            // find median kt jet pt
            Int_t nktjets = jets.size();
            for(Int_t ijet = 0; ijet < nktjets; ijet++){

                if(!jets[ijet].has_area()) continue;
                Float_t kt_jet_area = jets[ijet].area();

                if(jets[ijet].constituents().size() == 0) continue;
                Int_t kt_jet_nparts = jets[ijet].constituents().size();

                if(jets[ijet].pt() == 0) continue;
                Float_t kt_jet_pt = jets[ijet].pt();
                
                Float_t pt_over_area = 1.0*(kt_jet_pt/kt_jet_area);
                Float_t pt_over_npart = 1.0*(kt_jet_pt/kt_jet_nparts);
                median_pt_over_area_calc_vec.push_back(pt_over_area);
                median_pt_over_npart_calc_vec.push_back(pt_over_npart);

            }
            // sort vectors
            std::sort(median_pt_over_npart_calc_vec.begin(), median_pt_over_npart_calc_vec.end());
            std::sort(median_pt_over_area_calc_vec.begin(), median_pt_over_area_calc_vec.end());
            median_pt_over_npart = 0;
            median_pt_over_area = 0;
            // find median
            if(median_pt_over_npart_calc_vec.size() > 0) median_pt_over_npart = median_pt_over_npart_calc_vec[median_pt_over_npart_calc_vec.size()/2];
            if(median_pt_over_area_calc_vec.size() > 0) median_pt_over_area = median_pt_over_area_calc_vec[median_pt_over_area_calc_vec.size()/2];
            //clear vectors
            jets.clear();
            median_pt_over_npart_calc_vec.clear();
            median_pt_over_area_calc_vec.clear();

            ////////////////////////////////////////////////////////////////////////////////
            // find antikt jets
            fastjet::JetDefinition antikt_jet_def(fastjet::antikt_algorithm, jetparams.jetparam);
            fastjet::ClusterSequenceArea antikt_cs(particles, antikt_jet_def, area_def);
            fastjet::ClusterSequenceArea antikt_cs_pythia(particles_pythia, antikt_jet_def, area_def);
            jets = antikt_cs.inclusive_jets(jetparams.jet_pt_min);
            jets_truth = antikt_cs_pythia.inclusive_jets(jetparams.jet_pt_min);
            // loop over truth antikt jets
            Int_t n_antikt_jets = jets.size();
            Int_t n_antikt_jets_truth = jets_truth.size();
            Int_t unmatchedJets = 0;
            for(Int_t ireco =0; ireco < n_antikt_jets; ireco++){

                Int_t reco_jet_is_matched = 1;
                // check if reco jet is in acceptance
                if(TMath::Abs(jets[ireco].eta()) > antikt_jet_abs_eta_max) continue;
                if(jets[ireco].pt() < jetparams.jet_pt_min) continue;
                if(jets[ireco].constituents().size() == 0) continue;

                Float_t dR_minimum = 0.101;
                Int_t matched_reco_index = -1;
                Int_t matched_pythia_index = -1;
                Float_t dR_minimum_pythia = 0.101;
                
                // loop over pythia truth jets
                for(Int_t itruth =0; itruth<n_antikt_jets_truth; itruth++){
                    // check if reco jet is in acceptance
                    if(TMath::Abs(jets_truth[itruth].eta()) > antikt_jet_abs_eta_max) continue;
                    if(jets_truth[itruth].pt() < jetparams.jet_pt_min) continue;
                    if(jets_truth[itruth].constituents().size() == 0) continue;

                    Float_t deta = jets[ireco].eta() - jets_truth[itruth].eta();
                    Float_t dphi = jets[ireco].phi() - jets_truth[itruth].phi();
                    if(dphi < 0) dphi = -1.0*dphi;
                    if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                    Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
                    if(dR<dR_minimum){
                        dR_minimum = dR;
                        matched_pythia_index = itruth;
                    }
                }

                if(matched_pythia_index == -1){
                    reco_jet_is_matched = 0;
                }
                else if(matched_pythia_index != -1){ 
                    // cross check
                    for(Int_t ireco_check =0; ireco_check < n_antikt_jets; ireco_check++){
                        // check if truth jet is in acceptance
                        if(TMath::Abs(jets[ireco_check].eta()) > antikt_jet_abs_eta_max) continue;
                        if(jets[ireco_check].pt() < jetparams.jet_pt_min) continue;
                        if(jets[ireco_check].constituents().size() == 0) continue;
                        //if(!jets_truth[ireco_check].has_area()) continue;

                        Float_t deta_check = jets_truth[matched_pythia_index].eta() - jets[ireco_check].eta();
                        Float_t dphi_check = jets_truth[matched_pythia_index].phi() - jets[ireco_check].phi();
                        if(dphi_check < 0) dphi_check = -1.0*dphi_check;
                        if(dphi_check > TMath::Pi()) dphi_check = 2*TMath::Pi() - dphi_check;
                        Float_t dR_check = TMath::Sqrt(deta_check*deta_check + dphi_check*dphi_check);
                        if(dR_check<dR_minimum_pythia){
                            dR_minimum_pythia = dR_check;
                            matched_reco_index = ireco_check;
                        }
                    }

                    if(matched_reco_index != ireco){
                        reco_jet_is_matched = 0;
                    }
                }

                if(reco_jet_is_matched) continue;

                constituents.clear();
                jet_pt_pythia = 0;
                jet_pt_raw = jets[ireco].pt();
                jet_eta = jets[ireco].eta();
                jet_phi = jets[ireco].phi();
                jet_area = jets[ireco].area();
                jet_nparts = 0;
                jet_nparts_pythia = 0;
                Float_t jet_track_pt[8] = {0,0,0,0,0,0,0,0};
                Float_t temp_angularity = 0;
                constituents = sorted_by_pt(jets[ireco].constituents());
                Int_t n_constituents = constituents.size();
                for(Int_t ipart =0; ipart< n_constituents; ipart++){
                    if(constituents[ipart].user_index() == 1 || constituents[ipart].user_index() == 0){
                        Float_t deta = constituents[ipart].eta() - jet_eta;
                        Float_t dphi = constituents[ipart].phi() - jet_phi;
                        if(dphi < 0) dphi = -1.0*dphi;
                        if(dphi > TMath::Pi()) dphi = 2*TMath::Pi() - dphi;
                        Float_t dR = TMath::Sqrt(deta*deta + dphi*dphi);
                        temp_angularity += (dR/jetparams.jetparam)*(constituents[ipart].pt()/jet_pt_raw);
                        if(jet_nparts < 8){
                            jet_track_pt[jet_nparts] = constituents[ipart].pt();
                        }
                        jet_nparts++;
                    }
                    if(constituents[ipart].user_index() == 1){
                        jet_pt_pythia+= constituents[ipart].pt();
                        jet_nparts_pythia++;
                    }
                }

                jet_angularity = temp_angularity;
                jet_track_pt_0 = jet_track_pt[0];
                jet_track_pt_1 = jet_track_pt[1];
                jet_track_pt_2 = jet_track_pt[2];
                jet_track_pt_3 = jet_track_pt[3];
                jet_track_pt_4 = jet_track_pt[4];
                jet_track_pt_5 = jet_track_pt[5];
                jet_track_pt_6 = jet_track_pt[6];
                jet_track_pt_7 = jet_track_pt[7];

                // fill tree
                outTree->Fill();

                // clear vectors
                constituents.clear();
                nJets_total++;      
                unmatchedJets++;          
            }
            
            // clear vectors
            jets.clear();
            jets_truth.clear();
            particles.clear();
            particles_pythia.clear();
            nEvents_processed++;
            if(ievent%10000==0) cout<<"Event "<<ievent <<" contained " << unmatchedJets << " unmatched jets" << endl;
             
    }//end of event/cent loop


    fout->Write();
    fout->Close();
    tenngen_file.Close();
    pythia_file.Close();
    cout << "Total number of events processed: " << nEvents_processed << endl;
    cout << "Total number of unmatched jets: " << nJets_total << endl;

}

TString GetOutFilePath(const TString output_path, const Int_t collen, const Float_t jetparam, const TString mode){

    if(gSystem->AccessPathName(output_path)) gSystem->mkdir(output_path);
    TString output_dir = output_path;
    if(collen == 200) output_dir += "200GeV/";
    else if(collen == 2760) output_dir += "2760GeV/";
    if(gSystem->AccessPathName(output_dir)) gSystem->mkdir(output_dir);
    output_dir += Form("R0%.0f/", jetparam*10);
    if(gSystem->AccessPathName(output_dir)) gSystem->mkdir(output_dir);
    output_dir += "Raw/";
    if(gSystem->AccessPathName(output_dir)) gSystem->mkdir(output_dir);
    output_dir += Form("%s/", mode.Data());
    if(gSystem->AccessPathName(output_dir)) gSystem->mkdir(output_dir);
    return output_dir;

}

void Usage(){
    cout << "Usage: ./JetFinder <collen> <ptbin> <jetparam> <mode> <optional nevent> <verbose>" << endl;
    cout << "collen: 200 or 2760" << endl;
    cout << "ptbin: 0-24" << endl;
    cout << "jetparam: 0.2-0.6" << endl;
    cout << "mode: PP, FullJets, Match" << endl;
    cout << "nevent: -1 for all events" << endl;
    cout << "verbose: 0 for limited output, 1 for debug output" << endl;
    exit(0);
}

void InitJetAnalysis(Int_t argc, Char_t** argv,  struct Args &args){
    
    if(argc < 5) Usage();

    args.collen = atoi(argv[1]);
    args.ptbin = atoi(argv[2]);
    args.jetparam = atof(argv[3]);
    args.mode = argv[4];
    args.nEvents = -1;
    args.verbose = 0;
    if(argc > 5)  args.nEvents = atoi(argv[5]);
    if(argc > 6)  args.verbose = atoi(argv[6]);
    
    if(args.collen != 200 && args.collen != 2760){
        cout << "Invalid collen: " << args.collen << endl;
        Usage();
    }
    if(args.ptbin < 0 || args.ptbin > 24){
        cout << "Invalid ptbin: " << args.ptbin << endl;
        Usage();
    }
    if(args.jetparam < 0.2 || args.jetparam > 0.7){
        cout << "Invalid jetparam: " << args.jetparam << endl;
        Usage();
    }
    if(args.mode != "PP" && args.mode != "FullJets" && args.mode != "Match" && args.mode != "Missed" && args.mode != "Fake"){
        cout << "Invalid mode: " << args.mode << endl;
        Usage();
    }
    if(args.nEvents < -1){
        cout << "Invalid nEvents: " << args.nEvents << endl;
        Usage();
    }
    if(args.verbose != 0 && args.verbose != 1){
        cout << "Invalid verbose: " << args.verbose << endl;
        Usage();
    }

}

void Configure(struct Args args, struct Config &conf, const TString OUTPUT_DIR, const TString PYTHIA_DIR, const TString TENNGEN_DIR){

    
    conf.output_base_directory = GetOutFilePath(OUTPUT_DIR, args.collen, args.jetparam, args.mode);
    conf.pythia_base_directory = PYTHIA_DIR;
    conf.tenngen_base_directory = TENNGEN_DIR;
    conf.pythia_filename = Form("%s/signal/root-files/%dGeV_PP_ptbin%d.root", PYTHIA_DIR.Data(), args.collen, args.ptbin); 
    
    conf.n_cent_bins = args.collen == 200 ? 3 : 5;
    conf.collsys = args.collen == 200 ? "AuAu" : "PbPb";
    conf.cent_dir_names = args.collen == 200 ? std::vector<TString>({"0-10", "10-20", "20-40"}) : std::vector<TString>({"0-5", "5-10", "10-20", "20-30", "30-40"});
    conf.cent_event_start_fractions = args.collen == 200 ? std::vector<Float_t>({0.0, 0.25, 0.5}) : std::vector<Float_t>({0.0, 0.125, 0.25, 0.5, 0.75});
    conf.cent_event_end_fractions = args.collen == 200 ? std::vector<Float_t>({0.25, 0.5, 1.0}) : std::vector<Float_t>({0.125, 0.25, 0.5, 0.75, 1.0});

    conf.merged_output_filename =Form("%s%dGeV_%s_ptbin%d_cent0to40_R0%.0f.root", conf.output_base_directory.Data(), args.collen, conf.collsys.Data(), args.ptbin, args.jetparam*10);
    conf.pythia_jets_output_filename = Form("%s%dGeV_%s_ptbin%d_R0%.0f.root", conf.output_base_directory.Data(), args.collen, conf.collsys.Data(), args.ptbin, args.jetparam*10);
    conf.tenngen_filenames.clear();
    conf.output_filenames.clear();

    for (Int_t icent = 0; icent < conf.n_cent_bins; icent++){
        conf.tenngen_filenames.push_back(Form("%s/%s/%s/%dGeV_%s.root", conf.tenngen_base_directory.Data(), conf.collsys.Data(), conf.cent_dir_names.at(icent).Data(), args.collen, conf.collsys.Data()));
        conf.output_filenames.push_back(Form("%s%dGeV_%s_ptbin%d_cent%d_R0%.0f.root", conf.output_base_directory.Data(), args.collen, conf.collsys.Data(), args.ptbin, icent, args.jetparam*10));
    }

    if (args.verbose){
        cout << "Output Directory: " << conf.output_base_directory.Data() << endl;
        cout << "Pythia Directory: " << conf.pythia_base_directory.Data() << endl;
        cout << "Tenngen Directory: " << conf.tenngen_base_directory.Data() << endl;
        cout << "Pythia Filename: " << conf.pythia_filename.Data() << endl;
        cout << "Pythia Jets Output Filename: " << conf.pythia_jets_output_filename.Data() << endl;
        cout << "Merged Output Filename: " << conf.merged_output_filename.Data() << endl;
        for (Int_t icent = 0; icent < conf.n_cent_bins; icent++){
            cout << "Tenngen Filename: " << conf.tenngen_filenames.at(icent).Data() << endl;
            cout << "Output Filename: " << conf.output_filenames.at(icent).Data() << endl;
        }

        cout << "n_cent_bins: " << conf.n_cent_bins << endl;
        cout << "collsys: " << conf.collsys.Data() << endl;
        cout << "cent_dir_names: " << endl;
        for (Int_t icent = 0; icent < conf.n_cent_bins; icent++){
            cout << conf.cent_dir_names.at(icent).Data() << endl;
        }
       
    }

}

void SetJetFindingParameters(struct Args args, struct Config conf, struct JetFindingParameters &jetparams){

    jetparams.jetparam = args.jetparam;

    TFile pythia_file(conf.pythia_filename.Data());
    if(!pythia_file.IsOpen()){ cout << "PYTHIA file not valid" << endl; exit(1); }
    TTree *pythia_tree = (TTree*)pythia_file.Get("tree");
    if(!pythia_tree){ cout << "PYTHIA tree not valid" << endl; exit(1); }

    jetparams.total_pythia_events = pythia_tree->GetEntries();
    pythia_file.Close();
    jetparams.n_pythia_events = args.nEvents == -1 ? jetparams.total_pythia_events : args.nEvents;

    for(Int_t icent = 0; icent < conf.n_cent_bins; icent++){
        jetparams.n_mixed_events.push_back( (Int_t) (jetparams.n_pythia_events * (conf.cent_event_end_fractions.at(icent) - conf.cent_event_start_fractions.at(icent))) );
        jetparams.start_mixed_events.push_back( (Int_t) (jetparams.n_pythia_events * conf.cent_event_start_fractions.at(icent)) );
    }

    jetparams.particle_max_eta = 0.9;
    jetparams.jet_pt_min = 10.0;
    
}

void Debug(struct Args args, struct Config config){

    cout << "=====================================" << endl;
    cout << "Arguments:" << endl;
    cout << "nEvents: " << args.nEvents << endl;
    cout << "mode: " << args.mode.Data() << endl;
    cout << "collen: " << args.collen << endl;
    cout << "ptbin: " << args.ptbin << endl;
    cout << "jetparam: " << args.jetparam << endl;
    cout << "verbose: " << args.verbose << endl;
    cout << "=====================================" << endl;
    cout << "Config:" << endl;
    cout << "Output directory: " << config.output_base_directory.Data() << endl;
    cout << "Particle directory: " << config.pythia_base_directory.Data() << endl;
    cout << "TennGen directory: " << config.tenngen_base_directory.Data() << endl;
    cout << "n_cent_bins: " << config.n_cent_bins << endl;
    cout << "collsys: " << config.collsys.Data() << endl;
    cout << "cent_dir_names: ";
    for (auto &name : config.cent_dir_names) cout << name << " ";
    cout << endl;
    cout << "cent_event_fractions: ";
    for (auto &frac : config.cent_event_start_fractions) cout << frac << " ";
    cout << endl;
    for (auto &frac : config.cent_event_end_fractions) cout << frac << " ";
    cout << endl;
    cout << "pythia_filename: " << config.pythia_filename.Data() << endl;
    cout << "tenngen_filenames: ";
    for (auto &name : config.tenngen_filenames) cout << name << " ";
    cout << endl;
    cout << "output_tfile_names: ";
    for (auto &name : config.output_filenames) cout << name << " ";
    cout << endl;


}

Int_t main(Int_t argc, Char_t** argv){

    cout << "Starting Jet Analysis" << endl;
    // Initialize arguments
    Args args;
    InitJetAnalysis(argc, argv, args);

    // configure for different modes
    // data paths
    const TString output_directory = "/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/jet-trees/root-files/";
    const TString particle_directory = "/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees";
    const TString tenngen_directory = "/lustre/isaac/scratch/tmengel/jet-background-subtraction/simulation/particle-trees/background/root-files";
    Config config;
    Configure(args, config, output_directory, particle_directory, tenngen_directory);

    // Debug(args, config);

    // Set File parameters
    JetFindingParameters params;
    SetJetFindingParameters(args, config, params);

    if( args.mode == "PP" ) PythiaJets(config.pythia_filename, config.pythia_jets_output_filename, params);
    else
    {

        for (Int_t icent = 0; icent < config.n_cent_bins; icent++)
        {
            if( args.mode == "FullJets" ) Jets(config.pythia_filename, config.tenngen_filenames.at(icent), config.output_filenames.at(icent), icent,  params);
            else if( args.mode == "Match" ) JetMatch(config.pythia_filename, config.tenngen_filenames.at(icent), config.output_filenames.at(icent), icent, params);
            else if( args.mode == "Missed" ) MissedPythiaJets(config.pythia_filename, config.tenngen_filenames.at(icent), config.output_filenames.at(icent), icent, params);
            else if( args.mode == "Fake" ) UnmatchedJets(config.pythia_filename, config.tenngen_filenames.at(icent), config.output_filenames.at(icent), icent, params);
        }

        TChain *out_chain = new TChain("tree");
        for (auto &name : config.output_filenames) out_chain->Add(name.Data());
        out_chain->Merge(config.merged_output_filename.Data());
        for (auto &name : config.output_filenames) gSystem->Unlink(name.Data());
        out_chain->Delete();
    }
    return 0;


}
