#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TSystem.h>

#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace std;
using namespace ROOT;

struct Network_t{

    std::vector<TH2D*> w;
    std::vector<TH1D*> b;
    std::vector<TString> activations;
    std::vector<TString> inputs;
    std::vector<TString> targets;
    std::vector<TString> input_types;
    std::vector<TString> target_types;
    Int_t n_layers;
    Int_t n_inputs;
    Int_t n_targets;
    TString model_type;
    TString network_file;

};

void Compile(TString network_file, Network_t &network){

    network.network_file = network_file;
    TFile *fin = new TFile(network_file, "READ");
    if(!fin->IsOpen()){
        std::cout << "Error: could not open network file " << network_file << std::endl;
        exit(1);
    }


    TTree *arch_tree = (TTree*)fin->Get("arch");
    TTree *activation_tree = (TTree*)fin->Get("activations");
    TTree *input_tree = (TTree*)fin->Get("inputs");
    TTree *target_tree = (TTree*)fin->Get("targets");

    // architecture
    Int_t n_layers, n_inputs, n_targets;
    Char_t model_type_char[100];
    arch_tree->SetBranchAddress("n_layers", &n_layers);
    arch_tree->SetBranchAddress("n_inputs", &n_inputs);
    arch_tree->SetBranchAddress("n_targets", &n_targets);
    arch_tree->SetBranchAddress("model_type", &model_type_char);

    // activations
    Char_t acts_char[100];
    activation_tree->SetBranchAddress("activations", &acts_char);

    // inputs
    Char_t ins_char[100];
    Char_t in_types_char[100];
    input_tree->SetBranchAddress("inputs", &ins_char);
    input_tree->SetBranchAddress("input_types", &in_types_char);

    // targets
    Char_t target_char[100];
    Char_t target_types_char[100];
    target_tree->SetBranchAddress("targets", &target_char);
    target_tree->SetBranchAddress("target_types", &target_types_char);
   
    // read achitecture
    arch_tree->GetEntry(0);
    network.n_layers = n_layers;
    network.n_inputs = n_inputs;
    network.n_targets = n_targets;
    network.model_type = Form("%s", model_type_char);
    // read activations and weights
    for (int i = 0; i < n_layers; i++)
    {
        activation_tree->GetEntry(i);
        network.activations.push_back(Form("%s", acts_char));
        TH2D* w = (TH2D*)fin->Get(Form("weight%d", i));
        TH1D* b = (TH1D*)fin->Get(Form("bias%d", i));
        network.w.push_back(w);
        network.b.push_back(b);
    }
    // read inputs
    for (int i = 0; i < n_inputs; i++)
    {
        input_tree->GetEntry(i);
        network.inputs.push_back(Form("%s", ins_char));
        network.input_types.push_back(Form("%s", in_types_char));
    }
    // read targets
    // for (int i = 0; i < n_targets; i++)
    // {
        // target_tree->GetEntry(i);
        // network.targets.push_back(Form("%s", target_char));
        // network.target_types.push_back(Form("%s", target_types_char));
        // network.targets.push_back(Form("jet_pt_truth"));
        // network.target_types.push_back(Form("Float_t"));
    // }
    network.targets.push_back(Form("jet_pt_truth"));
    network.target_types.push_back(Form("Float_t"));
    // fin->Close();

}

void Summary(Network_t network){
    
    // cout for debugging
    std::cout << "====================" << std::endl;
    std::cout << "Network file: " << network.network_file.Data() << std::endl;
    std::cout << "Network architecture:" << std::endl;
    std::cout << "\tNumber of layers: " << network.n_layers << std::endl;
    std::cout << "\tNumber of inputs: " << network.n_inputs << std::endl;
    std::cout << "\tNumber of targets: " << network.n_targets << std::endl;
    std::cout << "\tModel type: " << network.model_type.Data() << std::endl;
    std::cout << "====================" << std::endl;
    std::cout << "\tActivations: ";
    for (int i = 0; i < network.activations.size(); i++)
    {
        std::cout << network.activations[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "\tInputs: ";
    for (int i = 0; i < network.inputs.size(); i++)
    {
        std::cout << network.inputs[i] << " (" << network.input_types[i] << ") ";
    }
    std::cout << std::endl;
    std::cout << "\tTargets: ";
    for (int i = 0; i < network.targets.size(); i++)
    {
        std::cout << network.targets[i] << " (" << network.target_types[i] << ") "; 
    }
    std::cout << std::endl;

    TCanvas *c = new TCanvas("c", "c", 1600, 500);
    // divide canvas into subpads for each layer
    c->Divide(network.w.size(), 2);
    // plot weights and biases
    for (int i = 0; i < network.w.size(); i++)
    {
        c->cd(i+1);
        network.w[i]->SetStats(0);
        network.w[i]->SetTitle(Form("Weights %d", i));
        network.w[i]->GetXaxis()->SetTitle("Input Dim");
        network.w[i]->GetYaxis()->SetTitle("Output Dim");
         network.w[i]->Draw("colz");
        c->Update();
        c->cd(i+1+network.w.size());
        network.b[i]->SetStats(0);
        network.b[i]->SetTitle(Form("Biases %d", i));
        network.b[i]->GetXaxis()->SetTitle("Output Dim");
        network.b[i]->GetYaxis()->SetTitle("Bias");
        network.b[i]->Draw();
        c->Update();

    }

    TString summary_name = network.network_file;
    summary_name.ReplaceAll(".root", "_summary.png");
    // remove path from filename
    summary_name.Replace(0, summary_name.Last('/')+1, "");
    TString dirpng = "weight-pngs/";
    TString png_name = Form("%s%s", dirpng.Data(), summary_name.Data());
    c->SaveAs(png_name.Data());


}

std::vector<Double_t> Layer(std::vector<Double_t> x, TH2D* w, TH1D* b, TString activation)
{
    std::vector<Double_t> y;
    for (int i = 0; i < w->GetNbinsY(); i++)
    {
        Double_t sum = 0;
        for (int j = 0; j < w->GetNbinsX(); j++)
        {
            sum += w->GetBinContent(j+1, i+1)*x[j];
        }
        sum += b->GetBinContent(i+1);
        if (activation == "relu")
        {
            if (sum < 0) sum = 0.0;
        }
        else if (activation == "linear")
        {
            // do nothing
        }
        else
        {
            std::cout << "Error: activation function " << activation.Data() << " not recognized" << std::endl;
            exit(1);
        }
        y.push_back(sum);
    }
    return y;
   
}

Float_t Predict(std::vector<Double_t> x, Network_t network)
{
    Float_t y;
    std::vector<Double_t> x_layer = x;
    for (int i = 0; i < network.w.size(); i++)
    {
        x_layer = Layer(x_layer, network.w[i], network.b[i], network.activations[i]);
    }
    //check if output is a single value
    if(x_layer.size() != 1){
        std::cout << "Error: output is not a single value" << std::endl;
        return -1;
    }
    y = x_layer[0]*1.0;
    return y;
}

void Evaluate(Network_t network, TString infile, TString outfile){

    TFile *f = new TFile(infile, "READ");
    if(!f->IsOpen()){
        std::cout << "Error: could not open file " << infile.Data() << std::endl;
        exit(1);
    }
    TTree *tin = (TTree*)f->Get("tree");

    // get floats and ints for passive variables
    std::vector<TString> used_variables;
    for (int i = 0; i < network.inputs.size(); i++)
    {
        used_variables.push_back(network.inputs[i]);
    }
    for (int i = 0; i < network.targets.size(); i++)
    {
        used_variables.push_back(network.targets[i]);
    }

    Int_t n_passive_ints = 0;
    Int_t n_passive_floats = 0;
    std::vector<TString> passive_ints_names;
    std::vector<TString> passive_float_names;

    Int_t found_target = 0;
    TObjArray *branches = tin->GetListOfBranches();
    for (Int_t i = 0; i < branches->GetEntries(); i++){
        TBranch *branch = (TBranch*)branches->At(i);
        TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
        TString branch_name = branch->GetName();
        TString leaf_type = leaf->GetTypeName();
        // check if variable is used
        Int_t save = 1;
        for (Int_t j = 0; j < used_variables.size(); j++){
            if (branch_name == used_variables[j]){
                save = 0;
            }
        }
        for(Int_t j = 0; j < network.targets.size(); j++){
            if(branch_name == network.targets[j]){
                found_target = 1;
            }
        }
        // save passive variables
        if(save){
            if (leaf_type == "Int_t"){
                passive_ints_names.push_back(branch_name);
                n_passive_ints++;
            }
            else if (leaf_type == "Float_t"){
                passive_float_names.push_back(branch_name);
                n_passive_floats++;
            }
        }
    }


    // get number of ints and floats 
    Int_t n_input_ints = 0;
    Int_t n_input_floats = 0;
    Int_t int_input_index[network.inputs.size()];
    for (int i = 0; i < network.input_types.size(); i++)
    {
        if (network.input_types[i] == "Int_t")
        {
            n_input_ints++;
            int_input_index[i] = 1;
        }
        else
        {
            n_input_floats++;
            int_input_index[i] = 0;
        }
    }

    // variable arrays
    Int_t int_inputs[n_input_ints];
    Float_t float_inputs[n_input_floats];
    Float_t float_targets[network.targets.size()];
    Int_t int_passive[n_passive_ints];
    Float_t float_passive[n_passive_floats];
    
    // set branch addresses for passive variables
    for (int i = 0; i < passive_ints_names.size(); i++)
    {
        tin->SetBranchAddress(passive_ints_names[i], &int_passive[i]);
    }
    for (int i = 0; i < passive_float_names.size(); i++)
    {
        tin->SetBranchAddress(passive_float_names[i], &float_passive[i]);
    }

    // set branch address for targets
    if(found_target){
        
        for (int i = 0; i < network.targets.size(); i++)
        {
            tin->SetBranchAddress(network.targets[i], &float_targets[i]);
        }
    }

    // set branch address for inputs
    for(int i = 0; i < network.inputs.size(); i++)
    {
        if (int_input_index[i] == 1)
        {
            tin->SetBranchAddress(network.inputs[i], &int_inputs[i]);
        }
        else
        {
            tin->SetBranchAddress(network.inputs[i], &float_inputs[i]);
        }
    }

    // create output tmp file
    TString tmpfile = outfile;
    tmpfile.ReplaceAll(".root", "_tmp.root");
    TFile *fout = new TFile(tmpfile, "RECREATE");
    TTree *tout = new TTree("tree", "tree");

    // create output branches for passive variables
    for (int i = 0; i < passive_ints_names.size(); i++)
    {
        tout->Branch(passive_ints_names[i], &int_passive[i], Form("%s/I", passive_ints_names[i].Data()));
    }
    for (int i = 0; i < passive_float_names.size(); i++)
    {
        tout->Branch(passive_float_names[i], &float_passive[i], Form("%s/F", passive_float_names[i].Data()));
    }

    // create output branches for targets
    if(found_target){
        for (int i = 0; i < network.targets.size(); i++)
        {
            tout->Branch(network.targets[i], &float_targets[i], Form("%s/F", network.targets[i].Data()));
        }
    }

    // create output branches for inputs
    for(int i = 0; i < network.inputs.size(); i++)
    {
        if (int_input_index[i] == 1)
        {
            tout->Branch(network.inputs[i], &int_inputs[i], Form("%s/I", network.inputs[i].Data()));
        }
        else
        {
            tout->Branch(network.inputs[i], &float_inputs[i], Form("%s/F", network.inputs[i].Data()));
        }
    }

    Float_t prediction;
    TString prediction_varibale_name = Form("jet_pt_%s_corrected", network.model_type.Data());
    tout->Branch(prediction_varibale_name, &prediction, Form("%s/F", prediction_varibale_name.Data()));

    // loop over events in input file
    Int_t n_entries = tin->GetEntries();
    for (Int_t i = 0; i < n_entries; i++)
    {
        tin->GetEntry(i);
        std::vector<Double_t> x;
        for (int j = 0; j < network.inputs.size(); j++)
        {
            if (int_input_index[j] == 1)
            {
                x.push_back(int_inputs[j]*1.0);
            }
            else
            {
                x.push_back(float_inputs[j]*1.0);
            }
        }

        prediction = Predict(x, network);
        tout->Fill();
    }

    // write output file
    fout->Write();
    fout->Close();

    // close input file
    f->Close();

    // cp tmp file to input file
    TString command = Form("cp %s %s", tmpfile.Data(), outfile.Data());
    gSystem->Exec(command);

    // remove tmp file
    command = Form("rm %s", tmpfile.Data());
    gSystem->Exec(command);

}

// main function takes in network root file and input root file
Int_t main(Int_t argc, Char_t** argv){

    if(argc < 4){
        std::cout << "Usage: " << argv[0] << " <network.root> <input.root> <output.root>" << std::endl;
        return 1;
    }

    TString network_file = argv[1];
    TString input_file = argv[2];
    TString output_file = argv[3];

    Network_t network;
    Compile(network_file, network);
    Summary(network);

    Evaluate(network, input_file, output_file);
    
    return 0;

    
}


