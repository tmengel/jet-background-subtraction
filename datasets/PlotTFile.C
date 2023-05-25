// Author : Tanner Mengel
/*
    Creates plots for all variables in a TFile
    saves plots in a new TFile
*/

#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>

#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>

#include <TH1D.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace std;
using namespace ROOT;

void Usage(Int_t argc, Char_t** argv)
{
    cout << "Usage: " << argv[0] << " <input file>" << endl;
    cout << "Example: " << argv[0] << " /path/to/input/file.root" << endl;
    cout << "Program will create a new file with the same name as the input file, but with the extension _plots.root" << endl;
}

struct Args{
    TString input_file;
};

struct TreeData
{
    TString tree_name;
    std::vector<TString> tree_branches;
    std::vector<TString> tree_branches_types;
    std::vector<Int_t> branch_is_int;
    Int_t n_entries;
    Int_t n_floats, n_ints;

};

Args GetArgs(Int_t argc, Char_t** argv)
{
    Args args;
    if (argc != 2)
    {
        Usage(argc, argv);
        exit(1);
    }
    args.input_file = argv[1];
    return args;
}

void SetTreeData(Args args, TreeData &tree_data)
{
    TFile *fin = new TFile(args.input_file, "READ");
    TTree *t = (TTree*)fin->Get("tree");
    tree_data.tree_name = Form("%s", t->GetName());
    tree_data.n_entries = t->GetEntries();
    TObjArray *branches = t->GetListOfBranches();
    for (Int_t i = 0; i < branches->GetEntries(); i++)
    {
        TBranch *branch = (TBranch*)branches->At(i);
        TLeaf *leaf = (TLeaf*)branch->GetListOfLeaves()->At(0);
        TString branch_name = branch->GetName();
        TString leaf_type = leaf->GetTypeName();

        tree_data.tree_branches.push_back(branch_name);
        tree_data.tree_branches_types.push_back(leaf_type);
        if (leaf_type == "Int_t")
        {
            tree_data.branch_is_int.push_back(1);
            tree_data.n_ints++;
        }
        else if (leaf_type == "Float_t")
        {
            tree_data.branch_is_int.push_back(0);
            tree_data.n_floats++;
        }
        else
        {
            cout << "Error: branch " << branch_name << " is neither Int_t nor Float_t" << endl;
            exit(1);
        }
    }
    fin->Close();
}


Int_t main(Int_t argc, Char_t** argv)
{
    Args args = GetArgs(argc, argv);
    TreeData tree_data;
    SetTreeData(args, tree_data);

    TFile *fin = new TFile(args.input_file, "READ");
    TTree *t = (TTree*)fin->Get("tree");

    TString output_file = args.input_file;
    output_file.ReplaceAll(".root", "_plots.root");
    TFile *fout = new TFile(output_file, "RECREATE");

    for (Int_t i = 0; i < tree_data.tree_branches.size(); i++)
    {
        TString branch_name = tree_data.tree_branches[i];
        TString branch_type = tree_data.tree_branches_types[i];
        Int_t is_int = tree_data.branch_is_int[i];

        if (is_int)
        {
            Int_t min = t->GetMinimum(branch_name);
            Int_t max = t->GetMaximum(branch_name);
            Int_t n_bins = max - min + 1;
            TH1D *h = new TH1D(Form("%s", branch_name.Data()), Form("%s", branch_name.Data()), n_bins, min, max);
            t->Draw(Form("%s>>%s", branch_name.Data(), branch_name.Data()));
            h->Write();
        }
        else
        {
            Float_t min = t->GetMinimum(branch_name);
            Float_t max = t->GetMaximum(branch_name);
            Int_t n_bins = 100;
            TH1D *h = new TH1D(Form("%s", branch_name.Data()), Form("%s", branch_name.Data()), n_bins, min, max);
            t->Draw(Form("%s>>%s", branch_name.Data(), branch_name.Data()));
            h->Write();
        }
    }

    fout->Close();
    fin->Close();
    return 0;
}

