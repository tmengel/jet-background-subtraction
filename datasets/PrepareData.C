// Author : Tanner Mengel
/*
    Reads in all pt-bins in input directory and creates a single output file with all pt-bins.
    Automatically reads all directories in input directory and creates a single output file within the input director for each sub-directory.
*/

#include <iostream>
#include <vector>

#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TTree.h>
#include <TChain.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TLegend.h>

using namespace std;
using namespace ROOT;
 
struct FileLists_t{
    TString base_directory;
    vector<TString> sub_directories;
    vector<vector<TString>> file_lists;
    vector<TString> output_filename;
    vector<TString> tree_names;
};

struct Args_t{
    TString base_directory;
    TString output_directory;
    Int_t verbose;
};

void Usage(Int_t argc, Char_t** argv)
{
    cout << "Usage: " << argv[0] << " -d <input base directory> -o <output base directory> -v <optional verbose>" << endl;
    cout << "- input base directory: directory containing all sub-directories with root files to be merged." << endl;
    cout << "- output base directory: directory to store output root files." << endl;
    cout << "Example: " << argv[0] << " -d /home/tmengel/analysis/jet-background-subtraction/simulation/jet-trees/root-files -o /home/tmengel/analysis/jet-background-subtraction/simulation/jet-trees/root-files-merged -v 0" << endl;
    cout << "--> This will read all files in the directory and sub-directories and create a single output file for each sub-directory." << endl;
    cout << "--> Produces output root files in <output base directory>/<sub-directory>/<output filename>.root" << endl;
    cout << "--> Output file name is the same as root files in sub-directory, minus the pt bin" << endl;
}

void ArgParse(Int_t argc, Char_t** argv, Args_t& args){
    if (argc < 5) {
        Usage(argc, argv);
        exit(1);
    }
    args.verbose = 0;
    for (Int_t i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-d") == 0) {
            i++;
            if (i >= argc) {
                Usage(argc, argv);
                exit(1);
            }
            args.base_directory = argv[i];
            if (args.base_directory.EndsWith("/")) {
                args.base_directory.Remove(args.base_directory.Length()-1);
            }
            // check if base directory exists
            if(gSystem->AccessPathName(args.base_directory)) {
                cout << "Base directory does not exist. Exiting" << endl;
                exit(1);
            }

        }
        else if (strcmp(argv[i], "-o") == 0) {
            i++;
            if (i >= argc) {
                Usage(argc, argv);
                exit(1);
            }
            args.output_directory = argv[i];
            if (args.output_directory.EndsWith("/")) {
                args.output_directory.Remove(args.output_directory.Length()-1);
            }
            // check if output directory exists
            if(gSystem->AccessPathName(args.output_directory)) {
                cout << "Output directory does not exist. Exiting" << endl;
                exit(1);
            }
        }
        else if (strcmp(argv[i], "-v") == 0) {
            i++;
            if (i >= argc) {
                Usage(argc, argv);
                exit(1);
            }
            args.verbose = atoi(argv[i]);
        }
        else {
            Usage(argc, argv);
            exit(1);
        }
    }

}

std::vector <TString> ListDirectory(TString directory)
{
    // List all files in directory
    void* dirp = gSystem->OpenDirectory(directory);
    const char* dir_file;
    vector <TString> file_list;
    while ((dir_file = gSystem->GetDirEntry(dirp))) {
        TString file_name = dir_file;
        if(file_name == "." || file_name == ".." || file_name == ".DS_Store") continue;
        file_list.push_back(Form("%s/%s", directory.Data(), file_name.Data()));
    }
    gSystem->FreeDirectory(dirp);
    return file_list;
}

void Init(struct FileLists_t& file_lists, struct Args_t args)
{
    // Initialize base directory
    file_lists.base_directory = args.base_directory;
    
    // Each base directory contains sub-directories with either more sub-directories or root files
    // Dig through base directory and find all sub-directories until we find root files

    TString current_directory = file_lists.base_directory;
    vector<TString> current_sub_directories = ListDirectory(current_directory);
    
    // find all sub-directories with root files
    while (current_sub_directories.size() > 0) {

        for (auto sub_directory : current_sub_directories) {
            // check if sub-directory contains root files
            vector<TString> sub_directory_files = ListDirectory(sub_directory);
            Bool_t has_root_files = false;
            for (auto file : sub_directory_files) {
                if (file.EndsWith(".root")) { // check for root files
                    has_root_files = true;
                }
                else if (!gSystem->AccessPathName(file)) { // check if file is a directory
                    current_sub_directories.push_back(file);
                }

            }

            // if sub-directory contains root files, add to list of sub-directories with root files to be merged
            if (has_root_files) {
                file_lists.sub_directories.push_back(sub_directory);
            }
            
            // remove sub-directory from list of sub-directories to check
            current_sub_directories.erase(remove(current_sub_directories.begin(), current_sub_directories.end(), sub_directory), current_sub_directories.end());
        }
        
    }

    // Get list of files in each sub-directory
    for (auto sub_directory : file_lists.sub_directories) {
        vector<TString> file_list;
        vector<TString> sub_directory_files = ListDirectory(sub_directory);
        TString output_filename = "";
        for (auto file : sub_directory_files) {
            if (file.EndsWith(".root")) { // check for root files
                file_list.push_back(file);
            }
            if (output_filename == "") {
                // get file base name
                TString filename = file;
                // remove .root and path
                filename.Remove(filename.Length()-5);
                filename.Remove(0, filename.Last('/')+1);
                // files have the form <energy>_<species>_ptbin<ptbin>_<radius> 
                // get locations of underscores
                vector<Int_t> underscore_locations;
                for (Int_t i = 0; i < filename.Length(); i++) {
                    if (filename[i] == '_') {
                        underscore_locations.push_back(i);
                    }
                }
                // remove ptbin<ptbin>
                filename.Remove(underscore_locations[1], underscore_locations[2]-underscore_locations[1]);
                // remove species
                filename.Remove(underscore_locations[0], underscore_locations[1]-underscore_locations[0]);
                // Get Energy
                TString energy = filename(0, underscore_locations[0]);
                // Get radius
                TString radius = filename(filename.Last('_')+1, filename.Length()-filename.Last('_')-1);
                // make output filename
                TString current_sub_directory = sub_directory;
                current_sub_directory.Remove(0, current_sub_directory.Last('/')+1);

                TString output_directory = args.output_directory;
                output_directory += Form("/%s", energy.Data());
                if (gSystem->AccessPathName(output_directory)) gSystem->mkdir(output_directory);
                output_directory += Form("/%s", radius.Data());
                if (gSystem->AccessPathName(output_directory)) gSystem->mkdir(output_directory);
                output_filename = Form("%s/%s_%s_%s.root", output_directory.Data(), energy.Data(), current_sub_directory.Data(), radius.Data());
            }
        }
        file_lists.file_lists.push_back(file_list);
        file_lists.output_filename.push_back(output_filename);
    }

    // Make sure same tree name is used in all files
    for(auto sub_directory : file_lists.sub_directories) {
        vector<TString> file_list = file_lists.file_lists[distance(file_lists.sub_directories.begin(), find(file_lists.sub_directories.begin(), file_lists.sub_directories.end(), sub_directory))];
        TString tree_name = "";
        for (auto file : file_list) {
            TFile* f = TFile::Open(file);
            TList* key_list = f->GetListOfKeys();
            for (auto key : *key_list) {
                TString key_name = key->GetName();
                if (key_name == "tree") {
                    if (tree_name == "") {
                        tree_name = key_name;
                    }
                    else if (tree_name != key_name) {
                        cout << "Error: different tree names found in files in sub-directory: " << sub_directory << endl;
                        exit(1);
                    }
                }
            }

            f->Close();
        }

        file_lists.tree_names.push_back(tree_name);
    }

    // Print out list of files to be merged
    if(args.verbose) {
        for (Int_t i = 0; i < file_lists.sub_directories.size(); i++) {
        cout << "Sub-directory: " << file_lists.sub_directories[i] << endl;
        cout << "Output file: " << file_lists.output_filename[i] << endl;
        cout << "Tree name: " << file_lists.tree_names[i] << endl;
        for (auto file : file_lists.file_lists[i]) {
            cout << file << endl;
        }
        cout << endl;
    }
    }
}

void Merge(struct FileLists_t file_lists){

    // Merge rootfiles in each sub-directory
    for (auto sub_directory : file_lists.sub_directories) {
        // get list of files in sub-directory
        vector<TString> file_list = file_lists.file_lists[distance(file_lists.sub_directories.begin(), find(file_lists.sub_directories.begin(), file_lists.sub_directories.end(), sub_directory))];
        // get output filename
        TString output_filename= file_lists.output_filename[distance(file_lists.sub_directories.begin(), find(file_lists.sub_directories.begin(), file_lists.sub_directories.end(), sub_directory))];
        // get tree name
        TString tree_name = file_lists.tree_names[distance(file_lists.sub_directories.begin(), find(file_lists.sub_directories.begin(), file_lists.sub_directories.end(), sub_directory))];

        cout << "Merging files in sub-directory: " << sub_directory << endl;
        cout << "Output file: " << output_filename << endl;
        cout << "Tree name: " << tree_name << endl;
        for (auto file : file_list) {
            cout << file << endl;
        }

        // create TChain
        TChain* chain = new TChain(tree_name.Data());
        for (auto file : file_list) {
            chain->Add(file);
        }

        // merge files
        chain->Merge(output_filename.Data());
        //cout << "Merged to file: " << output_filename << endl;
        
        // delete TChain
        delete chain;
    }

}

Int_t main(Int_t argc, Char_t** argv)
{
    // Parse arguments
    Args_t args;
    ArgParse(argc, argv, args);

    // Initialize variables
    FileLists_t file_lists;
    Init(file_lists, args);

    // Merge files
    Merge(file_lists);

}


    