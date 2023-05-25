#!/usr/bin/env python3
import warnings
import argparse

try:
    import pandas as pd
except:
    warnings.warn("Pandas not installed, please install it to use this function")
    

def ConvertTFileToHdf5(datafile, outputfile):
    import ROOT
    ROOT.ROOT.EnableImplicitMT()
    nThreads = ROOT.ROOT.GetThreadPoolSize()
    print("Using", nThreads, "threads")
    
    # get the root file
    rdf = ROOT.ROOT.RDataFrame("tree", datafile)
    variables = rdf.GetColumnNames()
    variables_for_strings = [str(x) for x in variables]
    
    # convert to pandas
    pdf = pd.DataFrame(rdf.AsNumpy(columns=variables_for_strings))
    ROOT.ROOT.DisableImplicitMT()
    
    # save to hdf5
    pdf.to_hdf(outputfile, key='df', mode='w')
    

def SampleHDF5(datafile, nevents=10000):
    pdf = pd.read_hdf(datafile)
    return pdf.sample(nevents).copy()
   
   
def FitPySrModel(datafile, outputfile, target, features):
    
    try:
        import pysr
        from pysr import PySRRegressor
    except:
        warnings.warn("PySR not installed, please install it to use this function")
        return None, None

    pdf = SampleHDF5(datafile,nevents=10000)
    fileout = outputfile.replace(".csv","")
    model = PySRRegressor(
            model_selection="best",  # Result is mix of simplicity+accuracy
            niterations=40, #40
            populations=30, #20
            maxdepth=15,
            binary_operators=["plus", "sub", "mult", "div"],
            unary_operators=["square"],
            select_k_features=3, # Select 5 features
            loss="L2DistLoss()", #mean_squared_error
            equation_file=f"{fileout}.csv",
            verbosity=0
        )
    
    pysr_target = "jet_pt_correction"
    pdf[pysr_target] = pdf[target] - pdf["jet_pt_raw"]
    
    X = pdf[features].values
    Y = pdf[pysr_target].values

    model.fit(X, Y,variable_names=features)
    return model.equation_file_ , model.get_best()["sympy_format"]


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--datafile", required=True, help="Path to data file", default=None)
    parser.add_argument("-o", "--outputfile", required=True, help="Path to output file", default=None)
    parser.add_argument("-m", "--model", required=True, help="Model to use", default="SNN")
    args = parser.parse_args()
    
    
    
    
    
    target = None 
    features = []
    pysr_output = f"equations-files/{args.outputfile}.csv"
    
    if args.model == "SNN":
        target = "jet_pt_snn_corrected"
        features = ["jet_pt_raw","jet_nparts","jet_area","jet_angularity","jet_track_pt_0"]
    
    elif args.model == "DNN":
        target = "jet_pt_dnn_corrected"
        # features = ["jet_pt_raw","jet_nparts","jet_area","jet_angularity","jet_track_pt_0"]
        features = ["jet_pt_raw","jet_nparts","jet_area","jet_angularity","jet_track_pt_0",
                        "jet_track_pt_1","jet_track_pt_2","jet_track_pt_3","jet_track_pt_4",
                        "jet_track_pt_5","jet_track_pt_6","jet_track_pt_7"]
    
    
    mod, mod_symbolic = FitPySrModel(datafile=args.datafile, outputfile=pysr_output, target=target, features=features)