# root utils for loading and writing ROOT files
# reading trees and histograms

from __future__ import absolute_import, division, print_function

import warnings

__all__ = ['LoadTFile', 'WriteTFile', 'ReadTreeKeys', 'GetTFileConfig']

try:
    import ROOT
    from ROOT import TFile
    from ROOT import RDataFrame as rd
    import pandas as pd
    
except ImportError as e:
    warnings.warn('Could not import - ' + str(e))


def LoadTFile(path, tree="tree", vars=None):
    ROOT.ROOT.EnableImplicitMT()
    if vars is not None:
        rdf = rd(tree, path).AsNumpy(vars)
    else:
        rdf = rd(tree, path).AsNumpy()
    cols = rdf.keys()
    ROOT.ROOT.DisableImplicitMT()
    return pd.DataFrame(rdf, columns=cols)

def WriteTFile(pdf, path, tree="tree"):
    ROOT.ROOT.EnableImplicitMT()   
    data = {key: pdf[key].values for key in pdf.columns}
    rdf = ROOT.RDF.MakeNumpyDataFrame(data)
    rdf.Snapshot(tree, path)
    ROOT.ROOT.DisableImplicitMT()
    return

def ReadTreeKeys(path, tree="tree"):
    ROOT.ROOT.EnableImplicitMT()
    rdf = rd(tree, path)
    cols = rdf.GetColumnNames()
    ROOT.ROOT.DisableImplicitMT()
    return cols

def GetTFileConfig(path):
    config = {}
    
    ROOT.ROOT.EnableImplicitMT()
    try:
        file = TFile(path)
    except:
        print("Could not open file: " + path)
        ROOT.ROOT.DisableImplicitMT()
        return {}
    
    config["path"] = path
    
    for key in file.GetListOfKeys():
        if key.GetClassName() == 'TTree':
            config['tree'] = key.GetName()
            break
    if 'tree' not in config:
        print("Could not find tree in file: " + path)
        ROOT.ROOT.DisableImplicitMT()
        return {}
    
    tree = file.Get(config['tree'])
    config['entries'] = tree.GetEntries()
    config['features'] = []
    config['types'] = []
    
    for branch in tree.GetListOfBranches():
        config['features'].append(branch.GetName())
        config['types'].append(branch.GetListOfLeaves()[0].GetTypeName())
    
    file.Close()
    
    ROOT.ROOT.DisableImplicitMT()
    return config
    

    
    



