paths = ["/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/200GeV/R02/200GeV_Match_R02.root",
       "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/200GeV/R04/200GeV_Match_R04.root",
       "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/200GeV/R06/200GeV_Match_R06.root",
       "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/2760GeV/R02/2760GeV_Match_R02.root",
       "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/2760GeV/R04/2760GeV_Match_R04.root",
       "/lustre/isaac/scratch/tmengel/jet-background-subtraction/datasets/2760GeV/R06/2760GeV_Match_R06.root"]

names = ['AuAu_R02', 'AuAu_R04', 'AuAu_R06', 'PbPb_R02', 'PbPb_R04', 'PbPb_R06']

import json
for path, name in zip(paths, names):
       Config = {}
       DataSetConfig = {}
       DataSetConfig['Name'] = name
       DataSetConfig['Path'] = path
       DataSetConfig['Overwrite'] = True
       DataSetConfig['Split'] = 0.5
              
       ModelConfig = {}
       ModelConfig['Name'] = f'{name}_SNN'
       ModelConfig['Type'] = 'snn'
       ModelConfig['Features'] = ["jet_pt_raw","jet_nparts","jet_angularity", "jet_area","jet_track_pt_0"]
       ModelConfig['Target'] = 'jet_pt_truth'
       ModelConfig['Retrain'] = True
       ModelConfig['Epochs'] = 50
       ModelConfig['Batch Size'] = 512
       ModelConfig['Optimizer'] = 'adam'
       ModelConfig['Loss'] = 'mse'
       
       Model2Config = {}
       Model2Config['Name'] = f'{name}_DNN'
       Model2Config['Type'] = 'dnn'
       Model2Config['Features'] = ["jet_pt_raw","jet_nparts","jet_angularity", "jet_area" ,
            "jet_track_pt_0", "jet_track_pt_1","jet_track_pt_2","jet_track_pt_3","jet_track_pt_4",
            "jet_track_pt_5","jet_track_pt_6","jet_track_pt_7"]
       Model2Config['Target'] = 'jet_pt_truth'
       Model2Config['Retrain'] = True
       Model2Config['Epochs'] = 50
       Model2Config['Batch Size'] = 512
       Model2Config['Optimizer'] = 'adam'
       Model2Config['Loss'] = 'mse'
       
       Config['Datasets'] = DataSetConfig
       Config['Models'] = {}
       Config['Models'].update({ModelConfig['Name']:ModelConfig})
       Config['Models'].update({Model2Config['Name']:Model2Config})
       
       with open(f'Configs/{name}.json', 'w') as f:
              json.dump(Config, f, indent=4)