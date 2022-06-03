cd ../src/pre-processed-data/jet-trees-formatted
cd ../../../data
echo $PWD

./FormatRawJetTrees 0.2 0 # R = 0.2 0-10%
./FormatRawJetTrees 0.2 1 # R = 0.2 20-40%
./FormatRawJetTrees 0.4 0 # R = 0.4 0-10%
./FormatRawJetTrees 0.4 1 # R = 0.4 20-40%