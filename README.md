# Information about the skimmer

## Setting up the skimmer
Before the skimmer can be used, the ROOT files containing recoil corrections need to be downloaded. To do so, follow these commands

```
cd ${CMSSW_BASE}/src
cmsenv
git clone https://github.com/CMS-HTT/RecoilCorrections.git  HTT-utilities/RecoilCorrections 
```

The skimmer also requires 

Then, to compile

```
cd ${CMSSW_BASE}/src/ETau_Skimmer/ROOT
scram b -j 8
```

## Skimming with Wisconsin Condor
The primary way to produce skimmed ntuples using this skimmer is through the use of the Wisconsin Condor system. The scripts in the `ROOT/scripts` directory are used for submitting jobs for processing with Condor.

`Skimminate.py` is used to skim single directory filled with FSA ntuples. Command-line arguments can be used to customize aspects of the job including the input directory and whether to apply recoil corrections. An example command is shown below

```
python Skimminate.py --job bkg --recoil Z --jobName chooseAName --samplename DYJets1 --sampledir /hdfs/store/user/tmitchel/input/DYJets1 
```

This command will skim all files in the directory `/hdfs/store/user/tmitchel/input/DYJets1` assuming they are from the dataset `DYJets1`. The corrections appropriate for a background MC sample will be applied and the output files will be stored in `/hdfs/store/user/tmitchel/chooseAName`.

The typical way to produce skims is using the `metaSubmitter.py` script which allows entire job types to be submitted at once. An example command to skim all embedded samples, with appropriate corrections, is shown below.

```
python metaSubmitter.py -j embed -p myName
```

The output files will be stored in `/hdfs/store/user/tmitchel/myName`.

