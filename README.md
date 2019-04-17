# Higgs to Tau Tau Skimming Code

This repository contains all of the code needed to skim either lepton+tau channel in 2016 or 2017. There is a single binary used to control the skimming for all year/lepton permutations and various headers containing the tree info for each permuation.

##### Table of Contents
[Input File Locations](#files) <br/>
[Output File Locations](#ofiles) <br/>
[Quick Start](#quickstart) <br/>
[Using Condor](#condor) <br/>

<a name="files"/>

## Input File Locations

Here are the locations of all currently used FSA ntuples
- 2016 ntuples
    - Background Monte Carlo: /hdfs/store/user/ndev/LFV_feb18_mc/
    - Signal Monte Carlo: /hdfs/store/user/truggles/SMHTT_signals_may30/
    - Data: /hdfs/store/user/ndev/LFV_reminiaod_feb18/
    - Embedded Electron: /hdfs/store/user/abdollah/MiniAOD_Embed_et_v5/
    - Embedded Muon: /hdfs/store/user/abdollah/MiniAOD_Embed_mt_v5/
    - ggH Only: /hdfs/store/user/truggles/SM-HTT_HTXS_ggH_aug31_v1/

- 2017 ntuples
    - Monte Carlo: /hdfs/store/user/caillol/SMHTT_2017_7nov/
    - Data: /hdfs/store/user/caillol/SMHTT2017_data_8nov/
    - Embedded: /hdfs/store/user/caillol/SMHTT2017_embedded_8nov/

- 2018 ntuples
    - Monte Carlo: /hdfs/store/user/caillol/TauID2018_19dec/
    - Data: /hdfs/store/user/caillol/TauID2018_data_19dec/
    - Embedded: Not Produced Yet

<a name="ofiles"/>

## Output File Locations

Here are the locations of the skims currently being used for studies
- 2016 skims:
    - /hdfs/store/user/tmitchel/etau2016_official_v3-skim
    - /hdfs/store/user/tmitchel/mutau2016_official_v3-skim

- 2017 ntuples
    - /hdfs/store/user/tmitchel/etau2017_official_v3-skim
    - /hdfs/store/user/tmitchel/mutau2017_official_v3-skim

- 2018 ntuples
    - /hdfs/store/user/tmitchel/etau2018_first-look_v2
    - /hdfs/store/user/tmitchel/mutau2018_first-look_v2

<a name="quickstart"/>

## Quick Start

This section is designed so that you can start producing skims by simply copy/pasting the following commands. Refer to other sections of the README for more detailed instructions, if needed.

1. Setup a new CMSSW release
    ```
    cmsrel CMSSW_9_4_0 && cd CMSSW_9_4_0/src && cmsenv
    ```
2. Clone all necessary repositories and get them setup
    - clone this repo
        ```
        git clone -b development ssh://git@gitlab.cern.ch:7999/KState-HEP-HTT/ltau_skimmer.git
        ```
    - get the files needed for recoil corrections
        ```
        git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections
        ```
    - get this super awesome project for reading JSON files in C++ code (thanks nlohmann!)
        ```
        cd ltau_skimmer/ROOT/src
        wget https://github.com/nlohmann/json/releases/download/v3.6.1/json.hpp
        cd -
        ```
    - now compile CMSSW things that need compiling
        ```
        cd $CMSSW_BASE/src
        scram b -j 8
        ```
3. Copy the provided fileMap.json to a place it can be used. This map is used to put the miniAOD filename into the output ROOT file to be used during Pileup reweighting for2017
    ```
    cd ltau_skimmer/ROOT/scripts
    cp fileMap.json $CMSSW_BASE/bin/$SCRAM_ARCH/fileMap.json
    ```
4. Submit skims to condor for a chosen year/lepton/job type
Submitting 2016 samples to be skimmed is done using a single python script to submit multiple jobs.
    ```
    voms-proxy-init --voms=cms --valid=48:00 # get certificate
    python submit2016.py -p mutau2016_stable_v3 -l mt -j sig
    ```
    This will farmout a job to skim all signal samples in the 2016 mutau channel. The output files will be saved to /hdfs/store/user/your_name/mutau2016_stable_v3 or whatever name you provide to the `-p` option. As you can see from the example, the `-l` flag is used to give the lepton type and the `-j` flag is used to pass the job type [sig, bkg1, bkg2, data, embed].

    The script `submit2017.py` can be used in the exact same way to submit jobs for 2017.

These scripts will be combined into a universal script later when I have some free time.

<a name="condor"/>

## Skimming with Wisconsin Condor
The primary way to produce skimmed ntuples using this skimmer is through the use of the Wisconsin Condor system. The scripts in the `ROOT/scripts` directory are used for submitting jobs for processing with Condor.

`Skimminate.py` is used to skim single directory filled with FSA ntuples. Command-line arguments can be used to customize aspects of the job including the input directory and whether to apply recoil corrections. An example command is shown below

```
python Skimminate.py --job bkg --recoil Z --jobName chooseAName --samplename DYJets1 --sampledir /hdfs/store/user/tmitchel/input/DYJets1 -l et -y 2016
```

This command will skim all files in the directory `/hdfs/store/user/tmitchel/input/DYJets1` assuming they are from the dataset `DYJets1`. The corrections appropriate for a background MC sample will be applied and the output files will be stored in `/hdfs/store/user/tmitchel/chooseAName`. The `-l` option tells the script for which channel we want to submit skims and the `-y` option provides the year.

The typical way to produce skims is using the `submit*.py` scripts which allows entire job types to be submitted at once. An example command to skim all embedded samples, with appropriate corrections, is shown below.

```
python submit2016.py -p myName -l et -j embed
```

The output files will be stored in `/hdfs/store/user/tmitchel/myName`.

