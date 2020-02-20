# Higgs to Tau Tau Skimming Code

This repository contains all of the code needed to skim either lepton+tau channel in 2016 or 2017. There is a single binary used to control the skimming for all year/lepton permutations and various headers containing the tree info for each permuation.

##### Table of Contents
[Input File Locations](#files) <br/>
[Quick Start](#quickstart) <br/>
[Using Condor](#condor) <br/>

<a name="files"/>

## Input File Locations

Here are the locations of all currently used FSA ntuples
- 2016 ntuples
    - Background Monte Carlo: `/hdfs/store/user/aloeliger/SMHTT_2016_20nov/`
    - Signal Monte Carlo: 
        - `/hdfs/store/user/abdollah/FSA_2016_AC_XMass/`
        - `/hdfs/store/user/aloeliger/SMHTT_2016_20nov_ggH2/`
    - Data: `/hdfs/store/user/aloeliger/SMHTT_2016_data_18Sep/`
    - Embedded: `/hdfs/store/user/aloeliger/SMHTT_2016_embedded_18sep/`

- 2017 ntuples
    - Background Monte Carlo: 
        - `/hdfs/store/user/tmitchel/SMHTT_2017_legacy_mc_v3p2_mt/`
        - `/hdfs/store/user/tmitchel/SMHTT_2017_legacy_mc_v3p2_et/`
    - Signal Monte Carlo: 
        - `/hdfs/store/user/abdollah/FSA_2017_AC_XMass/`
        - `/hdfs/store/user/tmitchel/SMHTT_2017_legacy_mc_v3p2_ggh-only/`
    - Data: `/hdfs/store/user/tmitchel/SMHTT_2017_legacy_data_v3/`
    - Embedded: `/hdfs/store/user/tmitchel/SMHTT_2017_legacy_embedded_v3/`

- 2018 ntuples (old global tag)
    - Background Monte Carlo: 
        - `/hdfs/store/user/caillol/SMHTT_2018_20nov_mc/`
        - `/hdfs/store/user/caillol/SMHTT_2018_20nov_highMem_etmt_mc/`
    - Signal Monte Carlo: `/hdfs/store/user/abdollah/FSA_2018_AC_XMass`/
    - Data: `/hdfs/store/user/caillol/SMHTT_2018_31oct_data/`
    - Embedded: `/hdfs/store/user/caillol/SMHTT_2018_20nov_embedded/`

- 2018 ntuples (new global tag)
    - all in progress

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
        git clone -b development git@github.com:tmitchel/LTau_skimmer.git
        ```
    - get the files needed for recoil corrections
        ```
        git clone https://github.com/CMS-HTT/RecoilCorrections.git HTT-utilities/RecoilCorrections
        ```
    - now compile CMSSW things that need compiling
        ```
        cd ${CMSSW_BASE}/src
        scram b -j 8
        ```
    - finally, download Tau ES files
        ```
        cd ${CMSSW_BASE}/src/LTau_skimmer
        source setup.sh
        ```
3. Submit skims to condor for a chosen year/lepton/job type
Submitting 2016 samples to be skimmed is done using a single python script to submit multiple jobs.
    ```
    voms-proxy-init --voms=cms --valid=48:00 # get certificate
    python submit2016.py -p mutau2016_legacy_v1 -l mt -j sig
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

