import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description="run on a directory containing directories containing the skimmed ntuples")
parser.add_argument('-p', '--prefix', action = 'store', help = 'name to prefix all directories')
parser.add_argument('-j', '--job', action = 'store', help = 'job type')
args = parser.parse_args()

bkg_samples_v2 = {
  "DYJets1_ext"  : ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v3_94X_mc2017_realistic', 'Z'],
}

bkg_samples = {

  "DYJets"       : ["DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1", 'Z'],
  'DYJets_ext1'  : ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v1'             , 'Z'],
  "DYJets1"      : ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1'                 , 'Z'],
  "DYJets2"      : ["DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1", 'Z'],
  "DYJets2_ext"  : ["DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v1", 'Z'],
  "DYJets3"      : ["DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1", 'Z'],
  "DYJets3_ext"  : ["DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v1", 'Z'],
  "DYJets4"      : ["DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2_94X_mc2017_realistic", 'Z'],
  "WJets"        : ["WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v2", 'W'],
  "WJets_ext1"   : ["WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v2", 'W'],
  "WJets1"       : ["W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v3", 'W'],
  "WJets2"       : ["W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v4", 'W'],
  "WJets3"       : ["W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1", 'W'],
  "WJets4"       : ["W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1", 'W'],
  "TTLep"        : ["TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_v14-v1", '0'],
  "TTHad"        : ["TTToHadronic_TuneCP5_13TeV-powheg-pythia8_v14-v2", '0'],
  "TTSemi"       : ["TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_v14-v1", '0'],
  "Tbar-tchan"   : ["ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v14-v2", '0'],
  "T-tchan"      : ["ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v14-v1", '0'],
  "Tbar-tW"      : ["ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v14-v2", '0'],
  "T-tW"         : ["ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v14-v2", '0'],
  "WW1l1nu2q"    : ["WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14-v1", '0'],
  "WW1l1nu2q_ext": ["WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14_ext1-v1", '0'],
  'WW2l2nu'      : ['WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'WW4q'         : ['WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  "ZZ2l2q"       : ["ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v1", '0'],
  "ZZ2l2nu"      : ["ZZTo2L2Nu_13TeV_powheg_pythia8_v14-v1", '0'],
  "EWKWMinus"    : ["EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v14-v1", '0'],
  "EWKWPlus"     : ["EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v14-v1", '0'],
  "EWKZ2l"       : ["EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8_v14-v1", '0'],
  "EWKZ2nu"      : ["EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8_v14-v1", '0'],
  "WZ3l1nu"      : ["WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_v14-v1", '0'],
  "WZ1l1nu2q"    : ["WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v2", '0'],
  "WZ1l3nu"      : ["WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2_v2", '0'],
  "WZ2l2q"       : ["WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v1", '0'],
  'WW'           : ['WW_TuneCP5_13TeV-pythia8_v14-v1', '0'],
  'WZ'           : ['WZ_TuneCP5_13TeV-pythia8_v14-v1', '0'],
  'ZZ'           : ['ZZ_TuneCP5_13TeV-pythia8_v14-v1', '0'],
}

sig_samples = {
  "ggH125_v1"  : ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'     , 'Z'],
  "ggH125_ext" : ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v14_ext1-v1', 'Z'],
  "VBF125"     : ['VBFHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'       , 'Z'],
  "WPlus125"   : ['WplusHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'     , 'Z'],
  "WMinus125"  : ['WminusHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'    , 'Z'],
  "ZH125"      : ['ZHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'         , 'Z'],
  "ttH125"     : ['ttHToTauTau_M125_TuneCP5_13TeV-powheg-pythia8_v14-v1' , 'Z'],

}

data_samples = {
    "datasE-B"     : ["data_SingleElectron_Run2017B-17Nov2017", '0'],
    "datasE-C"     : ["data_SingleElectron_Run2017C-17Nov2017", '0'],
    "datasE-D"     : ["data_SingleElectron_Run2017D-17Nov2017", '0'],
    "datasE-E"     : ["data_SingleElectron_Run2017E-17Nov2017", '0'],
    "datasE-F"     : ["data_SingleElectron_Run2017F-17Nov2017", '0'],
}

embedded_samples = {
  'embed-B' : ['EmbeddingRun2016B-v2', '0'],
  'embed-C' : ['EmbeddingRun2016C-v2', '0'],
  'embed-D' : ['EmbeddingRun2016D-v2', '0'],
  'embed-E' : ['EmbeddingRun2016E-v2', '0'],
  'embed-F' : ['EmbeddingRun2016F-v2', '0'],
  'embed-G' : ['EmbeddingRun2016G-v2', '0'],
  'embed-H' : ['EmbeddingRun2016H-v2', '0'],
}

prefix = args.prefix
jobType = args.job

mc_pref = "/hdfs/store/user/ymaravin/MC_2017/"
sig_pref = mc_pref
data_pref = "/hdfs/store/user/ymaravin/Data_2017/"
other_mc_pref = "/hdfs/store/user/tmitchel/FSA2017_allMC/MC_2017/"

samples = bkg_samples
pref = mc_pref
if args.job == 'sig':
  samples = sig_samples
  pref = sig_pref
elif args.job == 'data':
  samples = data_samples
  pref = data_pref
elif args.job == 'embed':
  samples = embedded_samples
  pref = embed_pref
elif args.job == 'bkg2':
  samples = bkg_samples_v2
  pref = other_mc_pref

for sample in sorted(samples.keys()):
  recoil = samples[sample][1]
  path = samples[sample][0]
  if 'ggHtoTauTau125' in sample:
    pref = ggH_pref
  subprocess.call('python Skimminate.py -sn %s -sd %s --jobName %s -j %s -r %s' % (sample, pref+path, prefix, jobType, recoil), shell=True)
