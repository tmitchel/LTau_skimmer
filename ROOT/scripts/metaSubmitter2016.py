import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description="run on a directory containing directories containing the skimmed ntuples")
parser.add_argument('-p', '--prefix', action = 'store', help = 'name to prefix all directories')
parser.add_argument('-l', '--lepton', action = 'store', help = 'which lepton (mt or et)')
parser.add_argument('-j', '--job', action = 'store', help = 'job type')
args = parser.parse_args()

bkg_samples = {
  "DYJets"       : ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v2", 'Z'],
  "DYJets_ext1"  : ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext2-v1", 'Z'],
  "DYJets1"      : ["DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1", 'Z'],
  "DYJets2"      : ["DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1", 'Z'],
  "DYJets3"      : ["DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1", 'Z'],
  "DYJets4"      : ["DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1", 'Z'],
  "WJets"        : ["WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1",  'W'],
  "WJets_ext1"   : ["WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext2-v1", 'W'],
  "WJets1"       : ["W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1", 'W'],
  "WJets2"       : ["W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1",  'W'],
  "WJets2_ext1"  : ["W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v1", 'W'],
  "WJets3"       : ["W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1",  'W'],
  "WJets3_ext1"  : ["W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v1", 'W'],
  "WJets4"       : ["W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1",  'W'],
  "WJets4_ext1"  : ["W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v1",  'W'],
  "WJets4_ext2"  : ["W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext2-v1", 'W'],
  "TT"           : ["TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v6-v1", '0'],
  "VV2l2nu"      : ["VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v6_ext1-v1", '0'],
  "Tbar-tchan"   : ["ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v6-v1", '0'],
  "T-tchan"      : ["ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v6-v1", '0'],
  "Tbar-tW"      : ["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v6_ext1-v1", '0'],
  "T-tW"         : ["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v6_ext1-v1", '0'],
  "WW1l1nu2q"    : ["WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1", '0'],
  "ZZ2l2q"       : ["ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1", '0'],
  "ZZ4l"         : ["ZZTo4L_13TeV-amcatnloFXFX-pythia8_v6_ext1-v1", '0'],
  "EWKWMinus"    : ["EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6-v1",  '0'],
  "EWKMinus_ext1": ["EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext1-v1",  '0'],
  "EWKMinus_ext2": ["EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext2-v1", '0'],
  "EWKWPlus"     : ["EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6-v1",  '0'],
  "EWKPlus_ext1" : ["EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext1-v1",  '0'],
  "EWKPlus_ext2" : ["EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v6_ext2-v1", '0'],
  "EWKZ2l"       : ["EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v6-v1",  '0'],
  "EWKZ2l_ext1"  : ["EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v6_ext1-v1",  '0'],
  "EWKZ2l_ext2"  : ["EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v6_ext2-v1", '0'],
  "EWKZ2nu"      : ["EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6-v1",  '0'],
  "EWKZ2nu_ext1" : ["EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6_ext1-v1",  '0'],
  "EWKZ2nu_ext2" : ["EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v6_ext2-v1", '0'],
  #"ggH_WW125"    : ["GluGluHToWWTo2L2Nu_M125_13TeV_powheg_pythia8_v6-v1", '0'],
  #"VBF_WW125"    : ["VBFHToWWTo2L2Nu_M125_13TeV_powheg_pythia8_v6-v1", '0'],
  #"WGLNu"        : ["WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext1-v1",  '0'],
  #"WGLNu_ext1"   : ["WGToLNuG_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v6_ext2-v1", '0'],
  #"WGstarEE"     : ["WGstarToLNuEE_012Jets_13TeV-madgraph_v6-v1", '0'],
  #"WGstarMuMu"   : ["WGstarToLNuMuMu_012Jets_13TeV-madgraph_v6-v1", '0'],
  #"WWW"          : ["WWW_4F_TuneCUETP8M1_13TeV-amcatnlo-pythia8_v6-v1", '0'],
  "WZ3l1nu"      : ["WZJToLLLNu_TuneCUETP8M1_13TeV-amcnlo-pythia8_v6-v1", '0'],
  "WZ1l1nu2q"    : ["WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v3", '0'],
  "WZ1l3nu"      : ["WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1", '0'],
  "WZ2l2q"       : ["WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v6-v1", '0'],
}

sig_samples = {
  "ggHtoTauTau125_v1"  : ["GluGluHToTauTau_M125_13TeV_powheg_pythia8_v6-v1", 'Z'],
  "ggHtoTauTau125_v2"  : ["GluGluHToTauTau_M125_13TeV_powheg_pythia8_v6-v2", 'Z'],
  "VBFHtoTauTau125" : ["VBFHToTauTau_M125_13TeV_powheg_pythia8_v6-v1", 'Z'],
  "WPlusHTauTau125" : ["WplusHToTauTau_M125_13TeV_powheg_pythia8_v6-v1", '0'],
  "WMinusHTauTau125": ["WminusHToTauTau_M125_13TeV_powheg_pythia8_v6-v1", '0'],
  "ZHTauTau125"     : ["ZHToTauTau_M125_13TeV_powheg_pythia8_v6-v1", '0'],
}

el_data_samples = {
    "datasE-B"     : ["data_SingleElectron_Run2016B_v1", '0'],
    "datasE-B_ext1": ["data_SingleElectron_Run2016B_v2", '0'],
    "datasE-C"     : ["data_SingleElectron_Run2016C", '0'],
    "datasE-D"     : ["data_SingleElectron_Run2016D", '0'],
    "datasE-E"     : ["data_SingleElectron_Run2016E", '0'],
    "datasE-F"     : ["data_SingleElectron_Run2016F", '0'],
    "datasE-G"     : ["data_SingleElectron_Run2016G", '0'],
    "datasE-H"     : ["data_SingleElectron_Run2016H_v2", '0'],
    "datasE-H_ext1": ["data_SingleElectron_Run2016H_v3", '0'],
}

mu_data_samples = {
    "datasMu-B"     : ["data_SingleMuon_Run2015B_v1", '0'],
    "datasMu-B_ext1": ["data_SingleMuon_Run2016B_v2", '0'],
    "datasMu-C"     : ["data_SingleMuon_Run2016C", '0'],
    "datasMu-D"     : ["data_SingleMuon_Run2016D", '0'],
    "datasMu-E"     : ["data_SingleMuon_Run2016E", '0'],
    "datasMu-F"     : ["data_SingleMuon_Run2016F", '0'],
    "datasMu-G"     : ["data_SingleMuon_Run2016G", '0'],
    "datasMu-H"     : ["data_SingleMuon_Run2016H_v2", '0'],
    "datasMu-H_ext1": ["data_SingleMuon_Run2016H_v3", '0'],
}

el_embedded_samples = {
  'embedEl-B' : ['EmbeddingRun2016B-v2', '0'],
  'embedEl-C' : ['EmbeddingRun2016C-v2', '0'],
  'embedEl-D' : ['EmbeddingRun2016D-v2', '0'],
  'embedEl-E' : ['EmbeddingRun2016E-v2', '0'],
  'embedEl-F' : ['EmbeddingRun2016F-v2', '0'],
  'embedEl-G' : ['EmbeddingRun2016G-v2', '0'],
  'embedEl-H' : ['EmbeddingRun2016H-v2', '0'],
}

mu_embedded_samples = {
  'embedMu-B' : ['EmbeddingRun2016B-v2', '0'],
  'embedMu-C' : ['EmbeddingRun2016C-v2', '0'],
  'embedMu-D' : ['EmbeddingRun2016D-v2', '0'],
  'embedMu-E' : ['EmbeddingRun2016E-v2', '0'],
  'embedMu-F' : ['EmbeddingRun2016F-v2', '0'],
  'embedMu-G' : ['EmbeddingRun2016G-v2', '0'],
  'embedMu-H' : ['EmbeddingRun2016H-v2', '0'],
}

prefix = args.prefix
jobType = args.job

bkg_pref = "/hdfs/store/user/ndev/LFV_feb18_mc/"
sig_pref = "/hdfs/store/user/truggles/SMHTT_signals_may30/"
data_pref = "/hdfs/store/user/ndev/LFV_reminiaod_feb18/"
el_embed_pref = '/hdfs/store/user/abdollah/MiniAOD_Embed_et_v5/'
mu_embed_pref = '/hdfs/store/user/abdollah/MiniAOD_Embed_mt_v5/'
ggH_pref = '/hdfs/store/user/truggles/SM-HTT_HTXS_ggH_aug31_v1/'

samples = bkg_samples
pref = bkg_pref
if args.job == 'sig':
  pref = sig_pref
  samples = sig_samples
elif args.job == 'data':
  pref = data_pref
  if args.lepton == 'mt':
    samples = mu_data_samples
  elif args.lepton == 'et':
    samples = el_data_samples
elif args.job == 'embed':
  if args.lepton == 'mt':
    samples = mu_embedded_samples
    pref = mu_embed_pref
  elif args.lepton == 'et':
    samples = el_embedded_samples
    pref = el_embed_pref

if args.lepton == 'mt':
  lep = 'MT'
elif args.lepton == 'et':
  lep = 'ET'

for sample in sorted(samples.keys()):
  recoil = samples[sample][1]
  path = samples[sample][0]
  if 'ggHtoTauTau125' in sample:
    pref = ggH_pref
  subprocess.call('python Skimminate.py -sn %s -sd %s --jobName %s -j %s -r %s -l %s2016' % (sample, pref+path, prefix, jobType, recoil, lep), shell=True)
