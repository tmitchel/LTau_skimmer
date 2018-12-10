import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description='run on a directory containing directories containing the skimmed ntuples')
parser.add_argument('-p', '--prefix', action = 'store', help = 'name to prefix all directories')
parser.add_argument('-j', '--job', action = 'store', help = 'job type')
parser.add_argument('-l', '--lepton', action = 'store', help = 'which lepton (mt or et)')
args = parser.parse_args()

bkg_samples = {
  'DYJets_lowM'  : ['DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1', 'Z'],
  'DYJets'       : ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1', 'Z'],
  'DYJets_ext1'  : ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v1'             , 'Z'],
  'DYJets1'      : ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1'                 , 'Z'],
  'DYJets1_v2'   : ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v2' , 'Z'],
  'DYJets1_ext1' : ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v1', 'Z'],
  'DYJets1_real' : ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v3_94X_mc2017_realistic', 'Z'],
  'DYJets2'      : ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1', 'Z'],
  'DYJets2_ext'  : ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v1', 'Z'],
  'DYJets3'      : ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1', 'Z'],
  'DYJets3_ext'  : ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v1', 'Z'],
  'DYJets4'      : ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1', 'Z'],
  'DYJets4_real' : ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2_94X_mc2017_realistic', 'Z'],
  'WJets'        : ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v2', 'W'],
  'WJets_ext1'   : ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14_ext1-v2', 'W'],
  'WJets1'       : ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v3', 'W'],
  'WJets2'       : ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v4', 'W'],
  'WJets3'       : ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1', 'W'],
  'WJets4'       : ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1', 'W'],
  'TTLep'        : ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'TTHad'        : ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'TTHad_v2'     : ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_v14-v2', '0'],
  'TTSemi'       : ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'TTSemi_v2'    : ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_v14-v2', '0'],
  'Tbar-tchan'   : ['ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v14-v1', '0'],
  'Tbar-tchan_v2': ['ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v14-v2', '0'],
  'T-tchan'      : ['ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v14-v1', '0'],
  'Tbar-tW'      : ['ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v14-v2', '0'],
  'T-tW'         : ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'T-tW_v2'      : ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v14-v2', '0'],
  'VV'           : ['VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v14-v1', '0'],
  'WW1l1nu2q'    : ['WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'WW1l1nu2q_ext': ['WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14_ext1-v1', '0'],
  'WW1l1nu2q_amc': ['WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v1', '0'],
  'WW2l2nu'      : ['WWTo2L2Nu_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'WW4q'         : ['WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8_v14-v1', '0'],
  'EWKWMinus'    : ['EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v14-v1', '0'],
  'EWKWPlus'     : ['EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v14-v1', '0'],
  'EWKZ2l'       : ['EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8_v14-v1', '0'],
  'EWKZ2nu'      : ['EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8_v14-v1', '0'],
  'WZ3l1nu'      : ['WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_v14-v1', '0'],
  'WZ1l1nu2q'    : ['WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v1', '0'],
  'WZ1l1nu2q_v2' : ['WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v2', '0'],
  'WZ1l3nu'      : ['WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2_v2', '0'],
  'WZ2l2q'       : ['WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v1', '0'],
  'ZZ2l2q'       : ['ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v14-v1', '0'],
  'ZZ2l2nu'      : ['ZZTo2L2Nu_13TeV_powheg_pythia8_v14-v1', '0'],
  'ZZ4l'         : ['ZZTo4L_13TeV_powheg_pythia8_v14-v1', '0'],
  'ZZ4l_v2'      : ['ZZTo4L_13TeV_powheg_pythia8_v14-v2', '0'],
  'ZZ4l_ext1'    : ['ZZTo4L_13TeV_powheg_pythia8_v14_ext1-v1', '0'],
}

sig_samples = {
  'ggH125_v1'  : ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'     , 'Z'],
  'ggH125_ext' : ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v14_ext1-v1', 'Z'],
  'VBF125'     : ['VBFHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'       , 'Z'],
  'WPlus125'   : ['WplusHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'     , 'Z'],
  'WMinus125'  : ['WminusHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'    , 'Z'],
  'ZH125'      : ['ZHToTauTau_M125_13TeV_powheg_pythia8_v14-v1'         , 'Z'],

}

el_data_samples = {
    'datasE-B'     : ['data_SingleElectron_Run2017B-31Mar2018', '0'],
    'datasE-C'     : ['data_SingleElectron_Run2017C-31Mar2018', '0'],
    'datasE-D'     : ['data_SingleElectron_Run2017D-31Mar2018', '0'],
    'datasE-E'     : ['data_SingleElectron_Run2017E-31Mar2018', '0'],
    'datasE-F'     : ['data_SingleElectron_Run2017F-31Mar2018', '0'],
}

mu_data_samples = {
    'datasMu-B'     : ['data_SingleMuon_Run2017B-31Mar2018', '0'],
    'datasMu-C'     : ['data_SingleMuon_Run2017C-31Mar2018', '0'],
    'datasMu-D'     : ['data_SingleMuon_Run2017D-31Mar2018', '0'],
    'datasMu-E'     : ['data_SingleMuon_Run2017E-31Mar2018', '0'],
    'datasMu-F'     : ['data_SingleMuon_Run2017F-31Mar2018', '0'],
}

el_embedded_samples = {
  'embedEl-B' : ['embedded_EmbeddingRun2017B_ElTauFinalState', '0'],
  'embedEl-C' : ['embedded_EmbeddingRun2017C_ElTauFinalState', '0'],
  'embedEl-D' : ['embedded_EmbeddingRun2017D_ElTauFinalState', '0'],
  'embedEl-E' : ['embedded_EmbeddingRun2017E_ElTauFinalState', '0'],
  'embedEl-F' : ['embedded_EmbeddingRun2017F_ElTauFinalState', '0'],
}

mu_embedded_samples = {
  'embedMu-B' : ['embedded_EmbeddingRun2017B_MuTauFinalState', '0'],
  'embedMu-C' : ['embedded_EmbeddingRun2017C_MuTauFinalState', '0'],
  'embedMu-D' : ['embedded_EmbeddingRun2017D_MuTauFinalState', '0'],
  'embedMu-E' : ['embedded_EmbeddingRun2017E_MuTauFinalState', '0'],
  'embedMu-F' : ['embedded_EmbeddingRun2017F_MuTauFinalState', '0'],
}

prefix = args.prefix
jobType = args.job

mc_pref = '/hdfs/store/user/caillol/SMHTT_2017_7nov/'
sig_pref = mc_pref
data_pref = '/hdfs/store/user/caillol/SMHTT2017_data_8nov/'
embed_pref = '/hdfs/store/user/caillol/SMHTT2017_embedded_8nov/'

samples = bkg_samples
pref = mc_pref
if args.job == 'sig':
  samples = sig_samples
  pref = sig_pref
elif args.job == 'data':
  pref = data_pref
  if args.lepton == 'mt':
    samples = mu_data_samples
  elif args.lepton == 'et':
    samples = el_data_samples
elif args.job == 'embed':
  pref = embed_pref
  if args.lepton == 'mt':
    samples = mu_embedded_samples
  elif args.lepton == 'et':
    samples = el_embedded_samples

if args.lepton == 'mt':
  lep = 'et'
elif args.lepton == 'et':
  lep = 'et'

for sample in sorted(samples.keys()):
  recoil = samples[sample][1]
  path = samples[sample][0]
  subprocess.call('python Skimminate.py -sn %s -sd %s --jobName %s -j %s -r %s -l %s -y %s' % (sample, pref+path, prefix, jobType, recoil, lep, '2017'), shell=True)
