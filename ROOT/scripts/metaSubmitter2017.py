import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(description='run on a directory containing directories containing the skimmed ntuples')
parser.add_argument('-p', '--prefix', action = 'store', help = 'name to prefix all directories')
parser.add_argument('-j', '--job', action = 'store', help = 'job type')
parser.add_argument('-l', '--lepton', action = 'store', help = 'which lepton (mt or et)')
args = parser.parse_args()

test_batch = {
  'DYJets1'      : ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v14-v1'                 , 'Z'],
    }

bkg_samples_batch1 = {
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
}
bkg_samples_batch2 = {
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

ac_samples = {
  'ggH_a1'      : ['JJHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'ggH_a3'      : ['JJHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'ggH_a3int'   : ['JJHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_a1'      : ['VBFHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_a2'      : ['VBFHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_a2int'   : ['VBFHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_a3'      : ['VBFHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_a3int'   : ['VBFHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_l1'      : ['VBFHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_l1int'   : ['VBFHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_l1zg'    : ['VBFHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'vbf_l1zgint' : ['VBFHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
  'wh_a1'       : ['WHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'wh_a2'       : ['WHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'wh_a2int'    : ['WHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v3', '0'],
  'wh_a3'       : ['WHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v3', '0'],
  'wh_a3int'    : ['WHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'wh_l1'       : ['WHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'wh_l1int'    : ['WHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_a1'       : ['ZHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_a2'       : ['ZHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_a2int'    : ['ZHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_a3'       : ['ZHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_a3int'    : ['ZHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_l1'       : ['ZHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_l1int'    : ['ZHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  'zh_l1zgint'  : ['ZHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  # 'ttHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  # 'ttHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
  # 'ttHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
} 

madgraph_samples = {
    "ggH_TwoJet_madgraph"        : ["GluGluToHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8__94X_mc2017_realistic_v14-v1", 'Z'],
    "ggH_PS_TwoJet_madgraph"     : ["GluGluToPseudoscalarHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8__94X_mc2017_realistic_v14-v1", 'Z'],
    "ggH_Maxmix_TwoJet_madgraph" : ["GluGluToMaxmixHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8__94X_mc2017_realistic_v14-v1", 'Z'],
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

bkg_pref = '/hdfs/store/user/caillol/SMHTT_2017_7nov/'
sig_pref = bkg_pref
ac_pref = '/hdfs/store/user/ymaravin/ac2017_v2/'
madgraph_pref = '/hdfs/store/user/ymaravin/mg2017_v1/'
data_pref = '/hdfs/store/user/caillol/SMHTT2017_data_8nov/'
embed_pref = '/hdfs/store/user/caillol/SMHTT2017_embedded_8nov/'

settings = {
  'sig': [sig_pref, sig_samples],
  'ac': [ac_pref, ac_samples],
  'madgraph': [madgraph_pref, madgraph_samples],
  'dataMu': [data_pref, mu_data_samples],
  'dataEl': [data_pref, el_data_samples],
  'embedMu': [embed_pref, mu_embedded_samples],
  'embedEl': [embed_pref, el_embedded_samples],
  'bkg1': [bkg_pref, bkg_samples_batch1],
  'bkg2': [bkg_pref, bkg_samples_batch2],
  'test': [bkg_pref, test_batch],
}

pref = settings[args.job][0]
samples = settings[args.job][1]

if args.lepton == 'mt':
  lep = 'mt'
elif args.lepton == 'et':
  lep = 'et'

for sample in sorted(samples.keys()):
  recoil = samples[sample][1]
  path = samples[sample][0]
  subprocess.call('python Skimminate.py -sn %s -sd %s --jobName %s -j %s -r %s -l %s -y %s' % (sample, pref+path, prefix, jobType, recoil, lep, '2017'), shell=True)
