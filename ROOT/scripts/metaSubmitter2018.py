import os
import subprocess
from argparse import ArgumentParser

parser = ArgumentParser(
    description='run on a directory containing directories containing the skimmed ntuples')
parser.add_argument('-p', '--prefix', action='store',
                    help='name to prefix all directories')
parser.add_argument('-j', '--job', action='store', help='job type')
parser.add_argument('-l', '--lepton', action='store',
                    help='which lepton (mt or et)')
args = parser.parse_args()

bkg_samples_batch1 = {
    'DYJets': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'Z'],
    'DYJets1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'Z'],
    'DYJets1_v2': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'Z'],
    'DYJets2': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'Z'],
    'DYJets2_v2': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'Z'],
    'DYJets3': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'Z'],
    'DYJets3_v2': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'Z'],
    'DYJets4': ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'Z'],
    'WJets': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'W'],
    'WJets1': ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'W'],
    'WJets2': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'W'],
    'WJets2_v2': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'W'],
    'WJets3': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'W'],
    'WJets3_v2': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'W'],
    'WJets4': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v1', 'W'],
    'WJets4_v2': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_realistic_v15-v2', 'W'],
}
bkg_samples_batch2 = {
    'TTLep': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_realistic_v15-v1', '0'],
    'TTHad': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_realistic_v15-v1', '0'],
    'TTSemi': ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_realistic_v15-v1', '0'],
    'T-tW': ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_realistic_v15_ext1-v1', '0'],
    'WW': ['WW_TuneCP5_13TeV-pythia8_realistic_v15-v1', '0'],
    'WW_v2': ['WW_TuneCP5_13TeV-pythia8_realistic_v15-v2', '0'],
    'WZ': ['WZ_TuneCP5_13TeV-pythia8_realistic_v15-v1', '0'],
    'WZ_v3': ['WZ_TuneCP5_13TeV-pythia8_realistic_v15-v3', '0'],
    'ZZ': ['ZZ_TuneCP5_13TeV-pythia8_realistic_v15-v2', '0'],
    'EWKWMinus': ['EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_realistic_v15-v1', '0'],
    'EWKWPlus': ['EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_realistic_v15-v1', '0'],
}

sig_samples = {
    'ggH125_v1': ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_realistic_v15-v2', 'Z'],
    'VBF125': ['VBFHToTauTau_M125_13TeV_powheg_pythia8_realistic_v15_ext1-v1', 'Z'],
    'WPlus125': ['WplusHToTauTau_M125_13TeV_powheg_pythia8_realistic_v15-v2', 'Z'],
    'WMinus125': ['WminusHToTauTau_M125_13TeV_powheg_pythia8_realistic_v15-v2', 'Z'],
    'ZH125': ['ZHToTauTau_M125_13TeV_powheg_pythia8_realistic_v15-v2', 'Z'],
}

# ac_samples = {
#     'ggH_a1': ['JJHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'ggH_a3': ['JJHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'ggH_a3int': ['JJHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_a1': ['VBFHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_a2': ['VBFHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_a2int': ['VBFHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_a3': ['VBFHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_a3int': ['VBFHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_l1': ['VBFHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_l1int': ['VBFHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_l1zg': ['VBFHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'vbf_l1zgint': ['VBFHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
#     'wh_a1': ['WHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'wh_a2': ['WHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'wh_a2int': ['WHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v3', '0'],
#     'wh_a3': ['WHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v3', '0'],
#     'wh_a3int': ['WHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'wh_l1': ['WHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'wh_l1int': ['WHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_a1': ['ZHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_a2': ['ZHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_a2int': ['ZHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_a3': ['ZHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_a3int': ['ZHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_l1': ['ZHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_l1int': ['ZHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     'zh_l1zgint': ['ZHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     # 'ttHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     # 'ttHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
#     # 'ttHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
# }

el_data_samples = {
    'datasE-A': ['data_EGamma_Run2018A-17Sep2018', '0'],
    'datasE-B': ['data_EGamma_Run2018B-17Sep2018', '0'],
    'datasE-C': ['data_EGamma_Run2018C-17Sep2018', '0'],
    'datasE-D': ['data_EGamma_Run2018D-PromptReco', '0'],
}

mu_data_samples = {
    'datasMu-A': ['data_SingleMuon_Run2018A-17Sep2018', '0'],
    'datasMu-B': ['data_SingleMuon_Run2018B-17Sep2018', '0'],
    'datasMu-C': ['data_SingleMuon_Run2018C-17Sep2018', '0'],
    'datasMu-D': ['data_SingleMuon_Run2018D-PromptReco', '0'],
}

# el_embedded_samples = {
#     'embedEl-B': ['embedded_EmbeddingRun2017B_ElTauFinalState', '0'],
#     'embedEl-C': ['embedded_EmbeddingRun2017C_ElTauFinalState', '0'],
#     'embedEl-D': ['embedded_EmbeddingRun2017D_ElTauFinalState', '0'],
#     'embedEl-E': ['embedded_EmbeddingRun2017E_ElTauFinalState', '0'],
#     'embedEl-F': ['embedded_EmbeddingRun2017F_ElTauFinalState', '0'],
# }

# mu_embedded_samples = {
#     'embedMu-B': ['embedded_EmbeddingRun2017B_MuTauFinalState', '0'],
#     'embedMu-C': ['embedded_EmbeddingRun2017C_MuTauFinalState', '0'],
#     'embedMu-D': ['embedded_EmbeddingRun2017D_MuTauFinalState', '0'],
#     'embedMu-E': ['embedded_EmbeddingRun2017E_MuTauFinalState', '0'],
#     'embedMu-F': ['embedded_EmbeddingRun2017F_MuTauFinalState', '0'],
# }

prefix = args.prefix
jobType = args.job

bkg_pref = '/hdfs/store/user/caillol/TauID2018_19dec/'
sig_pref = bkg_pref
# ac_pref = '/hdfs/store/user/ymaravin/ac2017_v2/'
data_pref = '/hdfs/store/user/caillol/TauID2018_data_19dec/'
# embed_pref = '/hdfs/store/user/caillol/SMHTT2017_embedded_8nov/'

settings = {
    'sig': [sig_pref, sig_samples],
    # 'ac': [ac_pref, ac_samples],
    'dataMu': [data_pref, mu_data_samples],
    'dataEl': [data_pref, el_data_samples],
    # 'embedMu': [embed_pref, mu_embedded_samples],
    # 'embedEl': [embed_pref, el_embedded_samples],
    'bkg1': [bkg_pref, bkg_samples_batch1],
    'bkg2': [bkg_pref, bkg_samples_batch2],
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
    subprocess.call('python Skimminate.py -sn %s -sd %s --jobName %s -j %s -r %s -l %s -y %s' %
                    (sample, pref+path, prefix, jobType, recoil, lep, '2018'), shell=True)
