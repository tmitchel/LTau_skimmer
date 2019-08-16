import os
import condor_handler as ch
from argparse import ArgumentParser

parser = ArgumentParser(
    description='run on a directory containing directories containing the skimmed ntuples')
parser.add_argument('-p', '--prefix', action='store',
                    help='name to prefix all directories')
parser.add_argument('-j', '--job', action='store', help='job type')
parser.add_argument('-l', '--lepton', action='store',
                    help='which lepton (mt or et)')
args = parser.parse_args()

test_batch = {
    'DYJets1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
}

temp = {
    'WJets': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
}

bkg_samples_batch1 = {
    'DYJets1_v1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets1_v2': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'DYJets2_v1': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets2_v2': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'DYJets3_v1': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets3_v2': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'DYJets4': ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets_lowM': ['DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'DYJets': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],

    'WJets1': ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets2_v1': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'WJets2_v2': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets3_v1': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'WJets3_v2': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets4_v1': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'WJets4_v2': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
}
bkg_samples_batch2 = {
    'EWKWMinus': ['EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'EWKWPlus': ['EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'EWKZ2l': ['EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'EWKZ2nu': ['EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'Tbar-tchan': ['ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'T-tchan': ['ST_t-channel_top_5f_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'Tbar-tW': ['ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext1-v1', '0'],
    'T-tW': ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext1-v1', '0'],
    'TTLep': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'TTHad': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'TTSemi': ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WW_v1': ['WW_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WW_v2': ['WW_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'WZ_v1': ['WZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WZ_v2': ['WZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v3', '0'],
    'ZZ': ['ZZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
}

sig_samples = {
    'ggh125_powheg': ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_powheg': ['VBFHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15_ext1-v1', 'Z'],
    'wminus125_powheg': ['WminusHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wplus125_powheg': ['WplusHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_powheg': ['ZHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
}

# Needs to be updated
ac_samples = {
    'ggH_a1': ['JJHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'ggH_a3': ['JJHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'ggH_a3int': ['JJHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_a1': ['VBFHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_a2': ['VBFHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_a2int': ['VBFHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_a3': ['VBFHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_a3int': ['VBFHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_l1': ['VBFHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_l1int': ['VBFHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_l1zg': ['VBFHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf_l1zgint': ['VBFHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', 'Z'],
    'wh_a1': ['WHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'wh_a2': ['WHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'wh_a2int': ['WHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v3', '0'],
    'wh_a3': ['WHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v3', '0'],
    'wh_a3int': ['WHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'wh_l1': ['WHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'wh_l1int': ['WHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_a1': ['ZHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_a2': ['ZHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_a2int': ['ZHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_a3': ['ZHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_a3int': ['ZHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_l1': ['ZHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_l1int': ['ZHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
    'zh_l1zgint': ['ZHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8__94X_mc2017_realistic_v14-v2', '0'],
}

# Needs to be updated
madgraph_samples = {
    "ggH_TwoJet_madgraph": ["GluGluToHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8__94X_mc2017_realistic_v14-v1", 'Z'],
    "ggH_PS_TwoJet_madgraph": ["GluGluToPseudoscalarHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8__94X_mc2017_realistic_v14-v1", 'Z'],
    "ggH_Maxmix_TwoJet_madgraph": ["GluGluToMaxmixHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8__94X_mc2017_realistic_v14-v1", 'Z'],
}

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

# Need to be updated
el_embedded_samples = {
    'embedEl-A': ['embedded_EmbeddingRun2017B_ElTauFinalState', '0'],
    'embedEl-B': ['embedded_EmbeddingRun2017B_ElTauFinalState', '0'],
    'embedEl-C': ['embedded_EmbeddingRun2017C_ElTauFinalState', '0'],
    'embedEl-D': ['embedded_EmbeddingRun2017D_ElTauFinalState', '0'],
}

# Need to be updated
mu_embedded_samples = {
    'embedMu-A': ['embedded_EmbeddingRun2017B_MuTauFinalState', '0'],
    'embedMu-B': ['embedded_EmbeddingRun2017B_MuTauFinalState', '0'],
    'embedMu-C': ['embedded_EmbeddingRun2017C_MuTauFinalState', '0'],
    'embedMu-D': ['embedded_EmbeddingRun2017D_MuTauFinalState', '0'],
}

prefix = args.prefix
jobType = args.job

bkg_pref = '/hdfs/store/user/caillol/SMHTT_legacy_2018_240419/'
sig_pref = bkg_pref
data_pref = '/hdfs/store/user/caillol/SMHTT_legacy_2018_data_24042019/'
# Need updated
ac_pref = '/hdfs/store/user/ymaravin/ac2017_v2/'
madgraph_pref = '/hdfs/store/user/ymaravin/mg2017_v1/'
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
    'temp': [bkg_pref, temp]
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
    pref = settings[args.job][0]

    command = '$CMSSW_BASE/bin/$SCRAM_ARCH/uniSkim -d %s -j %s -r %s -y %s -l %s -i input_file.root -o \'$OUTPUT\'' % (
        pref+path, jobType, recoil, '2017', lep)
    ch.submit_command(command, prefix, pref+path, sample, use_input='-n ', dryrun=False)
