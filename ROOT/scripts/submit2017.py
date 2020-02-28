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
    'DYJets1_ext1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
}

bkg_samples = {
    'DYJets1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', 'Z'],
    'DYJets2': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v2', 'Z'],
    'DYJets3': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v2', 'Z'],
    'DYJets4_v1': ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'DYJets4_v2': ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_v2_94X_mc2017_realistic_v14-v2', 'Z'],
    'DYJets_v1': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'DYJets_v2': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],

    'WJets1_v1': ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', 'W'],
    'WJets1_v2': ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4', 'W'],
    'WJets2_v1': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4', 'W'],
    'WJets2_v2': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v5', 'W'],
    'WJets3': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'W'],
    'WJets4': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'W'],
    'WJets_v1': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'W'],
    'WJets_v2': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', 'W'],
    'WJets_v3': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2', 'W'],
    'EWKWMinus': ['EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'W'],
    'EWKWPlus': ['EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'W'],
    'EWKZ2l': ['EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'Z'],
    'EWKZ2nu': ['EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'Z'],

    'Tbar-tchan_v1': ['ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'Tbar-tchan_v2': ['ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'T-tchan': ['ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', '0'],
    'Tbar-tW': ['ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'T-tW_v1': ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'T-tW_v2': ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'TTLep': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', '0'],
    'TTHad': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],
    'TTSemi': ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', '0'],
    'WW_v1': ['WW_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'WW_v2': ['WW_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'WZ': ['WZ_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'ZZ': ['ZZ_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],

    'VV': ['VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'WW1l1nu2q': ['WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'WZ1l1nu2q_v1': ['WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'WZ1l1nu2q_v2': ['WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v2-PU2017_12Apr2018_old_pmx_94X_mc2017_realistic_v14-v1', '0'],
    'WZ1l3nu': ['WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'WZ2l2q': ['WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'WZ3l1nu': ['WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', '0'],
    'ZZ2l2q': ['ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'ZZ4l': ['ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],

    'vbf125_powheg': ['VBFHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', 'Z'],
    'wplus125_powheg': ['WplusHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'wminus125_powheg': ['WminusHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'zh125_powheg': ['ZHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
}

ggh_samples = {
    'ggh125_powheg': ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'Z'],
}

sig_samples = {
    'ggh125_madgraph_inc_a1-prod_nom-decay': ['GluGluToHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_two_a1-prod_nom-decay_v1': ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
    'ggh125_madgraph_two_a1-prod_nom-decay_v2': ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_two_a3int-prod_nom-decay': ['GluGluToMaxmixHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_two_a3-prod_nom-decay': ['GluGluToPseudoscalarHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_minlo_CP-down': ['GluGluToHToTauTauPlusTwoJets_M125_TuneCP5Down_PSweights_13TeV_powheg-minlo_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'ggh125_minlo_CP-nom': ['GluGluToHToTauTauPlusTwoJets_M125_TuneCP5_PSweights_13TeV_powheg-minlo_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'ggh125_minlo_CP-up': ['GluGluToHToTauTauPlusTwoJets_M125_TuneCP5Up_PSweights_13TeV_powheg-minlo_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_l1int-prod_nom-decay': ['VBFHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_l1-prod_nom-decay': ['VBFHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_l1zgint-prod_nom-decay': ['VBFHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_l1zg-prod_nom-decay': ['VBFHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a3int-prod_nom-decay': ['VBFHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a3-prod_nom-decay': ['VBFHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a2int-prod_nom-decay': ['VBFHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a2-prod_nom-decay': ['VBFHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a1-prod_nom-decay': ['VBFHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'wh125_JHU_l1int-prod_nom-decay': ['WHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_l1-prod_nom-decay': ['WHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a3int-prod_nom-decay': ['WHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a3-prod_nom-decay': ['WHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', '0'],
    'wh125_JHU_a2int-prod_nom-decay': ['WHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', '0'],
    'wh125_JHU_a2-prod_nom-decay': ['WHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a1-prod_nom-decay': ['WHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_l1int-prod_nom-decay': ['ZHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_l1-prod_nom-decay': ['ZHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_l1zgint-prod_nom-decay': ['ZHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_l1zg-prod_nom-decay': ['ZHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4', '0'],
    'zh125_JHU_a3int-prod_nom-decay': ['ZHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a3-prod_nom-decay': ['ZHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a2int-prod_nom-decay': ['ZHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a2-prod_nom-decay': ['ZHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a1-prod_nom-decay': ['ZHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
}

data_et_samples = {
    'datasE-B': ['data_SingleElectron_Run2017B-31Mar2018', '0'],
    'datasE-C': ['data_SingleElectron_Run2017C-31Mar2018', '0'],
    'datasE-D': ['data_SingleElectron_Run2017D-31Mar2018', '0'],
    'datasE-E': ['data_SingleElectron_Run2017E-31Mar2018', '0'],
    'datasE-F': ['data_SingleElectron_Run2017F-31Mar2018', '0'],
}

data_mt_samples = {
    'datasMu-B': ['data_SingleMuon_Run2017B-31Mar2018', '0'],
    'datasMu-C': ['data_SingleMuon_Run2017C-31Mar2018', '0'],
    'datasMu-D': ['data_SingleMuon_Run2017D-31Mar2018', '0'],
    'datasMu-E': ['data_SingleMuon_Run2017E-31Mar2018', '0'],
    'datasMu-F': ['data_SingleMuon_Run2017F-31Mar2018', '0'],
}

embed_et_samples = {
    'embedEl-B': ['embedded_EmbeddingRun2017B_ElTauFinalState', '0'],
    'embedEl-C': ['embedded_EmbeddingRun2017C_ElTauFinalState', '0'],
    'embedEl-D': ['embedded_EmbeddingRun2017D_ElTauFinalState', '0'],
    'embedEl-E': ['embedded_EmbeddingRun2017E_ElTauFinalState', '0'],
    'embedEl-F': ['embedded_EmbeddingRun2017F_ElTauFinalState', '0'],
}

embed_mt_samples = {
    'embedMu-B': ['embedded_EmbeddingRun2017B_MuTauFinalState', '0'],
    'embedMu-C': ['embedded_EmbeddingRun2017C_MuTauFinalState', '0'],
    'embedMu-D': ['embedded_EmbeddingRun2017D_MuTauFinalState', '0'],
    'embedMu-E': ['embedded_EmbeddingRun2017E_MuTauFinalState', '0'],
    'embedMu-F': ['embedded_EmbeddingRun2017F_MuTauFinalState', '0'],
}

prefix = args.prefix
jobType = args.job

ggh_pref = "/hdfs/store/user/tmitchel/SMHTT_2017_legacy_mc_v3p2_ggh-only/"
sig_pref = '/hdfs/store/user/abdollah/FSA_2017_AC_XMass/'
bkg_mt_pref = "/hdfs/store/user/tmitchel/SMHTT_2017_legacy_mc_v3p2_mt/"
bkg_et_pref = "/hdfs/store/user/tmitchel/SMHTT_2017_legacy_mc_v3p2_et/"
embed_pref = "/hdfs/store/user/tmitchel/SMHTT_2017_legacy_embedded_v3/"
data_pref = "/hdfs/store/user/tmitchel/SMHTT_2017_legacy_data_v3/"

settings = {
#   'signal_et': [sig_et_pref, sig_samples],
  'ggh': [ggh_pref, ggh_samples],
  'sig': [sig_pref, sig_samples],
  'bkg_mt': [bkg_mt_pref, bkg_samples],
  'bkg_et': [bkg_et_pref, bkg_samples],
  'dataMu': [data_pref, data_mt_samples],
  'dataEl': [data_pref, data_et_samples],
  'embedMu': [embed_pref, embed_mt_samples],
  'embedEl': [embed_pref, embed_et_samples],
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

    # base command
    command = '$CMSSW_BASE/bin/$SCRAM_ARCH/uniSkim -d {} -j {} -r {} -y {} -l {}'.format(
        pref+path, jobType, recoil, '2017', lep)
    
    # is this signal
    if job == 'sig' or job == 'ggh':
        command += ' -s '

    # boilerplate
    command +=  '-i input_file.root -o \'$OUTPUT\''

    ch.submit_command(command, prefix, pref+path, sample, use_input='-n ', dryrun=False)
