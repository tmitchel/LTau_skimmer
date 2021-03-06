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

bkg_samples_batch1 = {
    'DYJets1_v1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets1_v2': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'DYJets2_v1': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets2_v2': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'DYJets3_v1': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets3_v2': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'DYJets4': ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'DYJets': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],

    'WJets1': ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets2_v1': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'WJets2_v2': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets3_v1': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'WJets3_v2': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets4_v1': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'WJets4_v2': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],
    'WJets': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_-102X_upgrade2018_realistic_v15-v2', 'W'],

    'TTHad_v1': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'TTHad_v2': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext2-v2', '0'],

    'EWKWMinus': ['EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'EWKWPlus': ['EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', 'W'],
    'EWKZ2l': ['EWKZ2Jets_ZToLL_M-50_TuneCP5_PSweights_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'EWKZ2nu': ['EWKZ2Jets_ZToNuNu_TuneCP5_PSweights_13TeV-madgraph-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'Tbar-tchan': ['ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'T-tchan': ['ST_t-channel_top_5f_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'Tbar-tW': ['ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext1-v1', '0'],
    'T-tW': ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext1-v1', '0'],
    'WW_v1': ['WW_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WW_v2': ['WW_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'WZ': ['WZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v3', '0'],
    'ZZ': ['ZZ_TuneCP5_13TeV-pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],

    'VV': ['VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WW1l1nu2q': ['WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WZ1l1nu2q': ['WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WZ1l3nu': ['WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WZ2l2q': ['WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WZ3l1nu_v1': ['WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'WZ3l1nu_v2': ['WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15_ext1-v2', '0'],
    'ZZ2l2q': ['ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'ZZ4l': ['ZZTo4L_TuneCP5_13TeV-amcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],

    'ggh125_powheg': ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_powheg': ['VBFHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15_ext1-v1', 'Z'],
    'wminus125_powheg': ['WminusHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wplus125_powheg': ['WplusHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_powheg': ['ZHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
}

bkg_samples_batch2 = {
    'TTLep': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'TTSemi_v1': ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15-v1', '0'],
    'TTSemi_v2': ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_-102X_upgrade2018_realistic_v15_ext3-v2', '0'],
}

sig_samples = {
    'ggh125_minlo_DownPS': ['GluGluToHToTauTauPlusTwoJets_M125_TuneCP5Down_PSweights_13TeV_powheg-minlo_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'ggh125_minlo': ['GluGluToHToTauTauPlusTwoJets_M125_TuneCP5_PSweights_13TeV_powheg-minlo_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'ggh125_minlo_UpPS': ['GluGluToHToTauTauPlusTwoJets_M125_TuneCP5Up_PSweights_13TeV_powheg-minlo_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'ggh125_madgraph_one_a3int_filtered': ['JJH0Mf05ph0ToTauTauPlusOneJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_mdagraph_one_a3int_unfiltered': ['JJH0Mf05ph0ToTauTauPlusOneJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_two_a3int_filtered': ['JJH0Mf05ph0ToTauTauPlusTwoJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_madgraph_two_a3int_unfiltered': ['JJH0Mf05ph0ToTauTauPlusTwoJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_zero_a3int_filtered': ['JJH0Mf05ph0ToTauTauPlusZeroJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_madgraph_zero_a3int_unfiltered': ['JJH0Mf05ph0ToTauTauPlusZeroJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_one_a3_filtered': ['JJH0MToTauTauPlusOneJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_mdagraph_one_a3_unfiltered': ['JJH0MToTauTauPlusOneJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_two_a3_filtered': ['JJH0MToTauTauPlusTwoJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_madgraph_two_a3_unfiltered': ['JJH0MToTauTauPlusTwoJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_zero_a3_filtered': ['JJH0MToTauTauPlusZeroJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_madgraph_zero_a3_unfiltered': ['JJH0MToTauTauPlusZeroJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_one_a1_filtered': ['JJH0PMToTauTauPlusOneJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_mdagraph_one_a1_unfiltered': ['JJH0PMToTauTauPlusOneJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_two_a1_filtered': ['JJH0PMToTauTauPlusTwoJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_madgraph_two_a1_unfiltered': ['JJH0PMToTauTauPlusTwoJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    'ggh125_madgraph_zero_a1_filtered': ['JJH0PMToTauTauPlusZeroJets_Filtered_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    # 'ggh125_madgraph_zero_a1_unfiltered': ['JJH0PMToTauTauPlusZeroJets_M125_TuneCP5_13TeV-mcatnloFXFX-pythia8_-102X_upgrade2018_realistic_v15-v1', 'Z'],
    
    'vbf125_JHU_l1-prod_nom-decay': ['VBFHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_l1zg-prod_nom-decay': ['VBFHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_l1zgint-prod_nom-decay': ['VBFHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_l1int-prod_nom-decay': ['VBFHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_a3-prod_nom-decay': ['VBFHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_a3int-prod_nom-decay': ['VBFHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_a2-prod_nom-decay': ['VBFHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_a2int-prod_nom-decay': ['VBFHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'vbf125_JHU_a1-prod_nom-decay': ['VBFHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', 'Z'],
    'wh125_JHU_l1-prod_nom-decay': ['WHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wh125_JHU_l1int-prod_nom-decay': ['WHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wh125_JHU_a3-prod_nom-decay': ['WHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wh125_JHU_a3int-prod_nom-decay': ['WHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wh125_JHU_a2-prod_nom-decay': ['WHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wh125_JHU_a2int-prod_nom-decay': ['WHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'wh125_JHU_a1-prod_nom-decay': ['WHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_l1-prod_nom-decay': ['ZHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_l1zg-prod_nom-decay': ['ZHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_l1int-prod_nom-decay': ['ZHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_l1zgint-prod_nom-decay': ['ZHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_a3-prod_nom-decay': ['ZHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_a3int-prod_nom-decay': ['ZHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_a2-prod_nom-decay': ['ZHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_a2int-prod_nom-decay': ['ZHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
    'zh125_JHU_a1-prod_nom-decay': ['ZHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_-102X_upgrade2018_realistic_v15-v2', '0'],
}

data_et_samples = {
    'datasE-A': ['data_EGamma_Run2018A-17Sep2018', '0'],
    'datasE-B': ['data_EGamma_Run2018B-17Sep2018', '0'],
    'datasE-C': ['data_EGamma_Run2018C-17Sep2018', '0'],
    'datasE-D': ['data_EGamma_Run2018D-PromptReco', '0'],
}

data_mt_samples = {
    'datasMu-A': ['data_SingleMuon_Run2018A-17Sep2018', '0'],
    'datasMu-B': ['data_SingleMuon_Run2018B-17Sep2018', '0'],
    'datasMu-C': ['data_SingleMuon_Run2018C-17Sep2018', '0'],
    'datasMu-D': ['data_SingleMuon_Run2018D-PromptReco', '0'],
}

embed_et_samples = {
    'embedEl-A': ['embedded_EmbeddingRun2018A_ElTauFinalState', '0'],
    'embedEl-B': ['embedded_EmbeddingRun2018B_ElTauFinalState', '0'],
    'embedEl-C': ['embedded_EmbeddingRun2018C_ElTauFinalState', '0'],
    'embedEl-D': ['embedded_EmbeddingRun2018D_ElTauFinalState', '0'],
}

embed_mt_samples = {
    'embedMu-A': ['embedded_EmbeddingRun2018A_MuTauFinalState', '0'],
    'embedMu-B': ['embedded_EmbeddingRun2018B_MuTauFinalState', '0'],
    'embedMu-C': ['embedded_EmbeddingRun2018C_MuTauFinalState', '0'],
    'embedMu-D': ['embedded_EmbeddingRun2018D_MuTauFinalState', '0'],
}


prefix = args.prefix
jobType = args.job

bkg_pref1 = '/hdfs/store/user/caillol/SMHTT_2018_20nov_mc/'
bkg_pref2 = '/hdfs/store/user/caillol/SMHTT_2018_20nov_highMem_etmt_mc/'
data_pref = '/hdfs/store/user/caillol/SMHTT_2018_31oct_data/'
embed_pref = '/hdfs/store/user/caillol/SMHTT_2018_20nov_embedded/'
sig_pref = '/hdfs/store/user/abdollah/FSA_2018_AC_XMass/'

settings = {
  'sig': [sig_pref, sig_samples],
  'bkg1': [bkg_pref1, bkg_samples_batch1],
  'bkg2': [bkg_pref2, bkg_samples_batch2],
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
    command = '$CMSSW_BASE/bin/$SCRAM_ARCH/uniSkim -d {} -j {} -r {} -y {} -l {} '.format(
        pref+path, jobType, recoil, '2018', lep)
    
    # is this signal
    if args.job == 'sig' or args.job == 'ggh':
        command += ' -s '

    # boilerplate
    command +=  ' -i input_file.root -o \'$OUTPUT\''

    ch.submit_command(command, prefix, pref+path, sample, use_input='-n ', dryrun=False)
