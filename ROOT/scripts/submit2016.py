import os
import subprocess
import condor_handler as ch
from argparse import ArgumentParser

parser = ArgumentParser(
    description="run on a directory containing directories containing the skimmed ntuples")
parser.add_argument('-p', '--prefix', action='store',
                    help='name to prefix all directories')
parser.add_argument('-l', '--lepton', action='store',
                    help='which lepton (mt or et)')
parser.add_argument('-j', '--job', action='store',
                    help='job type [dataMu, dataEl, sig, ggh, bkg1, bkg2, embedMu, embedEl]')
args = parser.parse_args()

test_batch = {
    "DYJets1": ["DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6-v1", 'Z'],
}

bkg_samples_batch1 = {
    "DYJets1": ["DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1", 'Z'],
    "DYJets2": ["DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", 'Z'],
    "DYJets3": ["DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", 'Z'],
    "DYJets4": ["DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", 'Z'],
    "DYJets_v1": ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", 'Z'],
    "DYJets_v2": ["DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2", 'Z'],

    "WJets1": ["W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1", 'W'],
    "WJets2_v1": ["W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", 'W'],
    "WJets2_v2": ["W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", 'W'],
    "WJets3_v1": ["W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", 'W'],
    "WJets3_v2": ["W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", 'W'],
    "WJets4_v1": ["W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", 'W'],
    "WJets4_v2": ["W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v1", 'W'],
    "WJets4_v3": ["W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", 'W'],
    "WJets_v1": ["WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2",  'W'],
    "WJets_v2": ["WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2",  'W'],

    "EWKWMinus_v1": ["EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2",  'W'],
    "EWKWMinus_v2": ["EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2",  'W'],
    "EWKWMinus_v3": ["EWKWMinus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2",  'W'],
    "EWKWPlus_v1": ["EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2",  'W'],
    "EWKWPlus_v2": ["EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2",  'W'],
    "EWKWPlus_v3": ["EWKWPlus2Jets_WToLNu_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2",  'W'],
    "EWKZ2l_v1": ["EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2",  'Z'],
    "EWKZ2l_v2": ["EWKZ2Jets_ZToLL_M-50_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2",  'Z'],
    "EWKZ2nu_v1": ["EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2",  'Z'],
    "EWKZ2nu_v2": ["EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2",  'Z'],
    "EWKZ2nu_v3": ["EWKZ2Jets_ZToNuNu_13TeV-madgraph-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext2-v2",  'Z'],

    "TT": ["TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1", '0'],
    "TT_evtgen": ['TT_TuneCUETP8M2T4_13TeV-powheg-pythia8-evtgen_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', '0'],
    "Tbar-tchan": ["ST_t-channel_antitop_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1", '0'],
    "T-tchan": ["ST_t-channel_top_4f_inclusiveDecays_13TeV-powhegV2-madspin-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1", '0'],
    "Tbar-tW": ["ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1", '0'],
    "T-tW": ["ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v1", '0'],

    "WW_v1": ["WW_TuneCUETP8M1_13TeV-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", "0"],
    "WW_v2": ["WW_TuneCUETP8M1_13TeV-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", "0"],
    "WZ_v1": ["WZ_TuneCUETP8M1_13TeV-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", "0"],
    "WZ_v2": ["WZ_TuneCUETP8M1_13TeV-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", "0"],
    "ZZ_v1": ["ZZ_TuneCUETP8M1_13TeV-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", "0"],
    "ZZ_v2": ["ZZ_TuneCUETP8M1_13TeV-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", "0"],

    "VV_v1": ["VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", '0'],
    "VV_v2": ["VVTo2L2Nu_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", '0'],
    "WW1l1nu2q": ["WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", '0'],
    "WZ1l1nu2q": ["WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", '0'],
    "WZ1l3nu": ["WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", '0'],
    "WZ2l2q": ["WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", '0'],
    "WZ3l1nu": ["WZTo3LNu_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2", '0'],
    "ZZ2l2q": ["ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1", '0'],
    "ZZ4l": ["ZZTo4L_13TeV-amcatnloFXFX-pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2", '0'],
 
    'vbf125_powheg' : ['VBFHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
    'wplus125_powheg' : ['WplusHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
    'wminus125_powheg' : ['WminusHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
    'zh125_powheg' : ['ZHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
}

ggh_samples = {
    'ggh125_powheg_v1' : ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
    'ggh125_powheg_v2' : ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
    'ggh125_powheg_v3' : ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v3', 'Z'],
}

sig_samples = {
   'ggh125_minlo_nnlops': ['GluGluHToTauTau_M-125_13TeV_powheg_MINLO_NNLOPS_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
   # 'ggh125_JHU_a1-prod_nom-decay' : ['GluGluH2JetsToTauTau_M125_13TeV_CPmixing_sm_JHU_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2', 'Z'],
   # 'ggh125_JHU_a3-prod_nom-decay' : ['GluGluH2JetsToTauTau_M125_13TeV_CPmixing_pseudoscalar_JHU_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext3-v2', 'Z'],
   # 'ggh125_JHU_a3int-prod_nom-decay' : ['GluGluH2JetsToTauTau_M125_13TeV_CPmixing_maxmix_JHU_v3-PUMoriond17_94X_mcRun2_asymptotic_v3_ext1-v2', 'Z'],
   'ggh125_minlo' : ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_powheg-minlo_pythia8_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_minlo_CUEDown' : ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_powheg-minlo_pythia8_CUETP8M1Down_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_minlo_CUEUp' : ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_powheg-minlo_pythia8_CUETP8M1Up_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_minlo_DownPS' : ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_powheg-minlo_pythia8_DownPS_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_minlo_UpPS' : ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_powheg-minlo_pythia8_UpPS_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_madgraph_inc_a1-prod_nom-decay' : ['GluGluToHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_tauola_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_madgraph_inc_a3-prod_nom-decay' : ['GluGluToPseudoscalarHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_tauola_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_madgraph_inc_a3int-prod_nom-decay' : ['GluGluToMaxmixHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_tauola_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_madgraph_two_a1-prod_nom-decay' : ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_tauola_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_madgraph_two_a3-prod_nom-decay' : ['GluGluToPseudoscalarHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_tauola_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],
   'ggh125_madgraph_two_a3int-prod_nom-decay' : ['GluGluToMaxmixHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_tauola_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', 'Z'],

   'vbf125_JHU_a1-prod_nom-decay' : ['VBFHiggs0PM_HToTauTau_M-125_13TeV-JHUGenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
   'vbf125_JHU_a2-prod_nom-decay' : ['VBFHiggs0PH_HToTauTau_M-125_13TeV-JHUGenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
   'vbf125_JHU_a3-prod_nom-decay' : ['VBFHiggs0M_HToTauTau_M-125_13TeV-JHUGenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
   'vbf125_JHU_a2int-prod_nom-decay' : ['VBFHiggs0PHf05ph0_HToTauTau_M-125_13TeV-JHUGenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
   'vbf125_JHU_a3int-prod_nom-decay' : ['VBFHiggs0Mf05ph0_HToTauTau_M-125_13TeV-JHUGenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],
   'vbf125_JHU_l1-prod_nom-decay': ['VBFHiggs0L1_HToTauTau_M-125_13TeV-JHUGenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v4', 'Z'],
   'vbf125_JHU_l1int-prod_nom-decay' : ['VBFHiggs0L1f05ph0_HToTauTau_M-125_13TeV-JHUGenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', 'Z'],

   'wh125_JHU_a1-prod_nom-decay' : ['WHiggs0PM_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'wh125_JHU_a2-prod_nom-decay' : ['WHiggs0PH_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'wh125_JHU_a3-prod_nom-decay' : ['WHiggs0M_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'wh125_JHU_a2int-prod_nom-decay' : ['WHiggs0PHfWH05ph0_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'wh125_JHU_a3int-prod_nom-decay' : ['WHiggs0MfWH05ph0_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'wh125_JHU_l1-prod_nom-decay': ['WHiggs0L1_HToTauTau_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'wh125_JHU_l1int-prod_nom-decay' : ['WHiggs0L1fWH05ph0_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   
   'zh125_JHU_a1-prod_nom-decay' : ['ZHiggs0PM_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'zh125_JHU_a2-prod_nom-decay' : ['ZHiggs0PH_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'zh125_JHU_a2int-prod_nom-decay': ['ZHiggs0PHfZH05ph0_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'zh125_JHU_a3-prod_nom-decay' : ['ZHiggs0M_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'zh125_JHU_a3int-prod_nom-decay' : ['ZHiggs0MfZH05ph0_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v1', '0'],
   'zh125_JHU_l1-prod_nom-decay' : ['ZHiggs0L1_Undecayed_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
   'zh125_JHU_l1int-prod_nom-decay' : ['ZHiggs0L1fZH05ph0_HToTauTau_M-125_13TeV-JHUgenV6_v3-PUMoriond17_94X_mcRun2_asymptotic_v3-v2', '0'],
}

data_et_samples = {
    "datasE-B_v1": ["data_SingleElectron_Run2016B_v1", '0'],
    "datasE-B_v2": ["data_SingleElectron_Run2016B_v2", '0'],
    "datasE-C": ["data_SingleElectron_Run2016C", '0'],
    "datasE-D": ["data_SingleElectron_Run2016D", '0'],
    "datasE-E": ["data_SingleElectron_Run2016E", '0'],
    "datasE-F": ["data_SingleElectron_Run2016F", '0'],
    "datasE-G": ["data_SingleElectron_Run2016G", '0'],
    "datasE-H": ["data_SingleElectron_Run2016H", '0'],
}

data_mt_samples = {
    "datasMu-B_v1": ["data_SingleMuon_Run2016B_v1", '0'],
    "datasMu-B_v2": ["data_SingleMuon_Run2016B_v2", '0'],
    "datasMu-C": ["data_SingleMuon_Run2016C", '0'],
    "datasMu-D": ["data_SingleMuon_Run2016D", '0'],
    "datasMu-E": ["data_SingleMuon_Run2016E", '0'],
    "datasMu-F": ["data_SingleMuon_Run2016F", '0'],
    "datasMu-G": ["data_SingleMuon_Run2016G", '0'],
    "datasMu-H": ["data_SingleMuon_Run2016H", '0'],
}

embed_et_samples = {
    'embedEl-B': ['embedded_EmbeddingRun2016B_ElTauFinalState', '0'],
    'embedEl-C': ['embedded_EmbeddingRun2016C_ElTauFinalState', '0'],
    'embedEl-D': ['embedded_EmbeddingRun2016D_ElTauFinalState', '0'],
    'embedEl-E': ['embedded_EmbeddingRun2016E_ElTauFinalState', '0'],
    'embedEl-F': ['embedded_EmbeddingRun2016F_ElTauFinalState', '0'],
    'embedEl-G': ['embedded_EmbeddingRun2016G_ElTauFinalState', '0'],
    'embedEl-H': ['embedded_EmbeddingRun2016H_ElTauFinalState', '0'],
}

embed_mt_samples = {
    'embedMu-B': ['embedded_EmbeddingRun2016B_MuTauFinalState', '0'],
    'embedMu-C': ['embedded_EmbeddingRun2016C_MuTauFinalState', '0'],
    'embedMu-D': ['embedded_EmbeddingRun2016D_MuTauFinalState', '0'],
    'embedMu-E': ['embedded_EmbeddingRun2016E_MuTauFinalState', '0'],
    'embedMu-F': ['embedded_EmbeddingRun2016F_MuTauFinalState', '0'],
    'embedMu-G': ['embedded_EmbeddingRun2016G_MuTauFinalState', '0'],
    'embedMu-H': ['embedded_EmbeddingRun2016H_MuTauFinalState', '0'],
}

prefix = args.prefix
jobType = args.job

bkg_pref = '/hdfs/store/user/aloeliger/SMHTT_2016_20nov/'
ggh_pref = '/hdfs/store/user/aloeliger/SMHTT_2016_20nov_ggH2/'
data_pref = '/hdfs/store/user/aloeliger/SMHTT_2016_data_18Sep/'
embed_pref = '/hdfs/store/user/aloeliger/SMHTT_2016_embedded_18sep/'
sig_pref = '/hdfs/store/user/abdollah/FSA_2016_AC_XMass/'

settings = {
  'sig': [sig_pref, sig_samples],
  'ggh': [ggh_pref, ggh_samples],
  'bkg_all': [bkg_pref, bkg_samples_batch1],
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
        pref+path, jobType, recoil, '2016', lep)
    
    # is this signal
    if args.job == 'sig' or args.job == 'ggh':
        command += ' -s '

    # boilerplate
    command +=  '-i input_file.root -o \'$OUTPUT\''

    ch.submit_command(command, prefix, pref+path, sample)
