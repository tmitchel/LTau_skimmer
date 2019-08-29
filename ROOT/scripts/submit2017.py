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

bkg_samples_batch1 = {
    # 'DYJets1_v1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
    'DYJets1': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', 'Z'],
    # 'DYJets1_old_pmx_v2': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_old_pmx_94X_mc2017_realistic_v14-v2', 'Z'],
    # 'DYJets1_v3': ['DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_v3_94X_mc2017_realistic_v14_ext1-v2', 'Z'],
    # 'DYJets2_v1': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    # 'DYJets2_ext_v1': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
    'DYJets2': ['DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v2', 'Z'],
    # 'DYJets3_v1': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    # 'DYJets3_ext_v1': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
    'DYJets3': ['DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14_ext1-v2', 'Z'],
    'DYJets4_v1': ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'DYJets4_v2': ['DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_v2_94X_mc2017_realistic_v14-v2', 'Z'],
    # 'DYJets_lowM_v1': ['DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    # 'DYJets_lowM_ext_v1': ['DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
    # 'DYJets_lowM_ext_v2': ['DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2', 'Z'],
    'DYJets_v1': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'DYJets_v2': ['DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],

    'WJets1_v1': ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', 'W'],
    'WJets1_v2': ['W1JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4', 'W'],
    'WJets2_v1': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4', 'W'],
    'WJets2_v2': ['W2JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v5', 'W'],
    'WJets3': ['W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'W'],
    # 'WJets4_v1': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'W'],
    'WJets4': ['W4JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'W'],
    'WJets_v1': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'W'],
    'WJets_v2': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', 'W'],
    'WJets_v3': ['WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2', 'W'],
}
bkg_samples_batch2 = {
    'EWKWMinus': ['EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'EWKWMinus': ['EWKWMinus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],
    'EWKWPlus': ['EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'EWKWPlus': ['EWKWPlus2Jets_WToLNu_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],
    # 'EWKZ2l_v1': ['EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'EWKZ2l': ['EWKZ2Jets_ZToLL_M-50_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],
    # 'EWKZ2nu_v1': ['EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'EWKZ2nu': ['EWKZ2Jets_ZToNuNu_TuneCP5_13TeV-madgraph-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],

    'Tbar-tchan_v1': ['ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'Tbar-tchan_v2': ['ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    # 'T-tchan_v1': ['ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'T-tchan': ['ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', '0'],
    'Tbar-tW': ['ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'T-tW_v1': ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'T-tW_v2': ['ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    # 'TTLep_v1': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    # 'TTLep_v2': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'TTLep': ['TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', '0'],
    # 'TTHad_v1': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'TTHad': ['TTToHadronic_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],
    # 'TTSemi_v2': ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'TTSemi': ['TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', '0'],
    'WW_v1': ['WW_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'WW_v2': ['WW_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'WZ': ['WZ_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    # 'ZZ_v1': ['ZZ_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'ZZ': ['ZZ_TuneCP5_13TeV-pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', '0'],
}

sig_samples_p1 = {
    'ggh125_powheg_v1': ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_powheg_v2': ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
    'ggh125_madgraph_inc_maxmix_decay': ['GluGluToHToTauTauMaxmixDecay_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_two_a1-prod_nom-decay_v1': ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_two_a1-prod_nom-decay_v2': ['GluGluToHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1', 'Z'],
    'ggh125_madgraph_two_a3-prod_nom-decay': ['GluGluToPseudoscalarHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_inc_a1-prod_ps-decay': ['GluGluToHToTauTauPseudoscalarDecay_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_inc_a3int-prod_nom-decay': ['GluGluToMaxmixHToTauTauPseudoscalarDecay_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'ggh125_madgraph_inc_a3-prod_nom-decay': ['GluGluToPseudoscalarHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_JHU_a3-prod_ps-decay': ['JJHiggs0MToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'ggh125_JHU_a3-prod_nom-decay': ['JJHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'ggh125_JHU_a3int-prod_nom-decay': ['JJHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'ggh125_JHU_a1-prod_nom-decay': ['JJHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_madgraph': ['VBFHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'vbf125_JHU_l1-prod_nom-decay': ['VBFHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_l1int-prod_nom-decay': ['VBFHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a3int-prod_nom-decay': ['VBFHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'wh125_JHU_l1-prod_nom-decay': ['WHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a2-prod_nom-decay': ['WHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a2int-prod_nom-decay': ['WHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', '0'],
    'wh125_JHU_a1-prod_maxmix-decay': ['WHiggs0PMToTauTauMaxmixDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a1-prod_ps-decay': ['WHiggs0PMToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', '0'],
    'wh125_JHU_a1-prod_nom-decay': ['WHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wplus125_powheg': ['WplusHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'zh125_JHU_l1-prod_nom-decay': ['ZHiggs0L1ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_l1zg-prod_nom-decay': ['ZHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v4', '0'],
    'zh125_JHU_l1zgint-prod_nom-decay': ['ZHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a3-prod_ps-decay': ['ZHiggs0MToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a3-prod_nom-decay': ['ZHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a2-prod_nom-decay': ['ZHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a2int-prod_nom-decay': ['ZHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a1-prod_maxmix-decay': ['ZHiggs0PMToTauTauMaxmixDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a1-prod_ps-decay': ['ZHiggs0PMToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', '0'],
    }

sig_samples_p2 = {
    'ggh125_powheg_v3': ['GluGluHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_powheg_v1': ['VBFHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'vbf125_powheg_v2': ['VBFHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1', 'Z'],
    'vbf125_JHU_l1zg-prod_nom-decay': ['VBFHiggs0L1ZgToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_l1zgint-prod_nom-decay': ['VBFHiggs0L1Zgf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a3-prod_ps-decay': ['VBFHiggs0MToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a3-prod_nom-decay': ['VBFHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a2-prod_nom-decay': ['VBFHiggs0PHToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a2int-prod_nom-decay': ['VBFHiggs0PHf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a1-prod_nom-decay': ['VBFHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'wh125_JHU_l1int-prod_nom-decay': ['WHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a3-prod_ps-decay': ['WHiggs0MToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wh125_JHU_a3-prod_nom-decay': ['WHiggs0MToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v3', '0'],
    'wh125_JHU_a3int-prod_nom-decay': ['WHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'wminus125_powheg': ['WminusHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'zh125_powheg': ['ZHToTauTau_M125_13TeV_powheg_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', '0'],
    'zh125_JHU_l1int-prod_nom-decay': ['ZHiggs0L1f05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a3int-prod_nom-decay': ['ZHiggs0Mf05ph0ToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
    'zh125_JHU_a1-prod_nom-decay': ['ZHiggs0PMToTauTau_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', '0'],
}

sig_samples_p3 = {
    'ggh125_madgraph_inc_nominal': ['GluGluToHToTauTau_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'],
    'ggh125_madgraph_two_a3int-prod_nom-decay': ['GluGluToMaxmixHToTauTauPlusTwoJets_M125_13TeV_amcatnloFXFX_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1', 'Z'], 
    'ggh125_JHU_a1-prod_maxmix-decay': ['JJHiggs0PMToTauTauMaxmixDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'ggh125_JHU_a1-prod_ps-decay': ['JJHiggs0PMToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a1-prod_maxmix-decay': ['VBFHiggs0PMToTauTauMaxmixDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
    'vbf125_JHU_a1-prod_ps-decay': ['VBFHiggs0PMToTauTauPseudoscalarDecay_M125_13TeV_JHUGenV7011_pythia8_v2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2', 'Z'],
}

el_data_samples = {
    'datasE-B': ['data_SingleElectron_Run2017B-31Mar2018', '0'],
    'datasE-C': ['data_SingleElectron_Run2017C-31Mar2018', '0'],
    'datasE-D': ['data_SingleElectron_Run2017D-31Mar2018', '0'],
    'datasE-E': ['data_SingleElectron_Run2017E-31Mar2018', '0'],
    'datasE-F': ['data_SingleElectron_Run2017F-31Mar2018', '0'],
}

mu_data_samples = {
    'datasMu-B': ['data_SingleMuon_Run2017B-31Mar2018', '0'],
    'datasMu-C': ['data_SingleMuon_Run2017C-31Mar2018', '0'],
    'datasMu-D': ['data_SingleMuon_Run2017D-31Mar2018', '0'],
    'datasMu-E': ['data_SingleMuon_Run2017E-31Mar2018', '0'],
    'datasMu-F': ['data_SingleMuon_Run2017F-31Mar2018', '0'],
}

el_embedded_samples = {
    'embedEl-B': ['embedded_EmbeddingRun2017B_ElTauFinalState', '0'],
    'embedEl-C': ['embedded_EmbeddingRun2017C_ElTauFinalState', '0'],
    'embedEl-D': ['embedded_EmbeddingRun2017D_ElTauFinalState', '0'],
    'embedEl-E': ['embedded_EmbeddingRun2017E_ElTauFinalState', '0'],
    'embedEl-F': ['embedded_EmbeddingRun2017F_ElTauFinalState', '0'],
}

mu_embedded_samples = {
    'embedMu-B': ['embedded_EmbeddingRun2017B_MuTauFinalState', '0'],
    'embedMu-C': ['embedded_EmbeddingRun2017C_MuTauFinalState', '0'],
    'embedMu-D': ['embedded_EmbeddingRun2017D_MuTauFinalState', '0'],
    'embedMu-E': ['embedded_EmbeddingRun2017E_MuTauFinalState', '0'],
    'embedMu-F': ['embedded_EmbeddingRun2017F_MuTauFinalState', '0'],
}

prefix = args.prefix
jobType = args.job

bkg_pref = '/hdfs/store/user/tmitchel/SMHTT_2017_legacy_mc_v1/'
data_pref = '/hdfs/store/user/tmitchel/SMHTT_2017_legacy_data_v1/'
embed_pref = '/hdfs/store/user/tmitchel/SMHTT_2017_embedded_v1/'
sig_p1_pref = '/hdfs/store/user/senka/SMHTT_2017_Ntuples_valid_part1/'
sig_p2_pref = '/hdfs/store/user/senka/SMHTT_2017_Ntuples_valid_part2/'
sig_p3_pref = '/hdfs/store/user/senka/SMHTT_2017_Ntuples_valid_part3/'

settings = {
    'sig_p1': [sig_p1_pref, sig_samples_p1],
    'sig_p2': [sig_p2_pref, sig_samples_p2],
    'sig_p3': [sig_p3_pref, sig_samples_p3],

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
    pref = settings[args.job][0]

    command = '$CMSSW_BASE/bin/$SCRAM_ARCH/uniSkim -d %s -j %s -r %s -y %s -l %s -i input_file.root -o \'$OUTPUT\'' % (
        pref+path, jobType, recoil, '2017', lep)
    ch.submit_command(command, prefix, pref+path, sample, use_input='-n ', dryrun=False)
