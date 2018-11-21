import os, sys, re

#datasets = "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets"

#datasets = "/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017RECOSIMstep_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets_ext1"

#datasets = "/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets1"

#datasets = "/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets1_ext1"

#datasets = "/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets2"

#datasets = "/DY2JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets2_ext"

#datasets = "/DY3JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets3"

#datasets = "/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_v2_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets4_real"

#datasets = "/DY4JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/DYJets4"

#datasets = "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/TTLep"

#datasets = "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/TTHad_v2"

#datasets = "/TTToHadronic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/TTHad"

#datasets = "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/TTSemi"

#datasets = "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/TTSemi_v2"

#datasets = "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/WJets"

#datasets = "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/WJets_ext1"

#datasets = "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/Tbar-tW"

#datasets = "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/T-tW_v2"

#datasets = "/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/T-tchan"

#datasets = "/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/Tbar-tchan_v2"

#datasets = "/WWTo4Q_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#test = "/hdfs/store/user/kkaadze/et_111418/WW4q"

#dir = "/hdfs/store/user/kkaadze/et_111418/"

#datasets = "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"
#dirName = "WW1l1nu2q_ext"

#datasets = "/WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#dirName = "WW1l1nu2q_amc"

#datasets = "/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#dirName = "WW1l1nu2q"

#datasets = "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
#dirName = "WZ1l1nu2q_v2"

#datasets = "/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_old_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#dirName = "WZ1l1nu2q"

#datasets = "/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#dirName = "WZ1l3nu"

#datasets = "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#dirName = "ZZ4l"

#datasets = "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM"
#dirName = "ZZ4l_ext1"

#datasets = "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#datasets = "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#datasets = "/ZZTo4L_13TeV_powheg_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
#dirName = "ZZ4l_v2"

dir = "/hdfs/store/user/caillol/SMHTT_2017_7nov/"
dirName = "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8_v14-v1"
#datasets = "/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM"
datasets = "/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/MINIAODSIM"
test = dir + dirName

data = datasets
output = os.popen("dasgoclient --query=\"file dataset=" + data + "\"")
files = output.readlines()

numbers = []
for file in files:
    elements = re.split("/", file)
    element = elements[len(elements) -1]
    numbers += [element[:-6],]

output = os.popen("ls " + test)
localFiles = output.readlines()

elements = []
for line in localFiles:
    elements += [re.split("make_ntuples_cfg-", line)[1][:-6],]

print "Local files:", len(elements), "in dbs:", len(numbers)

# ===================================

counter = 0
for el in numbers:
    if not el in elements:
        counter += 1
        print counter, el, "is not present locally"

counter = 0
for el in elements:
    if not el in numbers:
        counter += 1
        print counter, el, "is not in dataset"
