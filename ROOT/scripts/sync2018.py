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

sync = {
    'vbf125_powheg': ['VBFHToTauTau_M125_13TeV_powheg_pythia8_-102X_upgrade2018_realistic_v15_ext1-v1', 'Z'],
}

data_mt_samples = {
    'datasMu-A': ['data_SingleMuon_Run2018A-17Sep2018', '0'],
    'datasMu-B': ['data_SingleMuon_Run2018B-17Sep2018', '0'],
    'datasMu-C': ['data_SingleMuon_Run2018C-17Sep2018', '0'],
    'datasMu-D': ['data_SingleMuon_Run2018D-PromptReco', '0'],
}

prefix = args.prefix
jobType = args.job

bkg_pref1 = '/hdfs/store/user/caillol/SMHTT_2018_20nov_mc/'
bkg_pref2 = '/hdfs/store/user/caillol/SMHTT_2018_20nov_highMem_etmt_mc/'
data_pref = '/hdfs/store/user/caillol/SMHTT_2018_31oct_data/'
embed_pref = '/hdfs/store/user/caillol/SMHTT_2018_20nov_embedded/'
sig_pref = '/hdfs/store/user/abdollah/FSA_2018_AC_XMass/'

settings = {
  'sync': [bkg_pref1, sync],
  'dataMu': [data_pref, data_mt_samples],
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
        pref+path, jobType, recoil, '2018', 'sync')
    ch.submit_command(command, prefix, pref+path, sample, use_input='-n ', dryrun=False)
