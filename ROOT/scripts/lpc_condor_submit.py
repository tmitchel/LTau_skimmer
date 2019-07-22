import os
import sys
import pwd


def main(args):
    print "Begin submitting skims..."
    jobName = args.jobName
    sampledir = args.sampledir
    sample_name = os.path.basename(sampledir)
    print 'Processing samples from {} as {}'.format(sample_name, args.samplename)

    head_dir = '/uscmst1b_scratch/lpc1/3DayLifetime/{}/{}'.format(
        pwd.getpwuid(os.getuid())[0], jobName)
    if os.path.exists(head_dir):
        print 'Submission directory exists for {} {}.'.format(
            jobName, args.samplename)

    if sample_name == '':
        print "SAMPLE_NAME not defined, check for trailing '/' on sampledir path"
        return

    sample_dir = '{}/{}'.format(head_dir, args.samplename)

    exe_dir = '{}/executables'.format(sample_dir)
    os.system('mkdir -p {}'.format(exe_dir))

    config_dir = '{}/configs'.format(sample_dir)
    os.system('mkdir -p {}'.format(config_dir))

    fileList = [ifile for ifile in filter(None, os.popen(
        'xrdfs root://cmseos.fnal.gov/ ls '+sampledir).read().split('\n')) if '.root' in ifile]

    config_name = '{}/{}.jdl'.format(config_dir, args.samplename)
    condorConfig = '''universe = vanilla
Executable = overlord.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Output = sleep_\$(Cluster)_\$(Process).stdout
Error = sleep_\$(Cluster)_\$(Process).stderr
Log = sleep_\$(Cluster)_\$(Process).log
x509userproxy = $ENV(X509_USER_PROXY)
Queue {}
    '''.format(len(fileList))
    with open(config_name, 'w') as file:
        file.write(condorConfig)

    print 'Condor config has been written: {}'.format(config_name)

    overlord_name = '{}/{}.sh'.format(exe_dir, args.samplename)
    overloardScript = '''#!/bin/bash
let "sample=${{1}}"
bash {}_${{sample}}.sh
    '''.format(args.samplename)
    with open(overlord_name, 'w') as file:
        file.write(overloardScript)

    print 'Condor overlord has been written: {}'.format(overlord_name)


    bashScriptSetup = '''#!/bin/bash
source /cvmfs/cms.cern.ch/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
eval `scramv1 project CMSSW CMSSW_9_4_0`
cd CMSSW_9_4_0/src
eval `scramv1 runtime -sh`
xrdcp root://cmseos.fnal.gov//store/user/tmitchel/ltau_skimmer_packge.tar.gz .
tar xzf ltau_skimmer_packge.tar.gz
eval `scram b` \n
    '''

    i = 0
    for ifile in fileList:
        input_file = ifile
        output_file = '{}_{}.root'.format(args.samplename, i)

        # create the bash config script
        bash_name = '{}/{}.sh'.format(exe_dir, output_file.replace('.root', ''))
        bashScript = bashScriptSetup + '$CMSSW_BASE/bin/$SCRAM_ARCH/uniSkim -d {} -j {} -r {} -y {} -l {} -i {} -o {} \n'.format(
            args.samplename, args.job, args.recoil, args.year, args.lepton, input_file, output_file)
        bashScript += 'xrdcp {} root://cmseos.fnal.gov//store/user/{}/{}/{}'.format(
            output_file, pwd.getpwuid(os.getuid())[0], jobName, args.samplename)
        with open(bash_name, 'w') as file:
            file.write(bashScript)
        os.system('chmod +x {}'.format(bash_name))
        i += 1

    print 'All executables have been written.'
    if not args.dryrun:
        print 'Now submitting to condor...'
        os.system('condor_submit {}'.format(config_name))

    return


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the desired analyzer on FSA n-tuples")
    parser.add_argument('-dr', '--dryrun', action='store_true',
                        help='Create jobs but dont submit')
    parser.add_argument('-j', '--job', action='store', help='job type')
    parser.add_argument('-l', '--lepton', action='store', help='which lepton')
    parser.add_argument('-y', '--year', action='store', help='which year')
    parser.add_argument('-r', '--recoil', action='store', help='recoil type')
    parser.add_argument('-jn', '--jobName', nargs='?', type=str,
                        const='', help='Job Name for condor submission')
    parser.add_argument('-sn', '--samplename', nargs='?',
                        type=str, const='', help='Name of samples')
    parser.add_argument('-sd', '--sampledir', nargs='?',
                        type=str, const='', help='The Sample Input directory')
    args = parser.parse_args()
    main(args)
