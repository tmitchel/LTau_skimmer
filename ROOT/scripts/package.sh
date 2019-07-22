#!/bin/bash
pushd $CMSSW_BASE/src
tar --exclude="ltau_skimmer/.git" --exclude="HTT-utilities/.git" czf ltau_skimmer_packge.tar.gz ltau_skimmer HTT-utilities
xrdcp ltau_skimmer_packge.tar.gz root://cmseos.fnal.gov//store/user/tmitchel/ltau_skimmer_packge.tar.gz
rm ltau_skimmer_packge.tar.gz
popd