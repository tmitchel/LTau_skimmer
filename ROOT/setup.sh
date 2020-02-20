#!/bin/bash

pushd ${CMSSW_BASE}/src/LTau_skimmer/ROOT/data
wget https://github.com/cms-tau-pog/TauIDSFs/raw/master/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2016Legacy.root
wget https://github.com/cms-tau-pog/TauIDSFs/raw/master/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2017ReReco.root
wget https://github.com/cms-tau-pog/TauIDSFs/raw/master/data/TauFES_eta-dm_DeepTau2017v2p1VSe_2018ReReco.root
popd