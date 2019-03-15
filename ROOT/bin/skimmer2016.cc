// Copyright 2018 Tyler Mitchell

// general includes
#include <dirent.h>
#include <sys/types.h>
#include <fstream>
#include <string>
#include <vector>

// ROOT includes
#include "TChain.h"
#include "TFile.h"
#include "TH1F.h"
#include "TNamed.h"
#include "TTree.h"

// user includes
#include "ltau_skimmer/ROOT/src/CLParser.h"
#include "ltau_skimmer/ROOT/src/base_tree.h"
#include "ltau_skimmer/ROOT/src/etau_tree_2016.h"
#include "ltau_skimmer/ROOT/src/mutau_tree_2016.h"

static unsigned events(0);
int main(int argc, char *argv[]) {
    CLParser parser(argc, argv);
    std::string dir_name = parser.Option("-d");
    std::string job_type = parser.Option("-j");
    std::string inrecoil = parser.Option("-r");
    std::string lepton = parser.Option("-l");
    std::string ifile = parser.Option("-i");
    std::string ofile = parser.Option("-o");

    // recoil corrections
    int recoil(0);
    if (inrecoil.find("W") != std::string::npos) {
        recoil = 1;
    } else if (inrecoil.find("Z") != std::string::npos) {
        recoil = 2;
    }

    // set flag for MC or data
    bool isMC(true);
    if (job_type == "data") {
        isMC = false;
    }

    bool isEmbed(false);
    if (job_type == "embed") {
        isEmbed = true;
    }

    std::string treename;
    if (lepton == "et") {
        treename = "etau_tree";
    } else if (lepton == "mt") {
        truename = "mutau_tree";
    }

    RecoilCorrector recoilPFMetCorrector("HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root");
    TH1F *nevents = new TH1F("nevents", "N(events)", 2, 0.5, 2.5);
    TH1F *cutflow = new TH1F("cutflow", "cutflow", 10, 0.5, 10.5);

    auto open_file = new TFile(ifile.c_str(), "READ");
    auto ntuple = reinterpret_cast<TTree *>(open_file->Get((lepton+"/final/Ntuple").c_str());
    auto evt_count = reinterpret_cast<TH1F *>(open_file->Get((lepton+"/eventCount").c_str())->Clone());
    auto wt_count = reinterpret_cast<TH1F *>(open_file->Get((lepton+"/summedWeights").c_str())->Clone());

    nevents->SetBinContent(1, evt_count->Integral());
    nevents->SetBinContent(2, wt_count->Integral());

    auto fout = new TFile(ofile.c_str(), "RECREATE");

    TTree *newtree = new TTree(treename.c_str(), truename.c_str());
    base_tree *skimmer;
    if (lepton == "et") {
        skimmer = new etau_tree(ntuple, newtree, isMC, isEmbed, recoil);
    }
    skimmer->do_skimming(cutflow);
    auto skimmed_tree = skimmer->fill_tree(recoilPFMetCorrector);
    events += skimmed_tree->GetEntries();

    open_file->Close();
    fout->cd();
    nevents->Write();
    cutflow->Write();
    skimmed_tree->Write();
    fout->Close();

    std::cout << "\n"
                << events << " events saved in the tree.\n"
                << std::endl;
    return 0;
}