#ifndef ROOT_SRC_BASE_TREE_H_
#define ROOT_SRC_BASE_TREE_H_

#include <vector>

#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "TLorentzVector.h"
#include "TTree.h"

enum energy_scale { tes, efake, mfake };

class base_tree {
 protected:
  Int_t recoil;
  TTree *tree, *original;
  bool isMC, isEmbed;
  std::vector<Int_t> good_events;
  TLorentzVector ele, mu, tau, MET, MET_UESUp, MET_UESDown, MET_JESUp, MET_JESDown;

 public:
  // Member functions
  base_tree() {}
  base_tree(TTree* orig, TTree* itree, bool isMC, bool isEmbed, bool isSignal, Int_t rec);
  virtual ~base_tree() {}
  virtual void do_skimming(TH1F*) = 0;
  virtual void set_branches() = 0;
  virtual TTree* fill_tree(RecoilCorrector, MEtSys) = 0;
};

#endif  // ROOT_SRC_BASE_TREE_H_
