
#include <vector>

#include "TTree.h"

class base_tree {
 protected:
    TTree *tree, *original;
    bool isMC, isEmbed;
    std::vector<Int_t> good_events;
    TLorentzVector ele, mu, tau, MET, MET_UESUp, MET_UESDown, MET_JESUp, MET_JESDown;

 public:
    // Member functions
    base_tree(TTree* orig, TTree* itree, bool isMC, bool isEmbed, Int_t rec);
    virtual ~base_tree() {}
    void do_skimming(TH1F*);
    void set_branches();
    TTree* fill_tree(RecoilCorrector recoilPFMetCorrector);
}