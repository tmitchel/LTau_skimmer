// Copyright 2018 Tyler Mitchell

#ifndef ROOT_SRC_MUTAU_TREE_2017_H_
#define ROOT_SRC_MUTAU_TREE_2017_H_

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "../interface/mt_2017_input_branches.h"
#include "./base_tree.h"
#include "RecoilCorrector.h"
#include "TLorentzVector.h"
#include "TTree.h"

class mutau_tree2017 : public virtual base_tree {
 private:
  TTree *tree, *original;
  mt_2017_input_branches* in;
  bool isMC, isEmbed;
  std::vector<Int_t> good_events;
  TLorentzVector mu, tau, MET, MET_UESUp, MET_UESDown, MET_JESUp, MET_JESDown;

 public:
  // Member variables
  UInt_t Run, Lumi;
  Int_t recoil;
  Float_t placeholder;  // for all branches not present in 2018

  // // Constructed while running
  Int_t gen_match_1, gen_match_2, njets, nbtag, njetspt20;
  Float_t jetVeto20, jetVeto30, met, metphi, met_px, met_py, extraelec_veto, extramuon_veto, dilepton_veto, pfmetcorr_ex, pfmetcorr_ey;
  Float_t pfmetcorr_ex_UESUp, pfmetcorr_ey_UESUp, pfmetcorr_ex_UESDown, pfmetcorr_ey_UESDown, pfmetcorr_ex_JESUp, pfmetcorr_ey_JESUp,
      pfmetcorr_ex_JESDown, pfmetcorr_ey_JESDown;
  Float_t met_UESUp, met_UESDown, met_JESUp, met_JESDown, metphi_UESUp, metphi_UESDown, metphi_JESUp, metphi_JESDown;
  Float_t pt_1, eta_1, phi_1, m_1, e_1, px_1, py_1, pz_1, pt_2, eta_2, phi_2, m_2, e_2, px_2, py_2, pz_2;

  // Member functions
  mutau_tree2017(TTree* orig, TTree* itree, bool isMC, bool isEmbed, Int_t rec);
  virtual ~mutau_tree2017() {}
  void do_skimming(TH1F*);
  void set_branches();
  Float_t do_tes_met_corr(Float_t, Float_t, Float_t, Float_t, TLorentzVector&, TLorentzVector);
  TTree* fill_tree(RecoilCorrector recoilPFMetCorrector);
};

//////////////////////////////////////////////////////////////////
// Purpose: Initialize tree and original, then read branches    //
//          needed for skimming/sorting from Original           //
//////////////////////////////////////////////////////////////////
// Parameters                                                   //
//   - Original: A tree read from an input root file.           //
//               At this point, only the branches necessary     //
//               for skim selection and sorting are read        //
//   - itree: A newly constructed tree. This tree will be       //
//            filled for all events passing the skim selection  //
//////////////////////////////////////////////////////////////////
mutau_tree2017::mutau_tree2017(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec) : tree(itree),
                                                                                                    original(Original),
                                                                                                    in(new mt_2017_input_branches(Original)),
                                                                                                    isMC(IsMC),
                                                                                                    isEmbed(IsEmbed),
                                                                                                    recoil(rec) {
  original->SetBranchStatus("double*", 0);
  original->SetBranchStatus("Double*", 0);
  original->SetBranchStatus("Ele*", 0);
}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void mutau_tree2017::do_skimming(TH1F* cutflow) {
  // declare variables for sorting
  ULong64_t evt_now(0);
  ULong64_t evt_before(1);
  int best_evt(-1);
  std::pair<float, float> muCandidate, tauCandidate;

  std::cout << "Starting the skim..." << std::endl;

  Int_t nevt = (Int_t)original->GetEntries();
  for (auto ievt = 0; ievt < nevt; ievt++) {
    original->GetEntry(ievt);
    evt_now = in->evt;

    // TLorentzVector ele, tau;
    mu.SetPtEtaPhiM(in->mPt, in->mEta, in->mPhi, in->mMass);
    tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);

    // apply TES
    if (isMC && !isEmbed) {
      if (in->tZTTGenMatching == 5) {
        if (in->tDecayMode == 0) {
          tau *= 1.007;
        } else if (in->tDecayMode == 1) {
          tau *= 0.998;
        } else if (in->tDecayMode == 10) {
          tau *= 1.001;
        }
      } else if (in->tZTTGenMatching == 1 || in->tZTTGenMatching == 3) {
        if (in->tDecayMode == 0) {
          tau *= 1.003;
        } else if (in->tDecayMode == 1) {
          tau *= 1.036;
        }
      }
    } else if (isEmbed) {
      if (in->tZTTGenMatching == 5) {
        if (in->tDecayMode == 0) {
          tau *= 0.975;
        } else if (in->tDecayMode == 1) {
          tau *= (0.975 * 1.051);
        } else if (in->tDecayMode == 10) {
          tau *= (0.975 * 0.975 * 0.975);
        }
      }
    }

    float mu_pt_min(21.), tau_pt_min(23.);

    cutflow->Fill(1., 1.);
    // apply event selection

    auto Mu24 = in->IsoMu24Pass && in->mMatchesIsoMu24Path && in->mMatchesIsoMu24Filter;
    auto Mu27 = in->IsoMu27Pass && in->mMatchesIsoMu27Path && in->mMatchesIsoMu27Filter;
    auto Cross = in->Mu20Tau27Pass && in->mMatchesIsoMu20Tau27Filter && in->mMatchesIsoMu20Tau27Path && in->tMatchesIsoMu20Tau27Filter && in->tMatchesIsoMu20Tau27Path;

    if (Mu24 || Mu27 || Cross)
      cutflow->Fill(2., 1.);
    else
      continue;

    if (in->mPt > mu_pt_min && fabs(in->mEta) < 2.1 && fabs(in->mPVDZ) < 0.2 && fabs(in->mPVDXY) < 0.045)
      cutflow->Fill(3., 1.);  // electron kinematic selection
    else
      continue;

    bool goodglob = in->mIsGlobal && in->mNormalizedChi2 < 3 && in->mChi2LocalPosition < 12 && in->mTrkKink < 20;
    bool isMedium = in->mPFIDLoose && in->mValidFraction > 0.49 && in->mSegmentCompatibility > (goodglob ? 0.303 : 0.451);

    if (in->mPFIDMedium)
      cutflow->Fill(4., 1.);  // muon quality selection
    else
      continue;

    if (isMC || isEmbed || isMedium)
      cutflow->Fill(5., 1.);  // muon quality selection
    else
      continue;

    if (tau.Pt() > tau_pt_min && fabs(tau.Eta()) < 2.3 && fabs(in->tPVDZ) < 0.2)
      cutflow->Fill(6., 1.);  // tau kinematic selection
    else
      continue;

    if (in->tRerunMVArun2v2DBoldDMwLTVLoose && in->tDecayModeFinding > 0 && fabs(in->tCharge) < 2)
      cutflow->Fill(7., 1.);  // tau quality selection
    else
      continue;

    if (in->tAgainstMuonTight3 > 0.5 && in->tAgainstElectronVLooseMVA6 > 0.5)
      cutflow->Fill(8., 1.);  // tau against leptons
    else
      continue;

    if (in->muVetoZTTp001dxyzR0 < 2 && in->eVetoZTTp001dxyzR0 == 0 && in->dimuonVeto == 0)
      cutflow->Fill(9., 1.);  // vetos
    else
      continue;

    // implement new sorting per
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
    if (evt_now != evt_before) {  // new event, save the tau candidates
      // since it is new event, do we have the best entry to save? If yes, save it!
      if (best_evt > -1)
        good_events.push_back(best_evt);

      //  this is a new event, so the first tau pair is the best! :)
      best_evt = ievt;
      muCandidate = std::make_pair(in->mPt, in->mRelPFIsoDBDefaultR04);
      tauCandidate = std::make_pair(in->tPt, in->tRerunMVArun2v2DBoldDMwLTraw);
    } else {  // not a new event
      std::pair<float, float> currEleCandidate(in->mPt, in->mRelPFIsoDBDefaultR04);
      std::pair<float, float> currTauCandidate(in->tPt, in->tRerunMVArun2v2DBoldDMwLTraw);

      // clause 1, select the pair that has most isolated tau lepton 1
      if (currEleCandidate.second - muCandidate.second > 0.0001) best_evt = ievt;

      // check if the first tau is the same, and if so - move to clause 2
      if (fabs(currEleCandidate.second - muCandidate.second) < 0.0001) {
        // pick up  the pair with the highest pT of the first candidate
        if (currEleCandidate.first - muCandidate.first > 0.0001) best_evt = ievt;
        if (fabs(currEleCandidate.first - muCandidate.first) < 0.0001) {
          // same pT, same iso, move to clause 3
          if (currTauCandidate.second - tauCandidate.second > 0.0001) best_evt = ievt;
          if (fabs(currTauCandidate.second - tauCandidate.second) < 0.0001) {
            // same iso - pick the pair with the highest pT
            if (currTauCandidate.first - tauCandidate.first > 0.0001) best_evt = ievt;
          }  // tau2 has the same isolation
        }    // tau1 has the same pT
      }      // tau1 has the same isolation
    }        // not a new event
    evt_before = evt_now;
  }
  if (best_evt > -1)
    good_events.push_back(best_evt);
}

Float_t mutau_tree2017::do_tes_met_corr(Float_t decayMode, Float_t sf1, Float_t sf2, Float_t sf3, TLorentzVector& met, TLorentzVector tau) {
  if (decayMode == 0) {
    met = met + tau - sf1 * tau;
    return sf1;
  } else if (decayMode == 1) {
    met = met + tau - sf2 * tau;
    return sf2;
  } else if (decayMode == 10) {
    met = met + tau - sf3 * tau;
    return sf3;
  }
  return 1.;
}

//////////////////////////////////////////////////////////////////
// Purpose: Fill tree with variables from original and new      //
//          variables. Only events that have been skimmed,      //
//          sorted, and stored in the good_events vector will   //
//          be stored in the tree.                              //
//////////////////////////////////////////////////////////////////
// Return: The same TTree passed to the constructor and stored  //
//         in original, but now it is filled with good events   //
//////////////////////////////////////////////////////////////////
TTree* mutau_tree2017::fill_tree(RecoilCorrector recoilPFMetCorrector) {
  std::cout << "setting branches..." << std::endl;
  set_branches();  // get all the branches set up
  std::cout << "branches set." << std::endl;

  // loop through all events pasing skimming/sorting
  for (auto& ievt : good_events) {
    original->GetEntry(ievt);

    Run = in->run;
    Lumi = in->lumi;

    // convert from Float_t in FSA to Int_t for analyzer
    gen_match_1 = in->mZTTGenMatching;
    gen_match_2 = in->tZTTGenMatching;
    njets = in->jetVeto30WoNoisyJets;
    nbtag = in->bjetDeepCSVVeto20MediumWoNoisyJets;
    njetspt20 = in->jetVeto20WoNoisyJets;

    // TLorentzVector ele, tau;
    mu.SetPtEtaPhiM(in->mPt, in->mEta, in->mPhi, in->mMass);
    tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);

    met_px = in->type1_pfMetEt * cos(in->type1_pfMetPhi);
    met_py = in->type1_pfMetEt * sin(in->type1_pfMetPhi);

    extraelec_veto = in->eVetoZTTp001dxyzR0 > 1;
    extramuon_veto = in->muVetoZTTp001dxyzR0 > 0;
    dilepton_veto = in->dimuonVeto > 0;

    // TLorentzVector ele, tau;
    mu.SetPtEtaPhiM(in->mPt, in->mEta, in->mPhi, in->mMass);
    tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);
    MET.SetPtEtaPhiM(in->type1_pfMetEt, 0, in->type1_pfMetPhi, 0);
    MET_UESUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_UnclusteredEnUp, 0, in->type1_pfMet_shiftedPhi_UnclusteredEnUp, 0);
    MET_UESDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_UnclusteredEnDown, 0, in->type1_pfMet_shiftedPhi_UnclusteredEnDown, 0);
    MET_JESUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEnUp, 0, in->type1_pfMet_shiftedPhi_JetEnUp, 0);
    MET_JESDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEnDown, 0, in->type1_pfMet_shiftedPhi_JetEnDown, 0);

    pfmetcorr_ex = MET.Px();
    pfmetcorr_ey = MET.Py();
    pfmetcorr_ex_UESUp = MET_UESUp.Px();
    pfmetcorr_ey_UESUp = MET_UESUp.Py();
    pfmetcorr_ex_UESDown = MET_UESDown.Px();
    pfmetcorr_ey_UESDown = MET_UESDown.Py();
    pfmetcorr_ex_JESUp = MET_JESUp.Px();
    pfmetcorr_ey_JESUp = MET_JESUp.Py();
    pfmetcorr_ex_JESDown = MET_JESDown.Px();
    pfmetcorr_ey_JESDown = MET_JESDown.Py();

    if (recoil == 1) {
      recoilPFMetCorrector.CorrectByMeanResolution(
          MET.Px(),       // uncorrected type I pf met px (float)
          MET.Py(),       // uncorrected type I pf met py (float)
          in->genpX,      // generator Z/W/Higgs px (float)
          in->genpY,      // generator Z/W/Higgs py (float)
          in->vispX,      // generator visible Z/W/Higgs px (float)
          in->vispY,      // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets + 1,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets + 1,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESDown);  // corrected type I pf met py (float)

    } else if (recoil == 2) {
      recoilPFMetCorrector.CorrectByMeanResolution(
          MET.Px(),       // uncorrected type I pf met px (float)
          MET.Py(),       // uncorrected type I pf met py (float)
          in->genpX,      // generator Z/W/Higgs px (float)
          in->genpY,      // generator Z/W/Higgs py (float)
          in->vispX,      // generator visible Z/W/Higgs px (float)
          in->vispY,      // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets,      // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets,              // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30WoNoisyJets,              // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESDown);  // corrected type I pf met py (float)
    }

    MET.SetPxPyPzE(pfmetcorr_ex, pfmetcorr_ey, 0, sqrt(pfmetcorr_ex * pfmetcorr_ex + pfmetcorr_ey * pfmetcorr_ey));
    MET_UESUp.SetPxPyPzE(pfmetcorr_ex_UESUp, pfmetcorr_ey_UESUp, 0, sqrt(pfmetcorr_ex_UESUp * pfmetcorr_ex_UESUp + pfmetcorr_ey_UESUp * pfmetcorr_ey_UESUp));
    MET_UESDown.SetPxPyPzE(pfmetcorr_ex_UESDown, pfmetcorr_ey_UESDown, 0, sqrt(pfmetcorr_ex_UESDown * pfmetcorr_ex_UESDown + pfmetcorr_ey_UESDown * pfmetcorr_ey_UESDown));
    MET_JESUp.SetPxPyPzE(pfmetcorr_ex_JESUp, pfmetcorr_ey_JESUp, 0, sqrt(pfmetcorr_ex_JESUp * pfmetcorr_ex_JESUp + pfmetcorr_ey_JESUp * pfmetcorr_ey_JESUp));
    MET_JESDown.SetPxPyPzE(pfmetcorr_ex_JESDown, pfmetcorr_ey_JESDown, 0, sqrt(pfmetcorr_ex_JESDown * pfmetcorr_ex_JESDown + pfmetcorr_ey_JESDown * pfmetcorr_ey_JESDown));

    if (isMC && !isEmbed) {
      // met correction due to tau energy scale
      if (in->tZTTGenMatching == 5) {
        auto sf = do_tes_met_corr(in->tDecayMode, 1.007, 0.998, 1.001, MET, tau);
        do_tes_met_corr(in->tDecayMode, 1.007, 0.998, 1.001, MET_JESUp, tau);
        do_tes_met_corr(in->tDecayMode, 1.007, 0.998, 1.001, MET_JESDown, tau);
        do_tes_met_corr(in->tDecayMode, 1.007, 0.998, 1.001, MET_UESUp, tau);
        do_tes_met_corr(in->tDecayMode, 1.007, 0.998, 1.001, MET_UESDown, tau);
        tau *= sf;
      } else if (in->tZTTGenMatching == 1 || in->tZTTGenMatching == 3) {
        auto sf = do_tes_met_corr(in->tDecayMode, 1.003, 1.036, 1.00, MET, tau);
        do_tes_met_corr(in->tDecayMode, 1.003, 1.036, 1.00, MET_JESUp, tau);
        do_tes_met_corr(in->tDecayMode, 1.003, 1.036, 1.00, MET_JESDown, tau);
        do_tes_met_corr(in->tDecayMode, 1.003, 1.036, 1.00, MET_UESUp, tau);
        do_tes_met_corr(in->tDecayMode, 1.003, 1.036, 1.00, MET_UESDown, tau);
        tau *= sf;
      }
    } else if (isEmbed) {
      if (in->tZTTGenMatching == 5) {
        auto sf = do_tes_met_corr(in->tDecayMode, 0.975, 0.975 * 1.051, 0.975 * 0.975 * 0.975, MET, tau);
        do_tes_met_corr(in->tDecayMode, 0.975, 0.975 * 1.051, 0.975 * 0.975 * 0.975, MET_JESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.975, 0.975 * 1.051, 0.975 * 0.975 * 0.975, MET_JESDown, tau);
        do_tes_met_corr(in->tDecayMode, 0.975, 0.975 * 1.051, 0.975 * 0.975 * 0.975, MET_UESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.975, 0.975 * 1.051, 0.975 * 0.975 * 0.975, MET_UESDown, tau);
        tau *= sf;
      }
    }

    met = MET.Pt();
    metphi = MET.Phi();
    met_px = MET.Px();
    met_py = MET.Py();

    met_JESUp = MET_JESUp.Pt();
    met_JESDown = MET_JESDown.Pt();
    met_UESUp = MET_UESUp.Pt();
    met_UESDown = MET_UESDown.Pt();
    metphi_JESUp = MET_JESUp.Phi();
    metphi_JESDown = MET_JESDown.Phi();
    metphi_UESUp = MET_UESUp.Phi();
    metphi_UESDown = MET_UESDown.Phi();

    m_1 = mu.M();
    px_1 = mu.Px();
    py_1 = mu.Py();
    pz_1 = mu.Pz();
    e_1 = mu.E();
    pt_1 = mu.Pt();
    phi_1 = mu.Phi();
    eta_1 = mu.Eta();
    m_2 = tau.M();
    px_2 = tau.Px();
    py_2 = tau.Py();
    pz_2 = tau.Pz();
    e_2 = tau.E();
    pt_2 = tau.Pt();
    phi_2 = tau.Phi();
    eta_2 = tau.Eta();

    tree->Fill();
  }
  return tree;
}

//////////////////////////////////////////////////////////////////
// Purpose:                                                     //
//   - Read all branches we need from original. These branches  //
//          branches may be directly stored, be copied with a   //
//          new name, or used to construct new variables.       //
//   - Create branches in tree to store any variable we want    //
//////////////////////////////////////////////////////////////////
void mutau_tree2017::set_branches() {
  // output file branches
  tree->Branch("evt", &in->evt);
  tree->Branch("run", &Run);
  tree->Branch("lumi", &Lumi);
  tree->Branch("gen_match_1", &gen_match_1, "gen_match_1/I");
  tree->Branch("gen_match_2", &gen_match_2, "gen_match_2/I");
  tree->Branch("njets", &njets, "njets/I");
  tree->Branch("nbtag", &nbtag, "nbtag/I");
  tree->Branch("njetspt20", &njetspt20, "njetspt20/I");
  tree->Branch("vbfMass", &in->vbfMassWoNoisyJets, "vbfMass/F");
  tree->Branch("vbfMassWoNoisyJets", &in->vbfMassWoNoisyJets, "vbfMassWoNoisyJets/F");

  tree->Branch("mMatchesIsoMu20Tau27Path", &in->mMatchesIsoMu20Tau27Path, "mMatchesIsoMu20Tau27Path/F");
  tree->Branch("mMatchesIsoMu24Filter", &in->mMatchesIsoMu24Filter, "mMatchesIsoMu24Filter/F");
  tree->Branch("mMatchesIsoMu24Path", &in->mMatchesIsoMu24Path, "mMatchesIsoMu24Path/F");
  tree->Branch("mMatchesIsoMu27Filter", &in->mMatchesIsoMu27Filter, "mMatchesIsoMu27Filter/F");
  tree->Branch("mMatchesIsoMu27Path", &in->mMatchesIsoMu27Path, "mMatchesIsoMu27Path/F");
  tree->Branch("Mu20Tau27Pass", &in->Mu20Tau27Pass, "Mu20Tau27Pass/F");
  tree->Branch("IsoMu27Pass", &in->IsoMu27Pass, "IsoMu27Pass/F");
  tree->Branch("IsoMu24Pass", &in->IsoMu24Pass, "IsoMu24Pass/F");
  tree->Branch("tMatchesIsoMu20Tau27Filter", &in->tMatchesIsoMu20Tau27Filter, "tMatchesIsoMu20Tau27Filter/F");
  tree->Branch("tMatchesIsoMu20Tau27Path", &in->tMatchesIsoMu20Tau27Path, "tMatchesIsoMu20Tau27Path/F");

  tree->Branch("met_px", &met_px, "met_px/F");
  tree->Branch("met_py", &met_py, "met_py/F");
  tree->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/F");
  tree->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/F");
  tree->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/F");
  tree->Branch("pfmetcorr_ex", &pfmetcorr_ex, "pfmetcorr_ex/F");
  tree->Branch("pfmetcorr_ey", &pfmetcorr_ey, "pfmetcorr_ey/F");
  tree->Branch("genpX", &in->genpX, "genpX/F");
  tree->Branch("genpY", &in->genpY, "genpY/F");
  tree->Branch("vispX", &in->vispX, "vispX/F");
  tree->Branch("vispY", &in->vispY, "vispY/F");
  tree->Branch("met", &met, "met/F");
  tree->Branch("metphi", &metphi, "metphi/F");
  tree->Branch("met_px", &met_px, "met_px/F");
  tree->Branch("met_py", &met_py, "met_py/F");
  tree->Branch("pt_1", &pt_1, "pt_1/F");
  tree->Branch("eta_1", &eta_1, "eta_1/F");
  tree->Branch("phi_1", &phi_1, "phi_1/F");
  tree->Branch("m_1", &m_1, "m_1/F");
  tree->Branch("e_1", &e_1, "e_1/F");
  tree->Branch("px_1", &px_1, "px_1/F");
  tree->Branch("py_1", &py_1, "py_1/F");
  tree->Branch("pz_1", &pz_1, "pz_1/F");
  tree->Branch("dZ_1", &in->mPVDZ, "dZ_1/F");
  tree->Branch("d0_1", &in->mPVDXY, "d0_1/F");
  tree->Branch("q_1", &in->mCharge, "q_1/F");
  tree->Branch("iso_1", &in->mRelPFIsoDBDefaultR04, "iso_1/F");
  tree->Branch("pt_2", &pt_2, "pt_2/F");
  tree->Branch("eta_2", &eta_2, "eta_2/F");
  tree->Branch("phi_2", &phi_2, "phi_2/F");
  tree->Branch("m_2", &m_2, "m_2/F");
  tree->Branch("e_2", &e_2, "e_2/F");
  tree->Branch("px_2", &px_2, "px_2/F");
  tree->Branch("py_2", &py_2, "py_2/F");
  tree->Branch("pz_2", &pz_2, "pz_2/F");
  tree->Branch("dZ_2", &in->tPVDZ, "dZ_2/F");
  tree->Branch("q_2", &in->tCharge, "q_2/F");
  tree->Branch("iso_2", &in->tRerunMVArun2v2DBoldDMwLTraw, "iso_2/F");
  tree->Branch("decayModeFinding_2", &in->tDecayModeFinding, "decayModeFinding_2/F");
  tree->Branch("decayModeFindingNewDMs_2", &in->tDecayModeFindingNewDMs, "decayModeFindingNewDMs_2/F");
  tree->Branch("l2_decayMode", &in->tDecayMode, "l2_decayMode/F");
  tree->Branch("id_m_medium_1", &in->mPFIDMedium, "id_m_medium_1/F");

  tree->Branch("byVLooseIsolationMVArun2v1DBoldDMwLT_2", &in->tByVLooseIsolationMVArun2v1DBoldDMwLT, "byVLooseIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_2", &in->tByLooseIsolationMVArun2v1DBoldDMwLT, "byLooseIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_2", &in->tByMediumIsolationMVArun2v1DBoldDMwLT, "byMediumIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2", &in->tByTightIsolationMVArun2v1DBoldDMwLT, "byTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byVTightIsolationMVArun2v1DBoldDMwLT_2", &in->tByVTightIsolationMVArun2v1DBoldDMwLT, "byVTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byVVTightIsolationMVArun2v1DBoldDMwLT_2", &in->tByVVTightIsolationMVArun2v1DBoldDMwLT, "byVVTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVLoose", &in->tRerunMVArun2v2DBoldDMwLTVLoose, "tRerunMVArun2v2DBoldDMwLTVLoose/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTLoose", &in->tRerunMVArun2v2DBoldDMwLTLoose, "tRerunMVArun2v2DBoldDMwLTLoose/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTMedium", &in->tRerunMVArun2v2DBoldDMwLTMedium, "tRerunMVArun2v2DBoldDMwLTMedium/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTTight", &in->tRerunMVArun2v2DBoldDMwLTTight, "tRerunMVArun2v2DBoldDMwLTTight/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVTight", &in->tRerunMVArun2v2DBoldDMwLTVTight, "tRerunMVArun2v2DBoldDMwLTVTight/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVVTight", &in->tRerunMVArun2v2DBoldDMwLTVVTight, "tRerunMVArun2v2DBoldDMwLTVVTight/F");

  tree->Branch("againstElectronTightMVA6_2", &in->tAgainstElectronTightMVA6, "againstElectronTightMVA6_2/F");
  tree->Branch("againstElectronVLooseMVA6_2", &in->tAgainstElectronVLooseMVA6, "againstElectronVLooseMVA6_2/F");
  tree->Branch("againstMuonTight3_2", &in->tAgainstMuonLoose3, "againstMuonTight3_2/F");
  tree->Branch("againstMuonLoose3_2", &in->tAgainstMuonTight3, "againstMuonLoose3_2/F");

  tree->Branch("rho", &in->rho, "rho/F");
  tree->Branch("metcov00", &in->metcov00, "metcov00/F");
  tree->Branch("metcov01", &in->metcov01, "metcov01/F");
  tree->Branch("metcov10", &in->metcov10, "metcov10/F");
  tree->Branch("metcov11", &in->metcov11, "metcov11/F");
  tree->Branch("metcov00_v2", &in->metcov00_DESYlike, "metcov00_v2/F");
  tree->Branch("metcov01_v2", &in->metcov01_DESYlike, "metcov01_v2/F");
  tree->Branch("metcov10_v2", &in->metcov10_DESYlike, "metcov10_v2/F");
  tree->Branch("metcov11_v2", &in->metcov11_DESYlike, "metcov11_v2/F");
  tree->Branch("NUP", &in->NUP, "NUP/F");
  tree->Branch("genM", &in->genM, "genM/F");
  tree->Branch("genpT", &in->genpT, "genpT/F");
  tree->Branch("genEta", &in->genEta, "genEta/F");
  tree->Branch("numGenJets", &in->numGenJets, "numGenJets/F");
  tree->Branch("npu", &in->nTruePU, "npu/F");
  tree->Branch("npv", &in->nvtx, "npv/F");
  tree->Branch("genweight", &in->GenWeight, "genweight/F");
  tree->Branch("metSig", &in->metSig, "metSig/F");
  tree->Branch("Rivet_higgsPt", &in->Rivet_higgsPt, "Rivet_higgsPt/F");
  tree->Branch("Rivet_nJets30", &in->Rivet_nJets30, "Rivet_nJets30/F");

  tree->Branch("jpt_1", &in->j1ptWoNoisyJets, "jpt_1/F");
  tree->Branch("jeta_1", &in->j1etaWoNoisyJets, "jeta_1/F");
  tree->Branch("jphi_1", &in->j1phiWoNoisyJets, "jphi_1/F");
  tree->Branch("jcsv_1", &in->j1csvWoNoisyJets, "jcsv_1/F");
  tree->Branch("jpt_2", &in->j2ptWoNoisyJets, "jpt_2/F");
  tree->Branch("jeta_2", &in->j2etaWoNoisyJets, "jeta_2/F");
  tree->Branch("jphi_2", &in->j2phiWoNoisyJets, "jphi_2/F");
  tree->Branch("jcsv_2", &in->j2csvWoNoisyJets, "jcsv_2/F");
  tree->Branch("bpt_1", &in->jb1ptWoNoisyJets, "bpt_1/F");
  tree->Branch("beta_1", &in->jb1etaWoNoisyJets, "beta_1/F");
  tree->Branch("bphi_1", &in->jb1phiWoNoisyJets, "bphi_1/F");
  tree->Branch("bcsv_1", &in->jb1csvWoNoisyJets, "bcsv_1/F");
  tree->Branch("bflavor_1", &in->jb1hadronflavorWoNoisyJets, "bflavor_1/F");
  tree->Branch("bpt_2", &in->jb2ptWoNoisyJets, "bpt_2/F");
  tree->Branch("beta_2", &in->jb2etaWoNoisyJets, "beta_2/F");
  tree->Branch("bphi_2", &in->jb2phiWoNoisyJets, "bphi_2/F");
  tree->Branch("bcsv_2", &in->jb2csvWoNoisyJets, "bcsv_2/F");
  tree->Branch("bflavor_2", &in->jb2hadronflavorWoNoisyJets, "bflavor_2/F");

  tree->Branch("bjetDeepCSVVeto20Tight", &in->bjetDeepCSVVeto20Tight, "bjetDeepCSVVeto20Tight/F");
  tree->Branch("bjetDeepCSVVeto30Loose", &in->bjetDeepCSVVeto30Loose, "bjetDeepCSVVeto30Loose/F");
  tree->Branch("bjetDeepCSVVeto30Medium", &in->bjetDeepCSVVeto30Medium, "bjetDeepCSVVeto30Medium/F");
  tree->Branch("bjetDeepCSVVeto30Tight", &in->bjetDeepCSVVeto30Tight, "bjetDeepCSVVeto30Tight/F");

  tree->Branch("topQuarkPt1", &in->topQuarkPt1, "topQuarkPt1/F");
  tree->Branch("topQuarkPt2", &in->topQuarkPt2, "topQuarkPt2/F");

  tree->Branch("tZTTGenPt", &in->tZTTGenPt, "tZTTGenPt/F");
  tree->Branch("tZTTGenPhi", &in->tZTTGenPhi, "tZTTGenPhi/F");
  tree->Branch("tZTTGenEta", &in->tZTTGenEta, "tZTTGenEta/F");
  tree->Branch("tZTTGenDR", &in->tZTTGenDR, "tZTTGenDR/F");
  tree->Branch("tGenDecayMode", &in->tGenDecayMode, "tGenDecayMode/F");
  tree->Branch("tGenEnergy", &in->tGenEnergy, "tGenEnergy/F");
  tree->Branch("tGenEta", &in->tGenEta, "tGenEta/F");
  tree->Branch("tGenJetEta", &in->tGenJetEta, "tGenJetEta/F");
  tree->Branch("tGenJetPt", &in->tGenJetPt, "tGenJetPt/F");
  tree->Branch("tGenMotherEnergy", &in->tGenMotherEnergy, "tGenMotherEnergy/F");
  tree->Branch("tGenMotherEta", &in->tGenMotherEta, "tGenMotherEta/F");
  tree->Branch("tGenMotherPdgId", &in->tGenMotherPdgId, "tGenMotherPdgId/F");
  tree->Branch("tGenMotherPhi", &in->tGenMotherPhi, "tGenMotherPhi/F");
  tree->Branch("tGenMotherPt", &in->tGenMotherPt, "tGenMotherPt/F");
  tree->Branch("tGenPdgId", &in->tGenPdgId, "tGenPdgId/F");
  tree->Branch("tGenPhi", &in->tGenPhi, "tGenPhi/F");
  tree->Branch("tGenPt", &in->tGenPt, "tGenPt/F");
  tree->Branch("tGenStatus", &in->tGenStatus, "tGenStatus/F");
  tree->Branch("mGenCharge", &in->mGenCharge, "mGenCharge/F");
  tree->Branch("mGenDirectPromptTauDecayFinalState", &in->mGenDirectPromptTauDecayFinalState, "mGenDirectPromptTauDecayFinalState/F");
  tree->Branch("mGenEnergy", &in->mGenEnergy, "mGenEnergy/F");
  tree->Branch("mGenEta", &in->mGenEta, "mGenEta/F");
  tree->Branch("mGenIsPrompt", &in->mGenIsPrompt, "mGenIsPrompt/F");
  tree->Branch("mGenMotherPdgId", &in->mGenMotherPdgId, "mGenMotherPdgId/F");
  tree->Branch("mGenParticle", &in->mGenParticle, "mGenParticle/F");
  tree->Branch("mGenPdgId", &in->mGenPdgId, "mGenPdgId/F");
  tree->Branch("mGenPhi", &in->mGenPhi, "mGenPhi/F");
  tree->Branch("mGenPrompt", &in->mGenPrompt, "mGenPrompt/F");
  tree->Branch("mGenPromptTauDecay", &in->mGenPromptTauDecay, "mGenPromptTauDecay/F");
  tree->Branch("mGenPt", &in->mGenPt, "mGenPt/F");
  tree->Branch("mGenTauDecay", &in->mGenTauDecay, "mGenTauDecay/F");
  tree->Branch("mGenVZ", &in->mGenVZ, "mGenVZ/F");
  tree->Branch("mGenVtxPVMatch", &in->mGenVtxPVMatch, "mGenVtxPVMatch/F");

  // Flags
  tree->Branch("Flag_BadChargedCandidateFilter", &in->Flag_BadChargedCandidateFilter, "Flag_BadChargedCandidateFilter/F");
  tree->Branch("Flag_BadPFMuonFilter", &in->Flag_BadPFMuonFilter, "Flag_BadPFMuonFilter/F");
  tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &in->Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/F");
  tree->Branch("Flag_HBHENoiseFilter", &in->Flag_HBHENoiseFilter, "Flag_HBHENoiseFilter/F");
  tree->Branch("Flag_HBHENoiseIsoFilter", &in->Flag_HBHENoiseIsoFilter, "Flag_HBHENoiseIsoFilter/F");
  tree->Branch("Flag_badMuons", &in->Flag_badMuons, "Flag_badMuons/F");
  tree->Branch("Flag_duplicateMuons", &in->Flag_duplicateMuons, "Flag_duplicateMuons/F");
  tree->Branch("Flag_ecalBadCalibFilter", &in->Flag_ecalBadCalibFilter, "Flag_ecalBadCalibFilter/F");
  tree->Branch("Flag_eeBadScFilter", &in->Flag_eeBadScFilter, "Flag_eeBadScFilter/F");
  tree->Branch("Flag_globalSuperTightHalo2016Filter", &in->Flag_globalSuperTightHalo2016Filter, "Flag_globalSuperTightHalo2016Filter/F");
  tree->Branch("Flag_globalTightHalo2016Filter", &in->Flag_globalTightHalo2016Filter, "Flag_globalTightHalo2016Filter/F");
  tree->Branch("Flag_goodVertices", &in->Flag_goodVertices, "Flag_goodVertices/F");

  // Systematics
  tree->Branch("met_JESUp", &met_JESUp, "met_JESUp/F");
  tree->Branch("met_JESDown", &met_JESDown, "met_JESDown/F");
  tree->Branch("met_UESUp", &met_UESUp, "met_UESUp/F");
  tree->Branch("met_UESDown", &met_UESDown, "met_UESDown/F");
  tree->Branch("metphi_JESUp", &metphi_JESUp, "metphi_JESUp/F");
  tree->Branch("metphi_JESDown", &metphi_JESDown, "metphi_JESDown/F");
  tree->Branch("metphi_UESUp", &metphi_UESUp, "metphi_UESUp/F");
  tree->Branch("metphi_UESDown", &metphi_UESDown, "metphi_UESDown/F");

  tree->Branch("pfmetcorr_ex_UESUp", &pfmetcorr_ex_UESUp, "pfmetcorr_ex_UESUp/F");
  tree->Branch("pfmetcorr_ey_UESUp", &pfmetcorr_ey_UESUp, "pfmetcorr_ey_UESUp/F");
  tree->Branch("pfmetcorr_ex_UESDown", &pfmetcorr_ex_UESDown, "pfmetcorr_ex_UESDown/F");
  tree->Branch("pfmetcorr_ey_UESDown", &pfmetcorr_ey_UESDown, "pfmetcorr_ey_UESDown/F");
  tree->Branch("pfmetcorr_ex_JESUp", &pfmetcorr_ex_JESUp, "pfmetcorr_ex_JESUp/F");
  tree->Branch("pfmetcorr_ey_JESUp", &pfmetcorr_ey_JESUp, "pfmetcorr_ey_JESUp/F");
  tree->Branch("pfmetcorr_ex_JESDown", &pfmetcorr_ex_JESDown, "pfmetcorr_ex_JESDown/F");
  tree->Branch("pfmetcorr_ey_JESDown", &pfmetcorr_ey_JESDown, "pfmetcorr_ey_JESDown/F");
  tree->Branch("in->type1_pfMet_shiftedPt_UnclusteredEnUp", &in->type1_pfMet_shiftedPt_UnclusteredEnUp, "in->type1_pfMet_shiftedPt_UnclusteredEnUp/F");
  tree->Branch("in->type1_pfMet_shiftedPhi_UnclusteredEnUp", &in->type1_pfMet_shiftedPhi_UnclusteredEnUp, "in->type1_pfMet_shiftedPhi_UnclusteredEnUp/F");
  tree->Branch("in->type1_pfMet_shiftedPt_UnclusteredEnDown", &in->type1_pfMet_shiftedPt_UnclusteredEnDown, "in->type1_pfMet_shiftedPt_UnclusteredEnDown/F");
  tree->Branch("in->type1_pfMet_shiftedPhi_UnclusteredEnDown", &in->type1_pfMet_shiftedPhi_UnclusteredEnDown, "in->type1_pfMet_shiftedPhi_UnclusteredEnDown/F");
  tree->Branch("in->type1_pfMet_shiftedPt_JetEnUp", &in->type1_pfMet_shiftedPt_JetEnUp, "in->type1_pfMet_shiftedPt_JetEnUp/F");
  tree->Branch("in->type1_pfMet_shiftedPhi_JetEnUp", &in->type1_pfMet_shiftedPhi_JetEnUp, "in->type1_pfMet_shiftedPhi_JetEnUp/F");
  tree->Branch("in->type1_pfMet_shiftedPt_JetEnDown", &in->type1_pfMet_shiftedPt_JetEnDown, "in->type1_pfMet_shiftedPt_JetEnDown/F");
  tree->Branch("in->type1_pfMet_shiftedPhi_JetEnDown", &in->type1_pfMet_shiftedPhi_JetEnDown, "in->type1_pfMet_shiftedPhi_JetEnDown/F");

  tree->Branch("jetVeto30WoNoisyJets_JetEta0to3Down", &in->jetVeto30WoNoisyJets_JetEta0to3Down, "jetVeto30WoNoisyJets_JetEta0to3Down/F");
  tree->Branch("jetVeto30WoNoisyJets_JetEta0to3Up", &in->jetVeto30WoNoisyJets_JetEta0to3Up);
  tree->Branch("jetVeto30WoNoisyJets_JetEta0to5Down", &in->jetVeto30WoNoisyJets_JetEta0to5Down);
  tree->Branch("jetVeto30WoNoisyJets_JetEta0to5Up", &in->jetVeto30WoNoisyJets_JetEta0to5Up);
  tree->Branch("jetVeto30WoNoisyJets_JetEta3to5Down", &in->jetVeto30WoNoisyJets_JetEta3to5Down);
  tree->Branch("jetVeto30WoNoisyJets_JetEta3to5Up", &in->jetVeto30WoNoisyJets_JetEta3to5Up);
  tree->Branch("jetVeto30WoNoisyJets_JetRelativeBalDownWoNoisyJets", &in->jetVeto30WoNoisyJets_JetRelativeBalDownWoNoisyJets);
  tree->Branch("jetVeto30WoNoisyJets_JetRelativeBalUpWoNoisyJets", &in->jetVeto30WoNoisyJets_JetRelativeBalUpWoNoisyJets);
  tree->Branch("jetVeto30WoNoisyJets_JetRelativeSampleDown", &in->jetVeto30WoNoisyJets_JetRelativeSampleDown);
  tree->Branch("jetVeto30WoNoisyJets_JetRelativeSampleUp", &in->jetVeto30WoNoisyJets_JetRelativeSampleUp);
  tree->Branch("jetVeto30WoNoisyJets_JetTotalDown", &in->jetVeto30WoNoisyJets_JetTotalDown);
  tree->Branch("jetVeto30WoNoisyJets_JetTotalUp", &in->jetVeto30WoNoisyJets_JetTotalUp);
  tree->Branch("jetVeto30_JetAbsoluteFlavMapDown", &in->jetVeto30_JetAbsoluteFlavMapDown);
  tree->Branch("jetVeto30_JetAbsoluteFlavMapUp", &in->jetVeto30_JetAbsoluteFlavMapUp);
  tree->Branch("jetVeto30_JetAbsoluteMPFBiasDown", &in->jetVeto30_JetAbsoluteMPFBiasDown);
  tree->Branch("jetVeto30_JetAbsoluteMPFBiasUp", &in->jetVeto30_JetAbsoluteMPFBiasUp);
  tree->Branch("jetVeto30_JetAbsoluteScaleDown", &in->jetVeto30_JetAbsoluteScaleDown);
  tree->Branch("jetVeto30_JetAbsoluteScaleUp", &in->jetVeto30_JetAbsoluteScaleUp);
  tree->Branch("jetVeto30_JetAbsoluteStatDown", &in->jetVeto30_JetAbsoluteStatDown);
  tree->Branch("jetVeto30_JetAbsoluteStatUp", &in->jetVeto30_JetAbsoluteStatUp);
  tree->Branch("jetVeto30_JetClosureDown", &in->jetVeto30_JetClosureDown);
  tree->Branch("jetVeto30_JetClosureUp", &in->jetVeto30_JetClosureUp);
  tree->Branch("jetVeto30_JetEnDown", &in->jetVeto30_JetEnDown);
  tree->Branch("jetVeto30_JetFlavorQCDDown", &in->jetVeto30_JetFlavorQCDDown);
  tree->Branch("jetVeto30_JetFlavorQCDUp", &in->jetVeto30_JetFlavorQCDUp);
  tree->Branch("jetVeto30_JetFragmentationDown", &in->jetVeto30_JetFragmentationDown);
  tree->Branch("jetVeto30_JetFragmentationUp", &in->jetVeto30_JetFragmentationUp);
  tree->Branch("jetVeto30_JetPileUpDataMCDown", &in->jetVeto30_JetPileUpDataMCDown);
  tree->Branch("jetVeto30_JetPileUpDataMCUp", &in->jetVeto30_JetPileUpDataMCUp);
  tree->Branch("jetVeto30_JetPileUpPtBBDown", &in->jetVeto30_JetPileUpPtBBDown);
  tree->Branch("jetVeto30_JetPileUpPtBBUp", &in->jetVeto30_JetPileUpPtBBUp);
  tree->Branch("jetVeto30_JetPileUpPtEC1Down", &in->jetVeto30_JetPileUpPtEC1Down);
  tree->Branch("jetVeto30_JetPileUpPtEC1Up", &in->jetVeto30_JetPileUpPtEC1Up);
  tree->Branch("jetVeto30_JetPileUpPtEC2Down", &in->jetVeto30_JetPileUpPtEC2Down);
  tree->Branch("jetVeto30_JetPileUpPtEC2Up", &in->jetVeto30_JetPileUpPtEC2Up);
  tree->Branch("jetVeto30_JetPileUpPtHFDown", &in->jetVeto30_JetPileUpPtHFDown);
  tree->Branch("jetVeto30_JetPileUpPtHFUp", &in->jetVeto30_JetPileUpPtHFUp);
  tree->Branch("jetVeto30_JetPileUpPtRefDown", &in->jetVeto30_JetPileUpPtRefDown);
  tree->Branch("jetVeto30_JetPileUpPtRefUp", &in->jetVeto30_JetPileUpPtRefUp);
  tree->Branch("jetVeto30_JetRelativeBalDown", &in->jetVeto30_JetRelativeBalDown);
  tree->Branch("jetVeto30_JetRelativeBalUp", &in->jetVeto30_JetRelativeBalUp);
  tree->Branch("jetVeto30_JetRelativeFSRDown", &in->jetVeto30_JetRelativeFSRDown);
  tree->Branch("jetVeto30_JetRelativeFSRUp", &in->jetVeto30_JetRelativeFSRUp);
  tree->Branch("jetVeto30_JetRelativeJEREC1Down", &in->jetVeto30_JetRelativeJEREC1Down);
  tree->Branch("jetVeto30_JetRelativeJEREC1Up", &in->jetVeto30_JetRelativeJEREC1Up);
  tree->Branch("jetVeto30_JetRelativeJEREC2Down", &in->jetVeto30_JetRelativeJEREC2Down);
  tree->Branch("jetVeto30_JetRelativeJEREC2Up", &in->jetVeto30_JetRelativeJEREC2Up);
  tree->Branch("jetVeto30_JetRelativeJERHFDown", &in->jetVeto30_JetRelativeJERHFDown);
  tree->Branch("jetVeto30_JetRelativeJERHFUp", &in->jetVeto30_JetRelativeJERHFUp);
  tree->Branch("jetVeto30_JetRelativePtBBDown", &in->jetVeto30_JetRelativePtBBDown);
  tree->Branch("jetVeto30_JetRelativePtBBUp", &in->jetVeto30_JetRelativePtBBUp);
  tree->Branch("jetVeto30_JetRelativePtEC1Down", &in->jetVeto30_JetRelativePtEC1Down);
  tree->Branch("jetVeto30_JetRelativePtEC1Up", &in->jetVeto30_JetRelativePtEC1Up);
  tree->Branch("jetVeto30_JetRelativePtEC2Down", &in->jetVeto30_JetRelativePtEC2Down);
  tree->Branch("jetVeto30_JetRelativePtEC2Up", &in->jetVeto30_JetRelativePtEC2Up);
  tree->Branch("jetVeto30_JetRelativePtHFDown", &in->jetVeto30_JetRelativePtHFDown);
  tree->Branch("jetVeto30_JetRelativePtHFUp", &in->jetVeto30_JetRelativePtHFUp);
  tree->Branch("jetVeto30_JetRelativeSampleDown", &in->jetVeto30_JetRelativeSampleDown);
  tree->Branch("jetVeto30_JetRelativeSampleUp", &in->jetVeto30_JetRelativeSampleUp);
  tree->Branch("jetVeto30_JetRelativeStatECDown", &in->jetVeto30_JetRelativeStatECDown);
  tree->Branch("jetVeto30_JetRelativeStatECUp", &in->jetVeto30_JetRelativeStatECUp);
  tree->Branch("jetVeto30_JetRelativeStatFSRDown", &in->jetVeto30_JetRelativeStatFSRDown);
  tree->Branch("jetVeto30_JetRelativeStatFSRUp", &in->jetVeto30_JetRelativeStatFSRUp);
  tree->Branch("jetVeto30_JetRelativeStatHFDown", &in->jetVeto30_JetRelativeStatHFDown);
  tree->Branch("jetVeto30_JetRelativeStatHFUp", &in->jetVeto30_JetRelativeStatHFUp);
  tree->Branch("jetVeto30_JetSinglePionECALDown", &in->jetVeto30_JetSinglePionECALDown);
  tree->Branch("jetVeto30_JetSinglePionECALUp", &in->jetVeto30_JetSinglePionECALUp);
  tree->Branch("jetVeto30_JetSinglePionHCALDown", &in->jetVeto30_JetSinglePionHCALDown);
  tree->Branch("jetVeto30_JetSinglePionHCALUp", &in->jetVeto30_JetSinglePionHCALUp);
  tree->Branch("jetVeto30_JetTimePtEtaDown", &in->jetVeto30_JetTimePtEtaDown);
  tree->Branch("jetVeto30_JetTimePtEtaUp", &in->jetVeto30_JetTimePtEtaUp);
  tree->Branch("jetVeto30_JetTotalDown", &in->jetVeto30_JetTotalDown);
  tree->Branch("jetVeto30_JetTotalUp", &in->jetVeto30_JetTotalUp);

  tree->Branch("vbfMassWoNoisyJets_JetEta0to3Down", &in->vbfMassWoNoisyJets_JetEta0to3Down);
  tree->Branch("vbfMassWoNoisyJets_JetEta0to3Up", &in->vbfMassWoNoisyJets_JetEta0to3Up);
  tree->Branch("vbfMassWoNoisyJets_JetEta0to5Down", &in->vbfMassWoNoisyJets_JetEta0to5Down);
  tree->Branch("vbfMassWoNoisyJets_JetEta0to5Up", &in->vbfMassWoNoisyJets_JetEta0to5Up);
  tree->Branch("vbfMassWoNoisyJets_JetEta3to5Down", &in->vbfMassWoNoisyJets_JetEta3to5Down);
  tree->Branch("vbfMassWoNoisyJets_JetEta3to5Up", &in->vbfMassWoNoisyJets_JetEta3to5Up);
  tree->Branch("vbfMassWoNoisyJets_JetRelativeSampleDown", &in->vbfMassWoNoisyJets_JetRelativeSampleDown);
  tree->Branch("vbfMassWoNoisyJets_JetRelativeSampleUp", &in->vbfMassWoNoisyJets_JetRelativeSampleUp);
  tree->Branch("vbfMassWoNoisyJets_JetTotalDown", &in->vbfMassWoNoisyJets_JetTotalDown);
  tree->Branch("vbfMassWoNoisyJets_JetTotalUp", &in->vbfMassWoNoisyJets_JetTotalUp);
  tree->Branch("vbfMass_JetAbsoluteFlavMapDown", &in->vbfMass_JetAbsoluteFlavMapDown);
  tree->Branch("vbfMass_JetAbsoluteFlavMapUp", &in->vbfMass_JetAbsoluteFlavMapUp);
  tree->Branch("vbfMass_JetAbsoluteMPFBiasDown", &in->vbfMass_JetAbsoluteMPFBiasDown);
  tree->Branch("vbfMass_JetAbsoluteMPFBiasUp", &in->vbfMass_JetAbsoluteMPFBiasUp);
  tree->Branch("vbfMass_JetAbsoluteScaleDown", &in->vbfMass_JetAbsoluteScaleDown);
  tree->Branch("vbfMass_JetAbsoluteScaleUp", &in->vbfMass_JetAbsoluteScaleUp);
  tree->Branch("vbfMass_JetAbsoluteStatDown", &in->vbfMass_JetAbsoluteStatDown);
  tree->Branch("vbfMass_JetAbsoluteStatUp", &in->vbfMass_JetAbsoluteStatUp);
  tree->Branch("vbfMass_JetClosureDown", &in->vbfMass_JetClosureDown);
  tree->Branch("vbfMass_JetClosureUp", &in->vbfMass_JetClosureUp);
  tree->Branch("vbfMass_JetFlavorQCDDown", &in->vbfMass_JetFlavorQCDDown);
  tree->Branch("vbfMass_JetFlavorQCDUp", &in->vbfMass_JetFlavorQCDUp);
  tree->Branch("vbfMass_JetFragmentationDown", &in->vbfMass_JetFragmentationDown);
  tree->Branch("vbfMass_JetFragmentationUp", &in->vbfMass_JetFragmentationUp);
  tree->Branch("vbfMass_JetPileUpDataMCDown", &in->vbfMass_JetPileUpDataMCDown);
  tree->Branch("vbfMass_JetPileUpDataMCUp", &in->vbfMass_JetPileUpDataMCUp);
  tree->Branch("vbfMass_JetPileUpPtBBDown", &in->vbfMass_JetPileUpPtBBDown);
  tree->Branch("vbfMass_JetPileUpPtBBUp", &in->vbfMass_JetPileUpPtBBUp);
  tree->Branch("vbfMass_JetPileUpPtEC1Down", &in->vbfMass_JetPileUpPtEC1Down);
  tree->Branch("vbfMass_JetPileUpPtEC1Up", &in->vbfMass_JetPileUpPtEC1Up);
  tree->Branch("vbfMass_JetPileUpPtEC2Down", &in->vbfMass_JetPileUpPtEC2Down);
  tree->Branch("vbfMass_JetPileUpPtEC2Up", &in->vbfMass_JetPileUpPtEC2Up);
  tree->Branch("vbfMass_JetPileUpPtHFDown", &in->vbfMass_JetPileUpPtHFDown);
  tree->Branch("vbfMass_JetPileUpPtHFUp", &in->vbfMass_JetPileUpPtHFUp);
  tree->Branch("vbfMass_JetPileUpPtRefDown", &in->vbfMass_JetPileUpPtRefDown);
  tree->Branch("vbfMass_JetPileUpPtRefUp", &in->vbfMass_JetPileUpPtRefUp);
  tree->Branch("vbfMass_JetRelativeBalDown", &in->vbfMass_JetRelativeBalDown);
  tree->Branch("vbfMass_JetRelativeBalUp", &in->vbfMass_JetRelativeBalUp);
  tree->Branch("vbfMass_JetRelativeFSRDown", &in->vbfMass_JetRelativeFSRDown);
  tree->Branch("vbfMass_JetRelativeFSRUp", &in->vbfMass_JetRelativeFSRUp);
  tree->Branch("vbfMass_JetRelativeJEREC1Down", &in->vbfMass_JetRelativeJEREC1Down);
  tree->Branch("vbfMass_JetRelativeJEREC1Up", &in->vbfMass_JetRelativeJEREC1Up);
  tree->Branch("vbfMass_JetRelativeJEREC2Down", &in->vbfMass_JetRelativeJEREC2Down);
  tree->Branch("vbfMass_JetRelativeJEREC2Up", &in->vbfMass_JetRelativeJEREC2Up);
  tree->Branch("vbfMass_JetRelativeJERHFDown", &in->vbfMass_JetRelativeJERHFDown);
  tree->Branch("vbfMass_JetRelativeJERHFUp", &in->vbfMass_JetRelativeJERHFUp);
  tree->Branch("vbfMass_JetRelativePtBBDown", &in->vbfMass_JetRelativePtBBDown);
  tree->Branch("vbfMass_JetRelativePtBBUp", &in->vbfMass_JetRelativePtBBUp);
  tree->Branch("vbfMass_JetRelativePtEC1Down", &in->vbfMass_JetRelativePtEC1Down);
  tree->Branch("vbfMass_JetRelativePtEC1Up", &in->vbfMass_JetRelativePtEC1Up);
  tree->Branch("vbfMass_JetRelativePtEC2Down", &in->vbfMass_JetRelativePtEC2Down);
  tree->Branch("vbfMass_JetRelativePtEC2Up", &in->vbfMass_JetRelativePtEC2Up);
  tree->Branch("vbfMass_JetRelativePtHFDown", &in->vbfMass_JetRelativePtHFDown);
  tree->Branch("vbfMass_JetRelativePtHFUp", &in->vbfMass_JetRelativePtHFUp);
  tree->Branch("vbfMass_JetRelativeSampleDown", &in->vbfMass_JetRelativeSampleDown);
  tree->Branch("vbfMass_JetRelativeSampleUp", &in->vbfMass_JetRelativeSampleUp);
  tree->Branch("vbfMass_JetRelativeStatECDown", &in->vbfMass_JetRelativeStatECDown);
  tree->Branch("vbfMass_JetRelativeStatECUp", &in->vbfMass_JetRelativeStatECUp);
  tree->Branch("vbfMass_JetRelativeStatFSRDown", &in->vbfMass_JetRelativeStatFSRDown);
  tree->Branch("vbfMass_JetRelativeStatFSRUp", &in->vbfMass_JetRelativeStatFSRUp);
  tree->Branch("vbfMass_JetRelativeStatHFDown", &in->vbfMass_JetRelativeStatHFDown);
  tree->Branch("vbfMass_JetRelativeStatHFUp", &in->vbfMass_JetRelativeStatHFUp);
  tree->Branch("vbfMass_JetSinglePionECALDown", &in->vbfMass_JetSinglePionECALDown);
  tree->Branch("vbfMass_JetSinglePionECALUp", &in->vbfMass_JetSinglePionECALUp);
  tree->Branch("vbfMass_JetSinglePionHCALDown", &in->vbfMass_JetSinglePionHCALDown);
  tree->Branch("vbfMass_JetSinglePionHCALUp", &in->vbfMass_JetSinglePionHCALUp);
  tree->Branch("vbfMass_JetTimePtEtaDown", &in->vbfMass_JetTimePtEtaDown);
  tree->Branch("vbfMass_JetTimePtEtaUp", &in->vbfMass_JetTimePtEtaUp);
  tree->Branch("vbfMass_JetTotalDown", &in->vbfMass_JetTotalDown);
  tree->Branch("vbfMass_JetTotalUp", &in->vbfMass_JetTotalUp);

  // 2016 placeholders
  tree->Branch("amcatNLO_weight", &placeholder, "amcatNLO_weight/F");
  tree->Branch("matchIsoMu22eta2p1_1", &placeholder);
  tree->Branch("matchIsoTkMu22eta2p1_1", &placeholder);
  tree->Branch("matchIsoMu22_1", &placeholder);
  tree->Branch("matchIsoTkMu22_1", &placeholder);
  tree->Branch("matchIsoMu19Tau20_1", &placeholder);
  tree->Branch("filterIsoMu22eta2p1_1", &placeholder);
  tree->Branch("filterIsoTkMu22eta2p1_1", &placeholder);
  tree->Branch("filterIsoMu22_1", &placeholder);
  tree->Branch("filterIsoTkMu22_1", &placeholder);
  tree->Branch("filterIsoMu19Tau20_1", &placeholder);
  tree->Branch("passIsoMu22eta2p1", &placeholder);
  tree->Branch("passIsoTkMu22eta2p1", &placeholder);
  tree->Branch("passIsoMu22", &placeholder);
  tree->Branch("passIsoTkMu22", &placeholder);
  tree->Branch("passIsoMu19Tau20", &placeholder);
  tree->Branch("matchIsoMu19Tau20_2", &placeholder);
  tree->Branch("filterIsoMu19Tau20_2", &placeholder);
}

#endif  // ROOT_SRC_MUTAU_TREE_2017_H_
