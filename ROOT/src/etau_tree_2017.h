// Copyright 2018 Tyler Mitchell

#ifndef ROOT_SRC_ETAU_TREE_2017_H_
#define ROOT_SRC_ETAU_TREE_2017_H_

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "../interface/etau_input_branches.h"
#include "./base_tree.h"
#include "RecoilCorrector.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

class etau_tree2017 : public virtual base_tree {
 private:
  TTree *tree, *original;
  etau_input_branches* in;
  bool isMC, isEmbed;
  std::vector<Int_t> good_events;
  TLorentzVector ele, tau, MET, MET_UESUp, MET_UESDown, MET_JESUp, MET_JESDown;

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
  etau_tree2017(TTree* orig, TTree* itree, bool isMC, bool isEmbed, Int_t rec);
  virtual ~etau_tree2017() {}
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
etau_tree2017::etau_tree2017(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec) : tree(itree),
                                                                                                  original(Original),
                                                                                                  in(new etau_input_branches(Original)),
                                                                                                  isMC(IsMC),
                                                                                                  isEmbed(IsEmbed),
                                                                                                  recoil(rec) {}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void etau_tree2017::do_skimming(TH1F* cutflow) {
  // declare variables for sorting
  ULong64_t evt_now(0);
  ULong64_t evt_before(1);
  int best_evt(-1);
  std::pair<float, float> eleCandidate, tauCandidate;

  Int_t nevt = (Int_t)original->GetEntries();
  for (auto ievt = 0; ievt < nevt; ievt++) {
    original->GetEntry(ievt);
    evt_now = in->evt;

    // TLorentzVector ele, tau;
    ele.SetPtEtaPhiM(in->ePt, in->eEta, in->ePhi, in->eMass);
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

    float el_pt_min(25.), tau_pt_min(23.);

    cutflow->Fill(1., 1.);
    // apply event selection

    auto Ele27 = in->eMatchesEle27Filter && in->eMatchesEle27Path && in->Ele27WPTightPass;
    auto Ele32 = in->eMatchesEle32Filter && in->eMatchesEle32Path && in->Ele32WPTightPass;
    auto Ele35 = in->eMatchesEle35Filter && in->eMatchesEle35Path && in->Ele35WPTightPass;
    auto Cross = in->eMatchesEle24Tau30Filter && in->eMatchesEle24Tau30Path && in->Ele24LooseTau30Pass && in->tMatchesEle24Tau30Path && in->tMatchesEle24Tau30Filter;

    if (isEmbed || (Ele27 || Ele32 || Ele35 || Cross))
      cutflow->Fill(2., 1.);
    else
      continue;

    if (!isEmbed || (in->Ele27WPTightPass || in->Ele32WPTightPass || in->Ele35WPTightPass || in->Ele24LooseTau30Pass))
      cutflow->Fill(3., 1.);
    else
      continue;

    if (in->ePt > el_pt_min && fabs(in->eEta) < 2.1 && fabs(in->ePVDZ) < 0.2 && fabs(in->ePVDXY) < 0.045)
      cutflow->Fill(4., 1.);  // electron kinematic selection
    else
      continue;

    if (in->eMVANoisoWP90 && in->ePassesConversionVeto && in->eMissingHits < 2)
      cutflow->Fill(5., 1.);  // electron quality selection
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

    if (in->tAgainstMuonLoose3 > 0.5 && in->tAgainstElectronTightMVA6 > 0.5)
      cutflow->Fill(8., 1.);  // tau against leptons
    else
      continue;

    if (in->muVetoZTTp001dxyzR0 == 0 && in->eVetoZTTp001dxyzR0 < 2 && in->dielectronVeto == 0)
      cutflow->Fill(9., 1.);  // vetos
    else
      continue;

    // implement new sorting per
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
    if (evt_now != evt_before) {  // new event, save the tau candidates
      //   since it is new event, do we have the best entry to save? If yes, save it!
      if (best_evt > -1)
        good_events.push_back(best_evt);

      //  this is a new event, so the first tau pair is the best! :)
      best_evt = ievt;
      eleCandidate = std::make_pair(in->ePt, in->eIsoDB03);
      tauCandidate = std::make_pair(in->tPt, in->tRerunMVArun2v2DBoldDMwLTraw);
    } else {  // not a new event
      std::pair<float, float> currEleCandidate(in->ePt, in->eIsoDB03);
      std::pair<float, float> currTauCandidate(in->tPt, in->tRerunMVArun2v2DBoldDMwLTraw);

      // clause 1, select the pair that has most isolated tau lepton 1
      if (currEleCandidate.second - eleCandidate.second > 0.0001) best_evt = ievt;

      // check if the first tau is the same, and if so - move to clause 2
      if (fabs(currEleCandidate.second - eleCandidate.second) < 0.0001) {
        // pick up  the pair with the highest pT of the first candidate
        if (currEleCandidate.first - eleCandidate.first > 0.0001) best_evt = ievt;
        if (fabs(currEleCandidate.first - eleCandidate.first) < 0.0001) {
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

Float_t etau_tree2017::do_tes_met_corr(Float_t decayMode, Float_t sf1, Float_t sf2, Float_t sf3, TLorentzVector& met, TLorentzVector tau) {
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
TTree* etau_tree2017::fill_tree(RecoilCorrector recoilPFMetCorrector) {
  set_branches();  // get all the branches set up

  // loop through all events pasing skimming/sorting
  for (auto& ievt : good_events) {
    original->GetEntry(ievt);

    Run = in->run;
    Lumi = in->lumi;

    // convert from Float_t in FSA to Int_t for analyzer
    gen_match_1 = in->eZTTGenMatching;
    gen_match_2 = in->tZTTGenMatching;
    // njets = in->jetVeto30;
    njets = in->jetVeto30WoNoisyJets_JetEC2Down;  // temporary
    nbtag = in->bjetDeepCSVVeto20Medium_2017_DR0p5;
    njetspt20 = in->jetVeto20;

    // TLorentzVector ele, tau;
    ele.SetPtEtaPhiM(in->ePt, in->eEta, in->ePhi, in->eMass);
    tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);

    met_px = in->type1_pfMetEt * cos(in->type1_pfMetPhi);
    met_py = in->type1_pfMetEt * sin(in->type1_pfMetPhi);

    extraelec_veto = in->eVetoZTTp001dxyzR0 > 1;
    extramuon_veto = in->muVetoZTTp001dxyzR0 > 0;
    dilepton_veto = in->dielectronVeto > 0;

    // TLorentzVector ele, tau;
    ele.SetPtEtaPhiM(in->ePt, in->eEta, in->ePhi, in->eMass);
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
          MET.Px(),           // uncorrected type I pf met px (float)
          MET.Py(),           // uncorrected type I pf met py (float)
          in->genpX,          // generator Z/W/Higgs px (float)
          in->genpY,          // generator Z/W/Higgs py (float)
          in->vispX,          // generator visible Z/W/Higgs px (float)
          in->vispY,          // generator visible Z/W/Higgs py (float)
          in->jetVeto30 + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,       // corrected type I pf met px (float)
          pfmetcorr_ey);      // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30 + 1,    // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30 + 1,    // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30 + 1,      // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30 + 1,      // number of jets (hadronic jet multiplicity) (int)
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
          in->jetVeto30,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          in->jetVeto30,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          in->jetVeto30,          // number of jets (hadronic jet multiplicity) (int)
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

    m_1 = ele.M();
    px_1 = ele.Px();
    py_1 = ele.Py();
    pz_1 = ele.Pz();
    e_1 = ele.E();
    pt_1 = ele.Pt();
    phi_1 = ele.Phi();
    eta_1 = ele.Eta();
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
void etau_tree2017::set_branches() {
  // new branches
  tree->SetBranchAddress("run", &Run);
  tree->SetBranchAddress("lumi", &Lumi);
  tree->SetBranchAddress("gen_match_1", &gen_match_1);
  tree->SetBranchAddress("gen_match_2", &gen_match_2);
  tree->SetBranchAddress("njets", &njets);
  tree->SetBranchAddress("nbtag", &nbtag);
  tree->SetBranchAddress("njetspt20", &njetspt20);
  tree->SetBranchAddress("met_px", &met_px);
  tree->SetBranchAddress("met_py", &met_py);
  tree->SetBranchAddress("extraelec_veto", &extraelec_veto);
  tree->SetBranchAddress("extramuon_veto", &extramuon_veto);
  tree->SetBranchAddress("dilepton_veto", &dilepton_veto);
  tree->SetBranchAddress("pfmetcorr_ex", &pfmetcorr_ex);
  tree->SetBranchAddress("pfmetcorr_ey", &pfmetcorr_ey);
  tree->SetBranchAddress("pfmetcorr_ex_UESUp", &pfmetcorr_ex_UESUp);
  tree->SetBranchAddress("pfmetcorr_ey_UESUp", &pfmetcorr_ey_UESUp);
  tree->SetBranchAddress("pfmetcorr_ex_UESDown", &pfmetcorr_ex_UESDown);
  tree->SetBranchAddress("pfmetcorr_ey_UESDown", &pfmetcorr_ey_UESDown);
  tree->SetBranchAddress("pfmetcorr_ex_JESUp", &pfmetcorr_ex_JESUp);
  tree->SetBranchAddress("pfmetcorr_ey_JESUp", &pfmetcorr_ey_JESUp);
  tree->SetBranchAddress("pfmetcorr_ex_JESDown", &pfmetcorr_ex_JESDown);
  tree->SetBranchAddress("pfmetcorr_ey_JESDown", &pfmetcorr_ey_JESDown);
  tree->SetBranchAddress("met", &met);
  tree->SetBranchAddress("metphi", &metphi);
  tree->SetBranchAddress("met_px", &met_px);
  tree->SetBranchAddress("met_py", &met_py);
  tree->SetBranchAddress("met_JESUp", &met_JESUp);
  tree->SetBranchAddress("met_JESDown", &met_JESDown);
  tree->SetBranchAddress("met_UESUp", &met_UESUp);
  tree->SetBranchAddress("met_UESDown", &met_UESDown);
  tree->SetBranchAddress("metphi_JESUp", &metphi_JESUp);
  tree->SetBranchAddress("metphi_JESDown", &metphi_JESDown);
  tree->SetBranchAddress("metphi_UESUp", &metphi_UESUp);
  tree->SetBranchAddress("metphi_UESDown", &metphi_UESDown);
  tree->SetBranchAddress("m_1", &m_1);
  tree->SetBranchAddress("px_1", &px_1);
  tree->SetBranchAddress("py_1", &py_1);
  tree->SetBranchAddress("pz_1", &pz_1);
  tree->SetBranchAddress("e_1", &e_1);
  tree->SetBranchAddress("pt_1", &pt_1);
  tree->SetBranchAddress("phi_1", &phi_1);
  tree->SetBranchAddress("eta_1", &eta_1);
  tree->SetBranchAddress("m_2", &m_2);
  tree->SetBranchAddress("px_2", &px_2);
  tree->SetBranchAddress("py_2", &py_2);
  tree->SetBranchAddress("pz_2", &pz_2);
  tree->SetBranchAddress("e_2", &e_2);
  tree->SetBranchAddress("pt_2", &pt_2);
  tree->SetBranchAddress("phi_2", &phi_2);
  tree->SetBranchAddress("eta_2", &eta_2);

  // SVFit and MELA branches
  tree->SetBranchAddress("q_1", &in->eCharge);
  tree->SetBranchAddress("q_2", &in->tCharge);
  tree->SetBranchAddress("jeta_1", &in->j1eta);
  tree->SetBranchAddress("jphi_1", &in->j1phi);
  tree->SetBranchAddress("jeta_2", &in->j2eta);
  tree->SetBranchAddress("jphi_2", &in->j2phi);
  
  // copy the rest
  tree->SetBranchAddress("DoubleMediumHPSTau35Pass", &in->DoubleMediumHPSTau35Pass);
  tree->SetBranchAddress("DoubleMediumHPSTau35TightIDPass", &in->DoubleMediumHPSTau35TightIDPass);
  tree->SetBranchAddress("DoubleMediumHPSTau40Pass", &in->DoubleMediumHPSTau40Pass);
  tree->SetBranchAddress("DoubleMediumHPSTau40TightIDPass", &in->DoubleMediumHPSTau40TightIDPass);
  tree->SetBranchAddress("DoubleMediumTau35Pass", &in->DoubleMediumTau35Pass);
  tree->SetBranchAddress("DoubleMediumTau35TightIDPass", &in->DoubleMediumTau35TightIDPass);
  tree->SetBranchAddress("DoubleMediumTau40Pass", &in->DoubleMediumTau40Pass);
  tree->SetBranchAddress("DoubleMediumTau40TightIDPass", &in->DoubleMediumTau40TightIDPass);
  tree->SetBranchAddress("DoubleTightHPSTau35Pass", &in->DoubleTightHPSTau35Pass);
  tree->SetBranchAddress("DoubleTightHPSTau35TightIDPass", &in->DoubleTightHPSTau35TightIDPass);
  tree->SetBranchAddress("DoubleTightHPSTau40Pass", &in->DoubleTightHPSTau40Pass);
  tree->SetBranchAddress("DoubleTightHPSTau40TightIDPass", &in->DoubleTightHPSTau40TightIDPass);
  tree->SetBranchAddress("DoubleTightTau35Pass", &in->DoubleTightTau35Pass);
  tree->SetBranchAddress("DoubleTightTau35TightIDPass", &in->DoubleTightTau35TightIDPass);
  tree->SetBranchAddress("DoubleTightTau40Pass", &in->DoubleTightTau40Pass);
  tree->SetBranchAddress("DoubleTightTau40TightIDPass", &in->DoubleTightTau40TightIDPass);
  tree->SetBranchAddress("Ele24LooseHPSTau30Pass", &in->Ele24LooseHPSTau30Pass);
  tree->SetBranchAddress("Ele24LooseHPSTau30TightIDPass", &in->Ele24LooseHPSTau30TightIDPass);
  tree->SetBranchAddress("Ele24LooseTau30Pass", &in->Ele24LooseTau30Pass);
  tree->SetBranchAddress("Ele24LooseTau30TightIDPass", &in->Ele24LooseTau30TightIDPass);
  tree->SetBranchAddress("Ele27WPTightPass", &in->Ele27WPTightPass);
  tree->SetBranchAddress("Ele32WPTightPass", &in->Ele32WPTightPass);
  tree->SetBranchAddress("Ele35WPTightPass", &in->Ele35WPTightPass);
  tree->SetBranchAddress("Ele38WPTightPass", &in->Ele38WPTightPass);
  tree->SetBranchAddress("Ele40WPTightPass", &in->Ele40WPTightPass);
  tree->SetBranchAddress("EmbPtWeight", &in->EmbPtWeight);
  tree->SetBranchAddress("Eta", &in->Eta);
  tree->SetBranchAddress("Flag_BadChargedCandidateFilter", &in->Flag_BadChargedCandidateFilter);
  tree->SetBranchAddress("Flag_BadPFMuonFilter", &in->Flag_BadPFMuonFilter);
  tree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &in->Flag_EcalDeadCellTriggerPrimitiveFilter);
  tree->SetBranchAddress("Flag_HBHENoiseFilter", &in->Flag_HBHENoiseFilter);
  tree->SetBranchAddress("Flag_HBHENoiseIsoFilter", &in->Flag_HBHENoiseIsoFilter);
  tree->SetBranchAddress("Flag_badMuons", &in->Flag_badMuons);
  tree->SetBranchAddress("Flag_duplicateMuons", &in->Flag_duplicateMuons);
  tree->SetBranchAddress("Flag_ecalBadCalibFilter", &in->Flag_ecalBadCalibFilter);
  tree->SetBranchAddress("Flag_eeBadScFilter", &in->Flag_eeBadScFilter);
  tree->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &in->Flag_globalSuperTightHalo2016Filter);
  tree->SetBranchAddress("Flag_globalTightHalo2016Filter", &in->Flag_globalTightHalo2016Filter);
  tree->SetBranchAddress("Flag_goodVertices", &in->Flag_goodVertices);
  tree->SetBranchAddress("GenWeight", &in->GenWeight);
  tree->SetBranchAddress("Ht", &in->Ht);
  tree->SetBranchAddress("IsoMu24Pass", &in->IsoMu24Pass);
  tree->SetBranchAddress("IsoMu27Pass", &in->IsoMu27Pass);
  tree->SetBranchAddress("LT", &in->LT);
  tree->SetBranchAddress("Mass", &in->Mass);
  tree->SetBranchAddress("MassError", &in->MassError);
  tree->SetBranchAddress("MassErrord1", &in->MassErrord1);
  tree->SetBranchAddress("MassErrord2", &in->MassErrord2);
  tree->SetBranchAddress("MassErrord3", &in->MassErrord3);
  tree->SetBranchAddress("MassErrord4", &in->MassErrord4);
  tree->SetBranchAddress("Mt", &in->Mt);
  tree->SetBranchAddress("Mu20LooseHPSTau27Pass", &in->Mu20LooseHPSTau27Pass);
  tree->SetBranchAddress("Mu20LooseHPSTau27TightIDPass", &in->Mu20LooseHPSTau27TightIDPass);
  tree->SetBranchAddress("Mu20LooseTau27Pass", &in->Mu20LooseTau27Pass);
  tree->SetBranchAddress("Mu20LooseTau27TightIDPass", &in->Mu20LooseTau27TightIDPass);
  tree->SetBranchAddress("Mu50Pass", &in->Mu50Pass);
  tree->SetBranchAddress("NUP", &in->NUP);
  tree->SetBranchAddress("Phi", &in->Phi);
  tree->SetBranchAddress("Pt", &in->Pt);
  tree->SetBranchAddress("Rivet_VEta", &in->Rivet_VEta);
  tree->SetBranchAddress("Rivet_VPt", &in->Rivet_VPt);
  tree->SetBranchAddress("Rivet_errorCode", &in->Rivet_errorCode);
  tree->SetBranchAddress("Rivet_higgsEta", &in->Rivet_higgsEta);
  tree->SetBranchAddress("Rivet_higgsPt", &in->Rivet_higgsPt);
  tree->SetBranchAddress("Rivet_nJets25", &in->Rivet_nJets25);
  tree->SetBranchAddress("Rivet_nJets30", &in->Rivet_nJets30);
  tree->SetBranchAddress("Rivet_p4decay_VEta", &in->Rivet_p4decay_VEta);
  tree->SetBranchAddress("Rivet_p4decay_VPt", &in->Rivet_p4decay_VPt);
  tree->SetBranchAddress("Rivet_prodMode", &in->Rivet_prodMode);
  tree->SetBranchAddress("Rivet_stage0_cat", &in->Rivet_stage0_cat);
  tree->SetBranchAddress("Rivet_stage1_cat_pTjet25GeV", &in->Rivet_stage1_cat_pTjet25GeV);
  tree->SetBranchAddress("Rivet_stage1_cat_pTjet30GeV", &in->Rivet_stage1_cat_pTjet30GeV);
  tree->SetBranchAddress("Rivet_stage1p1_cat", &in->Rivet_stage1p1_cat);
  tree->SetBranchAddress("VBFDoubleLooseHPSTau20Pass", &in->VBFDoubleLooseHPSTau20Pass);
  tree->SetBranchAddress("VBFDoubleLooseTau20Pass", &in->VBFDoubleLooseTau20Pass);
  tree->SetBranchAddress("VBFDoubleMediumHPSTau20Pass", &in->VBFDoubleMediumHPSTau20Pass);
  tree->SetBranchAddress("VBFDoubleMediumTau20Pass", &in->VBFDoubleMediumTau20Pass);
  tree->SetBranchAddress("VBFDoubleTightHPSTau20Pass", &in->VBFDoubleTightHPSTau20Pass);
  tree->SetBranchAddress("VBFDoubleTightTau20Pass", &in->VBFDoubleTightTau20Pass);
  tree->SetBranchAddress("bjetDeepCSVVeto20Loose_2016_DR0p5", &in->bjetDeepCSVVeto20Loose_2016_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Loose_2017_DR0p5", &in->bjetDeepCSVVeto20Loose_2017_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Loose_2018_DR0p5", &in->bjetDeepCSVVeto20Loose_2018_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Medium_2016_DR0", &in->bjetDeepCSVVeto20Medium_2016_DR0);
  tree->SetBranchAddress("bjetDeepCSVVeto20Medium_2016_DR0p5", &in->bjetDeepCSVVeto20Medium_2016_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Medium_2017_DR0", &in->bjetDeepCSVVeto20Medium_2017_DR0);
  tree->SetBranchAddress("bjetDeepCSVVeto20Medium_2017_DR0p5", &in->bjetDeepCSVVeto20Medium_2017_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Medium_2018_DR0", &in->bjetDeepCSVVeto20Medium_2018_DR0);
  tree->SetBranchAddress("bjetDeepCSVVeto20Medium_2018_DR0p5", &in->bjetDeepCSVVeto20Medium_2018_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Tight_2016_DR0p5", &in->bjetDeepCSVVeto20Tight_2016_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Tight_2017_DR0p5", &in->bjetDeepCSVVeto20Tight_2017_DR0p5);
  tree->SetBranchAddress("bjetDeepCSVVeto20Tight_2018_DR0p5", &in->bjetDeepCSVVeto20Tight_2018_DR0p5);
  tree->SetBranchAddress("charge", &in->charge);
  tree->SetBranchAddress("dielectronVeto", &in->dielectronVeto);
  tree->SetBranchAddress("dimu9ele9Pass", &in->dimu9ele9Pass);
  tree->SetBranchAddress("dimuonVeto", &in->dimuonVeto);
  tree->SetBranchAddress("doubleE25Pass", &in->doubleE25Pass);
  tree->SetBranchAddress("doubleE33Pass", &in->doubleE33Pass);
  tree->SetBranchAddress("doubleE_23_12Pass", &in->doubleE_23_12Pass);
  tree->SetBranchAddress("doubleMuDZminMass3p8Pass", &in->doubleMuDZminMass3p8Pass);
  tree->SetBranchAddress("doubleMuDZminMass8Pass", &in->doubleMuDZminMass8Pass);
  tree->SetBranchAddress("doubleMuSingleEPass", &in->doubleMuSingleEPass);
  tree->SetBranchAddress("doubleTau35Pass", &in->doubleTau35Pass);
  tree->SetBranchAddress("doubleTauCmbIso35RegPass", &in->doubleTauCmbIso35RegPass);
  tree->SetBranchAddress("eCBIDLoose", &in->eCBIDLoose);
  tree->SetBranchAddress("eCBIDMedium", &in->eCBIDMedium);
  tree->SetBranchAddress("eCBIDTight", &in->eCBIDTight);
  tree->SetBranchAddress("eCBIDVeto", &in->eCBIDVeto);
  tree->SetBranchAddress("eCharge", &in->eCharge);
  tree->SetBranchAddress("eChargeIdLoose", &in->eChargeIdLoose);
  tree->SetBranchAddress("eChargeIdMed", &in->eChargeIdMed);
  tree->SetBranchAddress("eChargeIdTight", &in->eChargeIdTight);
  tree->SetBranchAddress("eComesFromHiggs", &in->eComesFromHiggs);
  tree->SetBranchAddress("eCorrectedEt", &in->eCorrectedEt);
  tree->SetBranchAddress("eE1x5", &in->eE1x5);
  tree->SetBranchAddress("eE2x5Max", &in->eE2x5Max);
  tree->SetBranchAddress("eE5x5", &in->eE5x5);
  tree->SetBranchAddress("eEcalIsoDR03", &in->eEcalIsoDR03);
  tree->SetBranchAddress("eEnergyError", &in->eEnergyError);
  tree->SetBranchAddress("eEnergyScaleDown", &in->eEnergyScaleDown);
  tree->SetBranchAddress("eEnergyScaleGainDown", &in->eEnergyScaleGainDown);
  tree->SetBranchAddress("eEnergyScaleGainUp", &in->eEnergyScaleGainUp);
  tree->SetBranchAddress("eEnergyScaleStatDown", &in->eEnergyScaleStatDown);
  tree->SetBranchAddress("eEnergyScaleStatUp", &in->eEnergyScaleStatUp);
  tree->SetBranchAddress("eEnergyScaleSystDown", &in->eEnergyScaleSystDown);
  tree->SetBranchAddress("eEnergyScaleSystUp", &in->eEnergyScaleSystUp);
  tree->SetBranchAddress("eEnergyScaleUp", &in->eEnergyScaleUp);
  tree->SetBranchAddress("eEnergySigmaDown", &in->eEnergySigmaDown);
  tree->SetBranchAddress("eEnergySigmaPhiDown", &in->eEnergySigmaPhiDown);
  tree->SetBranchAddress("eEnergySigmaPhiUp", &in->eEnergySigmaPhiUp);
  tree->SetBranchAddress("eEnergySigmaRhoDown", &in->eEnergySigmaRhoDown);
  tree->SetBranchAddress("eEnergySigmaRhoUp", &in->eEnergySigmaRhoUp);
  tree->SetBranchAddress("eEnergySigmaUp", &in->eEnergySigmaUp);
  tree->SetBranchAddress("eEta", &in->eEta);
  tree->SetBranchAddress("eGenCharge", &in->eGenCharge);
  tree->SetBranchAddress("eGenDirectPromptTauDecay", &in->eGenDirectPromptTauDecay);
  tree->SetBranchAddress("eGenEnergy", &in->eGenEnergy);
  tree->SetBranchAddress("eGenEta", &in->eGenEta);
  tree->SetBranchAddress("eGenIsPrompt", &in->eGenIsPrompt);
  tree->SetBranchAddress("eGenMotherPdgId", &in->eGenMotherPdgId);
  tree->SetBranchAddress("eGenParticle", &in->eGenParticle);
  tree->SetBranchAddress("eGenPdgId", &in->eGenPdgId);
  tree->SetBranchAddress("eGenPhi", &in->eGenPhi);
  tree->SetBranchAddress("eGenPrompt", &in->eGenPrompt);
  tree->SetBranchAddress("eGenPromptTauDecay", &in->eGenPromptTauDecay);
  tree->SetBranchAddress("eGenPt", &in->eGenPt);
  tree->SetBranchAddress("eGenTauDecay", &in->eGenTauDecay);
  tree->SetBranchAddress("eGenVZ", &in->eGenVZ);
  tree->SetBranchAddress("eGenVtxPVMatch", &in->eGenVtxPVMatch);
  tree->SetBranchAddress("eHadronicDepth1OverEm", &in->eHadronicDepth1OverEm);
  tree->SetBranchAddress("eHadronicDepth2OverEm", &in->eHadronicDepth2OverEm);
  tree->SetBranchAddress("eHadronicOverEM", &in->eHadronicOverEM);
  tree->SetBranchAddress("eHcalIsoDR03", &in->eHcalIsoDR03);
  tree->SetBranchAddress("eIP3D", &in->eIP3D);
  tree->SetBranchAddress("eIP3DErr", &in->eIP3DErr);
  tree->SetBranchAddress("eIsoDB03", &in->eIsoDB03);
  tree->SetBranchAddress("eJetArea", &in->eJetArea);
  tree->SetBranchAddress("eJetBtag", &in->eJetBtag);
  tree->SetBranchAddress("eJetDR", &in->eJetDR);
  tree->SetBranchAddress("eJetEtaEtaMoment", &in->eJetEtaEtaMoment);
  tree->SetBranchAddress("eJetEtaPhiMoment", &in->eJetEtaPhiMoment);
  tree->SetBranchAddress("eJetEtaPhiSpread", &in->eJetEtaPhiSpread);
  tree->SetBranchAddress("eJetHadronFlavour", &in->eJetHadronFlavour);
  tree->SetBranchAddress("eJetPFCISVBtag", &in->eJetPFCISVBtag);
  tree->SetBranchAddress("eJetPartonFlavour", &in->eJetPartonFlavour);
  tree->SetBranchAddress("eJetPhiPhiMoment", &in->eJetPhiPhiMoment);
  tree->SetBranchAddress("eJetPt", &in->eJetPt);
  tree->SetBranchAddress("eLowestMll", &in->eLowestMll);
  tree->SetBranchAddress("eMVAIsoWP80", &in->eMVAIsoWP80);
  tree->SetBranchAddress("eMVAIsoWP90", &in->eMVAIsoWP90);
  tree->SetBranchAddress("eMVAIsoWPHZZ", &in->eMVAIsoWPHZZ);
  tree->SetBranchAddress("eMVAIsoWPLoose", &in->eMVAIsoWPLoose);
  tree->SetBranchAddress("eMVANoisoWP80", &in->eMVANoisoWP80);
  tree->SetBranchAddress("eMVANoisoWP90", &in->eMVANoisoWP90);
  tree->SetBranchAddress("eMVANoisoWPLoose", &in->eMVANoisoWPLoose);
  tree->SetBranchAddress("eMass", &in->eMass);
  tree->SetBranchAddress("eMatchesEle24HPSTau30Filter", &in->eMatchesEle24HPSTau30Filter);
  tree->SetBranchAddress("eMatchesEle24HPSTau30Path", &in->eMatchesEle24HPSTau30Path);
  tree->SetBranchAddress("eMatchesEle24Tau30Filter", &in->eMatchesEle24Tau30Filter);
  tree->SetBranchAddress("eMatchesEle24Tau30Path", &in->eMatchesEle24Tau30Path);
  tree->SetBranchAddress("eMatchesEle25Filter", &in->eMatchesEle25Filter);
  tree->SetBranchAddress("eMatchesEle25Path", &in->eMatchesEle25Path);
  tree->SetBranchAddress("eMatchesEle27Filter", &in->eMatchesEle27Filter);
  tree->SetBranchAddress("eMatchesEle27Path", &in->eMatchesEle27Path);
  tree->SetBranchAddress("eMatchesEle32Filter", &in->eMatchesEle32Filter);
  tree->SetBranchAddress("eMatchesEle32Path", &in->eMatchesEle32Path);
  tree->SetBranchAddress("eMatchesEle35Filter", &in->eMatchesEle35Filter);
  tree->SetBranchAddress("eMatchesEle35Path", &in->eMatchesEle35Path);
  tree->SetBranchAddress("eMissingHits", &in->eMissingHits);
  tree->SetBranchAddress("eNearMuonVeto", &in->eNearMuonVeto);
  tree->SetBranchAddress("eNearestMuonDR", &in->eNearestMuonDR);
  tree->SetBranchAddress("eNearestZMass", &in->eNearestZMass);
  tree->SetBranchAddress("ePFChargedIso", &in->ePFChargedIso);
  tree->SetBranchAddress("ePFNeutralIso", &in->ePFNeutralIso);
  tree->SetBranchAddress("ePFPUChargedIso", &in->ePFPUChargedIso);
  tree->SetBranchAddress("ePFPhotonIso", &in->ePFPhotonIso);
  tree->SetBranchAddress("ePVDXY", &in->ePVDXY);
  tree->SetBranchAddress("ePVDZ", &in->ePVDZ);
  tree->SetBranchAddress("ePassesConversionVeto", &in->ePassesConversionVeto);
  tree->SetBranchAddress("ePhi", &in->ePhi);
  tree->SetBranchAddress("ePt", &in->ePt);
  tree->SetBranchAddress("eRelIso", &in->eRelIso);
  tree->SetBranchAddress("eRelPFIsoDB", &in->eRelPFIsoDB);
  tree->SetBranchAddress("eRelPFIsoRho", &in->eRelPFIsoRho);
  tree->SetBranchAddress("eRho", &in->eRho);
  tree->SetBranchAddress("eSCEnergy", &in->eSCEnergy);
  tree->SetBranchAddress("eSCEta", &in->eSCEta);
  tree->SetBranchAddress("eSCEtaWidth", &in->eSCEtaWidth);
  tree->SetBranchAddress("eSCPhi", &in->eSCPhi);
  tree->SetBranchAddress("eSCPhiWidth", &in->eSCPhiWidth);
  tree->SetBranchAddress("eSCPreshowerEnergy", &in->eSCPreshowerEnergy);
  tree->SetBranchAddress("eSCRawEnergy", &in->eSCRawEnergy);
  tree->SetBranchAddress("eSIP2D", &in->eSIP2D);
  tree->SetBranchAddress("eSIP3D", &in->eSIP3D);
  tree->SetBranchAddress("eSigmaIEtaIEta", &in->eSigmaIEtaIEta);
  tree->SetBranchAddress("eTrkIsoDR03", &in->eTrkIsoDR03);
  tree->SetBranchAddress("eVZ", &in->eVZ);
  tree->SetBranchAddress("eVetoHZZPt5", &in->eVetoHZZPt5);
  tree->SetBranchAddress("eVetoZTTp001dxyz", &in->eVetoZTTp001dxyz);
  tree->SetBranchAddress("eVetoZTTp001dxyzR0", &in->eVetoZTTp001dxyzR0);
  tree->SetBranchAddress("eZTTGenMatching", &in->eZTTGenMatching);
  tree->SetBranchAddress("e_t_DR", &in->e_t_DR);
  tree->SetBranchAddress("e_t_Mass", &in->e_t_Mass);
  tree->SetBranchAddress("e_t_doubleL1IsoTauMatch", &in->e_t_doubleL1IsoTauMatch);
  tree->SetBranchAddress("edeltaEtaSuperClusterTrackAtVtx", &in->edeltaEtaSuperClusterTrackAtVtx);
  tree->SetBranchAddress("edeltaPhiSuperClusterTrackAtVtx", &in->edeltaPhiSuperClusterTrackAtVtx);
  tree->SetBranchAddress("eeSuperClusterOverP", &in->eeSuperClusterOverP);
  tree->SetBranchAddress("eecalEnergy", &in->eecalEnergy);
  tree->SetBranchAddress("efBrem", &in->efBrem);
  tree->SetBranchAddress("etrackMomentumAtVtxP", &in->etrackMomentumAtVtxP);
  tree->SetBranchAddress("evt", &in->evt);
  tree->SetBranchAddress("genEta", &in->genEta);
  tree->SetBranchAddress("genHTT", &in->genHTT);
  tree->SetBranchAddress("genM", &in->genM);
  tree->SetBranchAddress("genMass", &in->genMass);
  tree->SetBranchAddress("genPhi", &in->genPhi);
  tree->SetBranchAddress("genpT", &in->genpT);
  tree->SetBranchAddress("genpX", &in->genpX);
  tree->SetBranchAddress("genpY", &in->genpY);
  tree->SetBranchAddress("isdata", &in->isdata);
  tree->SetBranchAddress("isembed", &in->isembed);
  tree->SetBranchAddress("j1csv", &in->j1csv);
  tree->SetBranchAddress("j1csvWoNoisyJets", &in->j1csvWoNoisyJets);
  tree->SetBranchAddress("j1eta", &in->j1eta);
  tree->SetBranchAddress("j1etaWoNoisyJets", &in->j1etaWoNoisyJets);
  tree->SetBranchAddress("j1hadronflavor", &in->j1hadronflavor);
  tree->SetBranchAddress("j1hadronflavorWoNoisyJets", &in->j1hadronflavorWoNoisyJets);
  tree->SetBranchAddress("j1phi", &in->j1phi);
  tree->SetBranchAddress("j1phiWoNoisyJets", &in->j1phiWoNoisyJets);
  tree->SetBranchAddress("j1pt", &in->j1pt);
  tree->SetBranchAddress("j1ptWoNoisyJets", &in->j1ptWoNoisyJets);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEC2Down", &in->j1ptWoNoisyJets_JetEC2Down);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEC2Up", &in->j1ptWoNoisyJets_JetEC2Up);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEta0to3Down", &in->j1ptWoNoisyJets_JetEta0to3Down);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEta0to3Up", &in->j1ptWoNoisyJets_JetEta0to3Up);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEta0to5Down", &in->j1ptWoNoisyJets_JetEta0to5Down);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEta0to5Up", &in->j1ptWoNoisyJets_JetEta0to5Up);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEta3to5Down", &in->j1ptWoNoisyJets_JetEta3to5Down);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetEta3to5Up", &in->j1ptWoNoisyJets_JetEta3to5Up);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetRelativeBalDown", &in->j1ptWoNoisyJets_JetRelativeBalDown);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetRelativeBalUp", &in->j1ptWoNoisyJets_JetRelativeBalUp);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetRelativeSampleDown", &in->j1ptWoNoisyJets_JetRelativeSampleDown);
  tree->SetBranchAddress("j1ptWoNoisyJets_JetRelativeSampleUp", &in->j1ptWoNoisyJets_JetRelativeSampleUp);
  tree->SetBranchAddress("j1pt_JetEC2Down", &in->j1pt_JetEC2Down);
  tree->SetBranchAddress("j1pt_JetEC2Up", &in->j1pt_JetEC2Up);
  tree->SetBranchAddress("j1pt_JetEta0to3Down", &in->j1pt_JetEta0to3Down);
  tree->SetBranchAddress("j1pt_JetEta0to3Up", &in->j1pt_JetEta0to3Up);
  tree->SetBranchAddress("j1pt_JetEta0to5Down", &in->j1pt_JetEta0to5Down);
  tree->SetBranchAddress("j1pt_JetEta0to5Up", &in->j1pt_JetEta0to5Up);
  tree->SetBranchAddress("j1pt_JetEta3to5Down", &in->j1pt_JetEta3to5Down);
  tree->SetBranchAddress("j1pt_JetEta3to5Up", &in->j1pt_JetEta3to5Up);
  tree->SetBranchAddress("j1pt_JetRelativeBalDown", &in->j1pt_JetRelativeBalDown);
  tree->SetBranchAddress("j1pt_JetRelativeBalUp", &in->j1pt_JetRelativeBalUp);
  tree->SetBranchAddress("j1pt_JetRelativeSampleDown", &in->j1pt_JetRelativeSampleDown);
  tree->SetBranchAddress("j1pt_JetRelativeSampleUp", &in->j1pt_JetRelativeSampleUp);
  tree->SetBranchAddress("j2csv", &in->j2csv);
  tree->SetBranchAddress("j2csvWoNoisyJets", &in->j2csvWoNoisyJets);
  tree->SetBranchAddress("j2eta", &in->j2eta);
  tree->SetBranchAddress("j2etaWoNoisyJets", &in->j2etaWoNoisyJets);
  tree->SetBranchAddress("j2hadronflavor", &in->j2hadronflavor);
  tree->SetBranchAddress("j2hadronflavorWoNoisyJets", &in->j2hadronflavorWoNoisyJets);
  tree->SetBranchAddress("j2phi", &in->j2phi);
  tree->SetBranchAddress("j2phiWoNoisyJets", &in->j2phiWoNoisyJets);
  tree->SetBranchAddress("j2pt", &in->j2pt);
  tree->SetBranchAddress("j2ptWoNoisyJets", &in->j2ptWoNoisyJets);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEC2Down", &in->j2ptWoNoisyJets_JetEC2Down);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEC2Up", &in->j2ptWoNoisyJets_JetEC2Up);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEta0to3Down", &in->j2ptWoNoisyJets_JetEta0to3Down);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEta0to3Up", &in->j2ptWoNoisyJets_JetEta0to3Up);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEta0to5Down", &in->j2ptWoNoisyJets_JetEta0to5Down);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEta0to5Up", &in->j2ptWoNoisyJets_JetEta0to5Up);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEta3to5Down", &in->j2ptWoNoisyJets_JetEta3to5Down);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetEta3to5Up", &in->j2ptWoNoisyJets_JetEta3to5Up);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetRelativeBalDown", &in->j2ptWoNoisyJets_JetRelativeBalDown);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetRelativeBalUp", &in->j2ptWoNoisyJets_JetRelativeBalUp);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetRelativeSampleDown", &in->j2ptWoNoisyJets_JetRelativeSampleDown);
  tree->SetBranchAddress("j2ptWoNoisyJets_JetRelativeSampleUp", &in->j2ptWoNoisyJets_JetRelativeSampleUp);
  tree->SetBranchAddress("j2pt_JetEC2Down", &in->j2pt_JetEC2Down);
  tree->SetBranchAddress("j2pt_JetEC2Up", &in->j2pt_JetEC2Up);
  tree->SetBranchAddress("j2pt_JetEta0to3Down", &in->j2pt_JetEta0to3Down);
  tree->SetBranchAddress("j2pt_JetEta0to3Up", &in->j2pt_JetEta0to3Up);
  tree->SetBranchAddress("j2pt_JetEta0to5Down", &in->j2pt_JetEta0to5Down);
  tree->SetBranchAddress("j2pt_JetEta0to5Up", &in->j2pt_JetEta0to5Up);
  tree->SetBranchAddress("j2pt_JetEta3to5Down", &in->j2pt_JetEta3to5Down);
  tree->SetBranchAddress("j2pt_JetEta3to5Up", &in->j2pt_JetEta3to5Up);
  tree->SetBranchAddress("j2pt_JetRelativeBalDown", &in->j2pt_JetRelativeBalDown);
  tree->SetBranchAddress("j2pt_JetRelativeBalUp", &in->j2pt_JetRelativeBalUp);
  tree->SetBranchAddress("j2pt_JetRelativeSampleDown", &in->j2pt_JetRelativeSampleDown);
  tree->SetBranchAddress("j2pt_JetRelativeSampleUp", &in->j2pt_JetRelativeSampleUp);
  tree->SetBranchAddress("jb1eta_2016", &in->jb1eta_2016);
  tree->SetBranchAddress("jb1eta_2017", &in->jb1eta_2017);
  tree->SetBranchAddress("jb1eta_2018", &in->jb1eta_2018);
  tree->SetBranchAddress("jb1hadronflavor_2016", &in->jb1hadronflavor_2016);
  tree->SetBranchAddress("jb1hadronflavor_2017", &in->jb1hadronflavor_2017);
  tree->SetBranchAddress("jb1hadronflavor_2018", &in->jb1hadronflavor_2018);
  tree->SetBranchAddress("jb1phi_2016", &in->jb1phi_2016);
  tree->SetBranchAddress("jb1phi_2017", &in->jb1phi_2017);
  tree->SetBranchAddress("jb1phi_2018", &in->jb1phi_2018);
  tree->SetBranchAddress("jb1pt_2016", &in->jb1pt_2016);
  tree->SetBranchAddress("jb1pt_2017", &in->jb1pt_2017);
  tree->SetBranchAddress("jb1pt_2018", &in->jb1pt_2018);
  tree->SetBranchAddress("jb2eta_2016", &in->jb2eta_2016);
  tree->SetBranchAddress("jb2eta_2017", &in->jb2eta_2017);
  tree->SetBranchAddress("jb2eta_2018", &in->jb2eta_2018);
  tree->SetBranchAddress("jb2hadronflavor_2016", &in->jb2hadronflavor_2016);
  tree->SetBranchAddress("jb2hadronflavor_2017", &in->jb2hadronflavor_2017);
  tree->SetBranchAddress("jb2hadronflavor_2018", &in->jb2hadronflavor_2018);
  tree->SetBranchAddress("jb2phi_2016", &in->jb2phi_2016);
  tree->SetBranchAddress("jb2phi_2017", &in->jb2phi_2017);
  tree->SetBranchAddress("jb2phi_2018", &in->jb2phi_2018);
  tree->SetBranchAddress("jb2pt_2016", &in->jb2pt_2016);
  tree->SetBranchAddress("jb2pt_2017", &in->jb2pt_2017);
  tree->SetBranchAddress("jb2pt_2018", &in->jb2pt_2018);
  tree->SetBranchAddress("jetVeto20", &in->jetVeto20);
  tree->SetBranchAddress("jetVeto20_JetEnDown", &in->jetVeto20_JetEnDown);
  tree->SetBranchAddress("jetVeto20_JetEnUp", &in->jetVeto20_JetEnUp);
  tree->SetBranchAddress("jetVeto30", &in->jetVeto30);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEC2Down", &in->jetVeto30WoNoisyJets_JetEC2Down);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEC2Up", &in->jetVeto30WoNoisyJets_JetEC2Up);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEta0to3Down", &in->jetVeto30WoNoisyJets_JetEta0to3Down);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEta0to3Up", &in->jetVeto30WoNoisyJets_JetEta0to3Up);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEta0to5Down", &in->jetVeto30WoNoisyJets_JetEta0to5Down);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEta0to5Up", &in->jetVeto30WoNoisyJets_JetEta0to5Up);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEta3to5Down", &in->jetVeto30WoNoisyJets_JetEta3to5Down);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetEta3to5Up", &in->jetVeto30WoNoisyJets_JetEta3to5Up);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetRelativeBalDownWoNoisyJets", &in->jetVeto30WoNoisyJets_JetRelativeBalDownWoNoisyJets);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetRelativeBalUpWoNoisyJets", &in->jetVeto30WoNoisyJets_JetRelativeBalUpWoNoisyJets);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetRelativeSampleDown", &in->jetVeto30WoNoisyJets_JetRelativeSampleDown);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetRelativeSampleUp", &in->jetVeto30WoNoisyJets_JetRelativeSampleUp);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetTotalDown", &in->jetVeto30WoNoisyJets_JetTotalDown);
  tree->SetBranchAddress("jetVeto30WoNoisyJets_JetTotalUp", &in->jetVeto30WoNoisyJets_JetTotalUp);
  tree->SetBranchAddress("jetVeto30_JetEC2Down", &in->jetVeto30_JetEC2Down);
  tree->SetBranchAddress("jetVeto30_JetEC2Up", &in->jetVeto30_JetEC2Up);
  tree->SetBranchAddress("jetVeto30_JetEnDown", &in->jetVeto30_JetEnDown);
  tree->SetBranchAddress("jetVeto30_JetEnUp", &in->jetVeto30_JetEnUp);
  tree->SetBranchAddress("jetVeto30_JetEta0to3Down", &in->jetVeto30_JetEta0to3Down);
  tree->SetBranchAddress("jetVeto30_JetEta0to3Up", &in->jetVeto30_JetEta0to3Up);
  tree->SetBranchAddress("jetVeto30_JetEta0to5Down", &in->jetVeto30_JetEta0to5Down);
  tree->SetBranchAddress("jetVeto30_JetEta0to5Up", &in->jetVeto30_JetEta0to5Up);
  tree->SetBranchAddress("jetVeto30_JetEta3to5Down", &in->jetVeto30_JetEta3to5Down);
  tree->SetBranchAddress("jetVeto30_JetEta3to5Up", &in->jetVeto30_JetEta3to5Up);
  tree->SetBranchAddress("jetVeto30_JetRelativeBalDown", &in->jetVeto30_JetRelativeBalDown);
  tree->SetBranchAddress("jetVeto30_JetRelativeBalUp", &in->jetVeto30_JetRelativeBalUp);
  tree->SetBranchAddress("jetVeto30_JetRelativeSampleDown", &in->jetVeto30_JetRelativeSampleDown);
  tree->SetBranchAddress("jetVeto30_JetRelativeSampleUp", &in->jetVeto30_JetRelativeSampleUp);
  tree->SetBranchAddress("jetVeto30_JetTotalDown", &in->jetVeto30_JetTotalDown);
  tree->SetBranchAddress("jetVeto30_JetTotalUp", &in->jetVeto30_JetTotalUp);
  tree->SetBranchAddress("lumi", &in->lumi);
  tree->SetBranchAddress("metSig", &in->metSig);
  tree->SetBranchAddress("metcov00", &in->metcov00);
  tree->SetBranchAddress("metcov01", &in->metcov01);
  tree->SetBranchAddress("metcov10", &in->metcov10);
  tree->SetBranchAddress("metcov11", &in->metcov11);
  tree->SetBranchAddress("mu12e23DZPass", &in->mu12e23DZPass);
  tree->SetBranchAddress("mu12e23Pass", &in->mu12e23Pass);
  tree->SetBranchAddress("mu23e12DZPass", &in->mu23e12DZPass);
  tree->SetBranchAddress("mu23e12Pass", &in->mu23e12Pass);
  tree->SetBranchAddress("mu8diele12DZPass", &in->mu8diele12DZPass);
  tree->SetBranchAddress("mu8diele12Pass", &in->mu8diele12Pass);
  tree->SetBranchAddress("mu8e23DZPass", &in->mu8e23DZPass);
  tree->SetBranchAddress("mu8e23Pass", &in->mu8e23Pass);
  tree->SetBranchAddress("muGlbIsoVetoPt10", &in->muGlbIsoVetoPt10);
  tree->SetBranchAddress("muVeto5", &in->muVeto5);
  tree->SetBranchAddress("muVetoZTTp001dxyz", &in->muVetoZTTp001dxyz);
  tree->SetBranchAddress("muVetoZTTp001dxyzR0", &in->muVetoZTTp001dxyzR0);
  tree->SetBranchAddress("nTruePU", &in->nTruePU);
  tree->SetBranchAddress("npNLO", &in->npNLO);
  tree->SetBranchAddress("numGenJets", &in->numGenJets);
  tree->SetBranchAddress("nvtx", &in->nvtx);
  tree->SetBranchAddress("processID", &in->processID);
  tree->SetBranchAddress("puppiMetEt", &in->puppiMetEt);
  tree->SetBranchAddress("puppiMetPhi", &in->puppiMetPhi);
  tree->SetBranchAddress("pvChi2", &in->pvChi2);
  tree->SetBranchAddress("pvDX", &in->pvDX);
  tree->SetBranchAddress("pvDY", &in->pvDY);
  tree->SetBranchAddress("pvDZ", &in->pvDZ);
  tree->SetBranchAddress("pvIsFake", &in->pvIsFake);
  tree->SetBranchAddress("pvIsValid", &in->pvIsValid);
  tree->SetBranchAddress("pvNormChi2", &in->pvNormChi2);
  tree->SetBranchAddress("pvRho", &in->pvRho);
  tree->SetBranchAddress("pvX", &in->pvX);
  tree->SetBranchAddress("pvY", &in->pvY);
  tree->SetBranchAddress("pvZ", &in->pvZ);
  tree->SetBranchAddress("pvndof", &in->pvndof);
  tree->SetBranchAddress("raw_pfMetEt", &in->raw_pfMetEt);
  tree->SetBranchAddress("raw_pfMetPhi", &in->raw_pfMetPhi);
  tree->SetBranchAddress("recoilDaught", &in->recoilDaught);
  tree->SetBranchAddress("recoilWithMet", &in->recoilWithMet);
  tree->SetBranchAddress("rho", &in->rho);
  tree->SetBranchAddress("run", &in->run);
  tree->SetBranchAddress("singleE25eta2p1TightPass", &in->singleE25eta2p1TightPass);
  tree->SetBranchAddress("singleIsoMu22Pass", &in->singleIsoMu22Pass);
  tree->SetBranchAddress("singleIsoMu22eta2p1Pass", &in->singleIsoMu22eta2p1Pass);
  tree->SetBranchAddress("singleIsoTkMu22Pass", &in->singleIsoTkMu22Pass);
  tree->SetBranchAddress("singleIsoTkMu22eta2p1Pass", &in->singleIsoTkMu22eta2p1Pass);
  tree->SetBranchAddress("singleMu19eta2p1LooseTau20Pass", &in->singleMu19eta2p1LooseTau20Pass);
  tree->SetBranchAddress("singleMu19eta2p1LooseTau20singleL1Pass", &in->singleMu19eta2p1LooseTau20singleL1Pass);
  tree->SetBranchAddress("tAgainstElectronLooseMVA6", &in->tAgainstElectronLooseMVA6);
  tree->SetBranchAddress("tAgainstElectronLooseMVA62018", &in->tAgainstElectronLooseMVA62018);
  tree->SetBranchAddress("tAgainstElectronMVA6Raw", &in->tAgainstElectronMVA6Raw);
  tree->SetBranchAddress("tAgainstElectronMVA6Raw2018", &in->tAgainstElectronMVA6Raw2018);
  tree->SetBranchAddress("tAgainstElectronMVA6category", &in->tAgainstElectronMVA6category);
  tree->SetBranchAddress("tAgainstElectronMVA6category2018", &in->tAgainstElectronMVA6category2018);
  tree->SetBranchAddress("tAgainstElectronMediumMVA6", &in->tAgainstElectronMediumMVA6);
  tree->SetBranchAddress("tAgainstElectronMediumMVA62018", &in->tAgainstElectronMediumMVA62018);
  tree->SetBranchAddress("tAgainstElectronTightMVA6", &in->tAgainstElectronTightMVA6);
  tree->SetBranchAddress("tAgainstElectronTightMVA62018", &in->tAgainstElectronTightMVA62018);
  tree->SetBranchAddress("tAgainstElectronVLooseMVA6", &in->tAgainstElectronVLooseMVA6);
  tree->SetBranchAddress("tAgainstElectronVLooseMVA62018", &in->tAgainstElectronVLooseMVA62018);
  tree->SetBranchAddress("tAgainstElectronVTightMVA6", &in->tAgainstElectronVTightMVA6);
  tree->SetBranchAddress("tAgainstElectronVTightMVA62018", &in->tAgainstElectronVTightMVA62018);
  tree->SetBranchAddress("tAgainstMuonLoose3", &in->tAgainstMuonLoose3);
  tree->SetBranchAddress("tAgainstMuonTight3", &in->tAgainstMuonTight3);
  tree->SetBranchAddress("tByCombinedIsolationDeltaBetaCorrRaw3Hits", &in->tByCombinedIsolationDeltaBetaCorrRaw3Hits);
  tree->SetBranchAddress("tByIsolationMVArun2v1DBdR03oldDMwLTraw", &in->tByIsolationMVArun2v1DBdR03oldDMwLTraw);
  tree->SetBranchAddress("tByIsolationMVArun2v1DBnewDMwLTraw", &in->tByIsolationMVArun2v1DBnewDMwLTraw);
  tree->SetBranchAddress("tByIsolationMVArun2v1DBoldDMwLTraw", &in->tByIsolationMVArun2v1DBoldDMwLTraw);
  tree->SetBranchAddress("tByLooseCombinedIsolationDeltaBetaCorr3Hits", &in->tByLooseCombinedIsolationDeltaBetaCorr3Hits);
  tree->SetBranchAddress("tByLooseIsolationMVArun2v1DBdR03oldDMwLT", &in->tByLooseIsolationMVArun2v1DBdR03oldDMwLT);
  tree->SetBranchAddress("tByLooseIsolationMVArun2v1DBnewDMwLT", &in->tByLooseIsolationMVArun2v1DBnewDMwLT);
  tree->SetBranchAddress("tByLooseIsolationMVArun2v1DBoldDMwLT", &in->tByLooseIsolationMVArun2v1DBoldDMwLT);
  tree->SetBranchAddress("tByMediumCombinedIsolationDeltaBetaCorr3Hits", &in->tByMediumCombinedIsolationDeltaBetaCorr3Hits);
  tree->SetBranchAddress("tByMediumIsolationMVArun2v1DBdR03oldDMwLT", &in->tByMediumIsolationMVArun2v1DBdR03oldDMwLT);
  tree->SetBranchAddress("tByMediumIsolationMVArun2v1DBnewDMwLT", &in->tByMediumIsolationMVArun2v1DBnewDMwLT);
  tree->SetBranchAddress("tByMediumIsolationMVArun2v1DBoldDMwLT", &in->tByMediumIsolationMVArun2v1DBoldDMwLT);
  tree->SetBranchAddress("tByPhotonPtSumOutsideSignalCone", &in->tByPhotonPtSumOutsideSignalCone);
  tree->SetBranchAddress("tByTightCombinedIsolationDeltaBetaCorr3Hits", &in->tByTightCombinedIsolationDeltaBetaCorr3Hits);
  tree->SetBranchAddress("tByTightIsolationMVArun2v1DBdR03oldDMwLT", &in->tByTightIsolationMVArun2v1DBdR03oldDMwLT);
  tree->SetBranchAddress("tByTightIsolationMVArun2v1DBnewDMwLT", &in->tByTightIsolationMVArun2v1DBnewDMwLT);
  tree->SetBranchAddress("tByTightIsolationMVArun2v1DBoldDMwLT", &in->tByTightIsolationMVArun2v1DBoldDMwLT);
  tree->SetBranchAddress("tByVLooseIsolationMVArun2v1DBdR03oldDMwLT", &in->tByVLooseIsolationMVArun2v1DBdR03oldDMwLT);
  tree->SetBranchAddress("tByVLooseIsolationMVArun2v1DBnewDMwLT", &in->tByVLooseIsolationMVArun2v1DBnewDMwLT);
  tree->SetBranchAddress("tByVLooseIsolationMVArun2v1DBoldDMwLT", &in->tByVLooseIsolationMVArun2v1DBoldDMwLT);
  tree->SetBranchAddress("tByVTightIsolationMVArun2v1DBdR03oldDMwLT", &in->tByVTightIsolationMVArun2v1DBdR03oldDMwLT);
  tree->SetBranchAddress("tByVTightIsolationMVArun2v1DBnewDMwLT", &in->tByVTightIsolationMVArun2v1DBnewDMwLT);
  tree->SetBranchAddress("tByVTightIsolationMVArun2v1DBoldDMwLT", &in->tByVTightIsolationMVArun2v1DBoldDMwLT);
  tree->SetBranchAddress("tByVVTightIsolationMVArun2v1DBdR03oldDMwLT", &in->tByVVTightIsolationMVArun2v1DBdR03oldDMwLT);
  tree->SetBranchAddress("tByVVTightIsolationMVArun2v1DBnewDMwLT", &in->tByVVTightIsolationMVArun2v1DBnewDMwLT);
  tree->SetBranchAddress("tByVVTightIsolationMVArun2v1DBoldDMwLT", &in->tByVVTightIsolationMVArun2v1DBoldDMwLT);
  tree->SetBranchAddress("tCharge", &in->tCharge);
  tree->SetBranchAddress("tChargedIsoPtSum", &in->tChargedIsoPtSum);
  tree->SetBranchAddress("tChargedIsoPtSumdR03", &in->tChargedIsoPtSumdR03);
  tree->SetBranchAddress("tComesFromHiggs", &in->tComesFromHiggs);
  tree->SetBranchAddress("tDecayMode", &in->tDecayMode);
  tree->SetBranchAddress("tDecayModeFinding", &in->tDecayModeFinding);
  tree->SetBranchAddress("tDecayModeFindingNewDMs", &in->tDecayModeFindingNewDMs);
  tree->SetBranchAddress("tDeepTau2017v1VSeraw", &in->tDeepTau2017v1VSeraw);
  tree->SetBranchAddress("tDeepTau2017v1VSjetraw", &in->tDeepTau2017v1VSjetraw);
  tree->SetBranchAddress("tDeepTau2017v1VSmuraw", &in->tDeepTau2017v1VSmuraw);
  tree->SetBranchAddress("tDpfTau2016v0VSallraw", &in->tDpfTau2016v0VSallraw);
  tree->SetBranchAddress("tDpfTau2016v1VSallraw", &in->tDpfTau2016v1VSallraw);
  tree->SetBranchAddress("tEta", &in->tEta);
  tree->SetBranchAddress("tFootprintCorrection", &in->tFootprintCorrection);
  tree->SetBranchAddress("tFootprintCorrectiondR03", &in->tFootprintCorrectiondR03);
  tree->SetBranchAddress("tGenCharge", &in->tGenCharge);
  tree->SetBranchAddress("tGenDecayMode", &in->tGenDecayMode);
  tree->SetBranchAddress("tGenEnergy", &in->tGenEnergy);
  tree->SetBranchAddress("tGenEta", &in->tGenEta);
  tree->SetBranchAddress("tGenJetEta", &in->tGenJetEta);
  tree->SetBranchAddress("tGenJetPt", &in->tGenJetPt);
  tree->SetBranchAddress("tGenMotherEnergy", &in->tGenMotherEnergy);
  tree->SetBranchAddress("tGenMotherEta", &in->tGenMotherEta);
  tree->SetBranchAddress("tGenMotherPdgId", &in->tGenMotherPdgId);
  tree->SetBranchAddress("tGenMotherPhi", &in->tGenMotherPhi);
  tree->SetBranchAddress("tGenMotherPt", &in->tGenMotherPt);
  tree->SetBranchAddress("tGenPdgId", &in->tGenPdgId);
  tree->SetBranchAddress("tGenPhi", &in->tGenPhi);
  tree->SetBranchAddress("tGenPt", &in->tGenPt);
  tree->SetBranchAddress("tGenStatus", &in->tGenStatus);
  tree->SetBranchAddress("tJetArea", &in->tJetArea);
  tree->SetBranchAddress("tJetBtag", &in->tJetBtag);
  tree->SetBranchAddress("tJetDR", &in->tJetDR);
  tree->SetBranchAddress("tJetEtaEtaMoment", &in->tJetEtaEtaMoment);
  tree->SetBranchAddress("tJetEtaPhiMoment", &in->tJetEtaPhiMoment);
  tree->SetBranchAddress("tJetEtaPhiSpread", &in->tJetEtaPhiSpread);
  tree->SetBranchAddress("tJetHadronFlavour", &in->tJetHadronFlavour);
  tree->SetBranchAddress("tJetPFCISVBtag", &in->tJetPFCISVBtag);
  tree->SetBranchAddress("tJetPartonFlavour", &in->tJetPartonFlavour);
  tree->SetBranchAddress("tJetPhiPhiMoment", &in->tJetPhiPhiMoment);
  tree->SetBranchAddress("tJetPt", &in->tJetPt);
  tree->SetBranchAddress("tL1IsoTauMatch", &in->tL1IsoTauMatch);
  tree->SetBranchAddress("tL1IsoTauPt", &in->tL1IsoTauPt);
  tree->SetBranchAddress("tLeadTrackPt", &in->tLeadTrackPt);
  tree->SetBranchAddress("tLooseDeepTau2017v1VSe", &in->tLooseDeepTau2017v1VSe);
  tree->SetBranchAddress("tLooseDeepTau2017v1VSjet", &in->tLooseDeepTau2017v1VSjet);
  tree->SetBranchAddress("tLooseDeepTau2017v1VSmu", &in->tLooseDeepTau2017v1VSmu);
  tree->SetBranchAddress("tLowestMll", &in->tLowestMll);
  tree->SetBranchAddress("tMass", &in->tMass);
  tree->SetBranchAddress("tMatchesDoubleMediumCombinedIsoTau35Path", &in->tMatchesDoubleMediumCombinedIsoTau35Path);
  tree->SetBranchAddress("tMatchesDoubleMediumHPSTau35Filter", &in->tMatchesDoubleMediumHPSTau35Filter);
  tree->SetBranchAddress("tMatchesDoubleMediumHPSTau35Path", &in->tMatchesDoubleMediumHPSTau35Path);
  tree->SetBranchAddress("tMatchesDoubleMediumHPSTau40Filter", &in->tMatchesDoubleMediumHPSTau40Filter);
  tree->SetBranchAddress("tMatchesDoubleMediumHPSTau40Path", &in->tMatchesDoubleMediumHPSTau40Path);
  tree->SetBranchAddress("tMatchesDoubleMediumIsoTau35Path", &in->tMatchesDoubleMediumIsoTau35Path);
  tree->SetBranchAddress("tMatchesDoubleMediumTau35Filter", &in->tMatchesDoubleMediumTau35Filter);
  tree->SetBranchAddress("tMatchesDoubleMediumTau35Path", &in->tMatchesDoubleMediumTau35Path);
  tree->SetBranchAddress("tMatchesDoubleMediumTau40Filter", &in->tMatchesDoubleMediumTau40Filter);
  tree->SetBranchAddress("tMatchesDoubleMediumTau40Path", &in->tMatchesDoubleMediumTau40Path);
  tree->SetBranchAddress("tMatchesDoubleTightHPSTau35Filter", &in->tMatchesDoubleTightHPSTau35Filter);
  tree->SetBranchAddress("tMatchesDoubleTightHPSTau35Path", &in->tMatchesDoubleTightHPSTau35Path);
  tree->SetBranchAddress("tMatchesDoubleTightHPSTau40Filter", &in->tMatchesDoubleTightHPSTau40Filter);
  tree->SetBranchAddress("tMatchesDoubleTightHPSTau40Path", &in->tMatchesDoubleTightHPSTau40Path);
  tree->SetBranchAddress("tMatchesDoubleTightTau35Filter", &in->tMatchesDoubleTightTau35Filter);
  tree->SetBranchAddress("tMatchesDoubleTightTau35Path", &in->tMatchesDoubleTightTau35Path);
  tree->SetBranchAddress("tMatchesDoubleTightTau40Filter", &in->tMatchesDoubleTightTau40Filter);
  tree->SetBranchAddress("tMatchesDoubleTightTau40Path", &in->tMatchesDoubleTightTau40Path);
  tree->SetBranchAddress("tMatchesEle24HPSTau30Filter", &in->tMatchesEle24HPSTau30Filter);
  tree->SetBranchAddress("tMatchesEle24HPSTau30Path", &in->tMatchesEle24HPSTau30Path);
  tree->SetBranchAddress("tMatchesEle24Tau30Filter", &in->tMatchesEle24Tau30Filter);
  tree->SetBranchAddress("tMatchesEle24Tau30Path", &in->tMatchesEle24Tau30Path);
  tree->SetBranchAddress("tMatchesIsoMu19Tau20Filter", &in->tMatchesIsoMu19Tau20Filter);
  tree->SetBranchAddress("tMatchesIsoMu19Tau20Path", &in->tMatchesIsoMu19Tau20Path);
  tree->SetBranchAddress("tMatchesIsoMu19Tau20SingleL1Filter", &in->tMatchesIsoMu19Tau20SingleL1Filter);
  tree->SetBranchAddress("tMatchesIsoMu19Tau20SingleL1Path", &in->tMatchesIsoMu19Tau20SingleL1Path);
  tree->SetBranchAddress("tMatchesIsoMu20HPSTau27Filter", &in->tMatchesIsoMu20HPSTau27Filter);
  tree->SetBranchAddress("tMatchesIsoMu20HPSTau27Path", &in->tMatchesIsoMu20HPSTau27Path);
  tree->SetBranchAddress("tMatchesIsoMu20Tau27Filter", &in->tMatchesIsoMu20Tau27Filter);
  tree->SetBranchAddress("tMatchesIsoMu20Tau27Path", &in->tMatchesIsoMu20Tau27Path);
  tree->SetBranchAddress("tMediumDeepTau2017v1VSe", &in->tMediumDeepTau2017v1VSe);
  tree->SetBranchAddress("tMediumDeepTau2017v1VSjet", &in->tMediumDeepTau2017v1VSjet);
  tree->SetBranchAddress("tMediumDeepTau2017v1VSmu", &in->tMediumDeepTau2017v1VSmu);
  tree->SetBranchAddress("tNChrgHadrIsolationCands", &in->tNChrgHadrIsolationCands);
  tree->SetBranchAddress("tNChrgHadrSignalCands", &in->tNChrgHadrSignalCands);
  tree->SetBranchAddress("tNGammaSignalCands", &in->tNGammaSignalCands);
  tree->SetBranchAddress("tNNeutralHadrSignalCands", &in->tNNeutralHadrSignalCands);
  tree->SetBranchAddress("tNSignalCands", &in->tNSignalCands);
  tree->SetBranchAddress("tNearestZMass", &in->tNearestZMass);
  tree->SetBranchAddress("tNeutralIsoPtSum", &in->tNeutralIsoPtSum);
  tree->SetBranchAddress("tNeutralIsoPtSumWeight", &in->tNeutralIsoPtSumWeight);
  tree->SetBranchAddress("tNeutralIsoPtSumWeightdR03", &in->tNeutralIsoPtSumWeightdR03);
  tree->SetBranchAddress("tNeutralIsoPtSumdR03", &in->tNeutralIsoPtSumdR03);
  tree->SetBranchAddress("tPVDXY", &in->tPVDXY);
  tree->SetBranchAddress("tPVDZ", &in->tPVDZ);
  tree->SetBranchAddress("tPhi", &in->tPhi);
  tree->SetBranchAddress("tPhotonPtSumOutsideSignalCone", &in->tPhotonPtSumOutsideSignalCone);
  tree->SetBranchAddress("tPhotonPtSumOutsideSignalConedR03", &in->tPhotonPtSumOutsideSignalConedR03);
  tree->SetBranchAddress("tPt", &in->tPt);
  tree->SetBranchAddress("tPuCorrPtSum", &in->tPuCorrPtSum);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTLoose", &in->tRerunMVArun2v2DBoldDMwLTLoose);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTMedium", &in->tRerunMVArun2v2DBoldDMwLTMedium);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTTight", &in->tRerunMVArun2v2DBoldDMwLTTight);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVLoose", &in->tRerunMVArun2v2DBoldDMwLTVLoose);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVTight", &in->tRerunMVArun2v2DBoldDMwLTVTight);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVVLoose", &in->tRerunMVArun2v2DBoldDMwLTVVLoose);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVVTight", &in->tRerunMVArun2v2DBoldDMwLTVVTight);
  tree->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTraw", &in->tRerunMVArun2v2DBoldDMwLTraw);
  tree->SetBranchAddress("tTightDeepTau2017v1VSe", &in->tTightDeepTau2017v1VSe);
  tree->SetBranchAddress("tTightDeepTau2017v1VSjet", &in->tTightDeepTau2017v1VSjet);
  tree->SetBranchAddress("tTightDeepTau2017v1VSmu", &in->tTightDeepTau2017v1VSmu);
  tree->SetBranchAddress("tTightDpfTau2016v0VSall", &in->tTightDpfTau2016v0VSall);
  tree->SetBranchAddress("tTightDpfTau2016v1VSall", &in->tTightDpfTau2016v1VSall);
  tree->SetBranchAddress("tVLooseDeepTau2017v1VSe", &in->tVLooseDeepTau2017v1VSe);
  tree->SetBranchAddress("tVLooseDeepTau2017v1VSjet", &in->tVLooseDeepTau2017v1VSjet);
  tree->SetBranchAddress("tVLooseDeepTau2017v1VSmu", &in->tVLooseDeepTau2017v1VSmu);
  tree->SetBranchAddress("tVTightDeepTau2017v1VSe", &in->tVTightDeepTau2017v1VSe);
  tree->SetBranchAddress("tVTightDeepTau2017v1VSjet", &in->tVTightDeepTau2017v1VSjet);
  tree->SetBranchAddress("tVTightDeepTau2017v1VSmu", &in->tVTightDeepTau2017v1VSmu);
  tree->SetBranchAddress("tVVLooseDeepTau2017v1VSe", &in->tVVLooseDeepTau2017v1VSe);
  tree->SetBranchAddress("tVVLooseDeepTau2017v1VSjet", &in->tVVLooseDeepTau2017v1VSjet);
  tree->SetBranchAddress("tVVLooseDeepTau2017v1VSmu", &in->tVVLooseDeepTau2017v1VSmu);
  tree->SetBranchAddress("tVVTightDeepTau2017v1VSe", &in->tVVTightDeepTau2017v1VSe);
  tree->SetBranchAddress("tVVTightDeepTau2017v1VSjet", &in->tVVTightDeepTau2017v1VSjet);
  tree->SetBranchAddress("tVVTightDeepTau2017v1VSmu", &in->tVVTightDeepTau2017v1VSmu);
  tree->SetBranchAddress("tVVVLooseDeepTau2017v1VSe", &in->tVVVLooseDeepTau2017v1VSe);
  tree->SetBranchAddress("tVVVLooseDeepTau2017v1VSjet", &in->tVVVLooseDeepTau2017v1VSjet);
  tree->SetBranchAddress("tVVVLooseDeepTau2017v1VSmu", &in->tVVVLooseDeepTau2017v1VSmu);
  tree->SetBranchAddress("tVZ", &in->tVZ);
  tree->SetBranchAddress("tZTTGenDR", &in->tZTTGenDR);
  tree->SetBranchAddress("tZTTGenEta", &in->tZTTGenEta);
  tree->SetBranchAddress("tZTTGenMatching", &in->tZTTGenMatching);
  tree->SetBranchAddress("tZTTGenPhi", &in->tZTTGenPhi);
  tree->SetBranchAddress("tZTTGenPt", &in->tZTTGenPt);
  tree->SetBranchAddress("tauVetoPt20Loose3HitsVtx", &in->tauVetoPt20Loose3HitsVtx);
  tree->SetBranchAddress("tauVetoPt20TightMVALTVtx", &in->tauVetoPt20TightMVALTVtx);
  tree->SetBranchAddress("topQuarkPt1", &in->topQuarkPt1);
  tree->SetBranchAddress("topQuarkPt2", &in->topQuarkPt2);
  tree->SetBranchAddress("tripleEPass", &in->tripleEPass);
  tree->SetBranchAddress("tripleMu10_5_5Pass", &in->tripleMu10_5_5Pass);
  tree->SetBranchAddress("tripleMu12_10_5Pass", &in->tripleMu12_10_5Pass);
  tree->SetBranchAddress("type1_pfMetEt", &in->type1_pfMetEt);
  tree->SetBranchAddress("type1_pfMetPhi", &in->type1_pfMetPhi);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnDown", &in->type1_pfMet_shiftedPhi_JetEnDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnUp", &in->type1_pfMet_shiftedPhi_JetEnUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEta0to3Down", &in->type1_pfMet_shiftedPhi_JetEta0to3Down);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEta0to3Up", &in->type1_pfMet_shiftedPhi_JetEta0to3Up);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEta0to5Down", &in->type1_pfMet_shiftedPhi_JetEta0to5Down);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEta0to5Up", &in->type1_pfMet_shiftedPhi_JetEta0to5Up);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEta3to5Down", &in->type1_pfMet_shiftedPhi_JetEta3to5Down);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetEta3to5Up", &in->type1_pfMet_shiftedPhi_JetEta3to5Up);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetRelativeBalDown", &in->type1_pfMet_shiftedPhi_JetRelativeBalDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetRelativeBalUp", &in->type1_pfMet_shiftedPhi_JetRelativeBalUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetRelativeSampleDown", &in->type1_pfMet_shiftedPhi_JetRelativeSampleDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetRelativeSampleUp", &in->type1_pfMet_shiftedPhi_JetRelativeSampleUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetResDown", &in->type1_pfMet_shiftedPhi_JetResDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetResUp", &in->type1_pfMet_shiftedPhi_JetResUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetTotalDown", &in->type1_pfMet_shiftedPhi_JetTotalDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_JetTotalUp", &in->type1_pfMet_shiftedPhi_JetTotalUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnDown", &in->type1_pfMet_shiftedPhi_UnclusteredEnDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnUp", &in->type1_pfMet_shiftedPhi_UnclusteredEnUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEnDown", &in->type1_pfMet_shiftedPt_JetEnDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEnUp", &in->type1_pfMet_shiftedPt_JetEnUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEta0to3Down", &in->type1_pfMet_shiftedPt_JetEta0to3Down);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEta0to3Up", &in->type1_pfMet_shiftedPt_JetEta0to3Up);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEta0to5Down", &in->type1_pfMet_shiftedPt_JetEta0to5Down);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEta0to5Up", &in->type1_pfMet_shiftedPt_JetEta0to5Up);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEta3to5Down", &in->type1_pfMet_shiftedPt_JetEta3to5Down);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetEta3to5Up", &in->type1_pfMet_shiftedPt_JetEta3to5Up);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetRelativeBalDown", &in->type1_pfMet_shiftedPt_JetRelativeBalDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetRelativeBalUp", &in->type1_pfMet_shiftedPt_JetRelativeBalUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetRelativeSampleDown", &in->type1_pfMet_shiftedPt_JetRelativeSampleDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetRelativeSampleUp", &in->type1_pfMet_shiftedPt_JetRelativeSampleUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetResDown", &in->type1_pfMet_shiftedPt_JetResDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetResUp", &in->type1_pfMet_shiftedPt_JetResUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetTotalDown", &in->type1_pfMet_shiftedPt_JetTotalDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_JetTotalUp", &in->type1_pfMet_shiftedPt_JetTotalUp);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnDown", &in->type1_pfMet_shiftedPt_UnclusteredEnDown);
  tree->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnUp", &in->type1_pfMet_shiftedPt_UnclusteredEnUp);
  tree->SetBranchAddress("vbfDeta", &in->vbfDeta);
  tree->SetBranchAddress("vbfJetVeto20", &in->vbfJetVeto20);
  tree->SetBranchAddress("vbfJetVeto30", &in->vbfJetVeto30);
  tree->SetBranchAddress("vbfMass", &in->vbfMass);
  tree->SetBranchAddress("vbfMassWoNoisyJets", &in->vbfMassWoNoisyJets);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEC2Down", &in->vbfMassWoNoisyJets_JetEC2Down);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEC2Up", &in->vbfMassWoNoisyJets_JetEC2Up);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEta0to3Down", &in->vbfMassWoNoisyJets_JetEta0to3Down);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEta0to3Up", &in->vbfMassWoNoisyJets_JetEta0to3Up);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEta0to5Down", &in->vbfMassWoNoisyJets_JetEta0to5Down);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEta0to5Up", &in->vbfMassWoNoisyJets_JetEta0to5Up);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEta3to5Down", &in->vbfMassWoNoisyJets_JetEta3to5Down);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetEta3to5Up", &in->vbfMassWoNoisyJets_JetEta3to5Up);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetRelativeBalDown", &in->vbfMassWoNoisyJets_JetRelativeBalDown);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetRelativeBalUp", &in->vbfMassWoNoisyJets_JetRelativeBalUp);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetRelativeSampleDown", &in->vbfMassWoNoisyJets_JetRelativeSampleDown);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetRelativeSampleUp", &in->vbfMassWoNoisyJets_JetRelativeSampleUp);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetTotalDown", &in->vbfMassWoNoisyJets_JetTotalDown);
  tree->SetBranchAddress("vbfMassWoNoisyJets_JetTotalUp", &in->vbfMassWoNoisyJets_JetTotalUp);
  tree->SetBranchAddress("vbfMass_JetEC2Down", &in->vbfMass_JetEC2Down);
  tree->SetBranchAddress("vbfMass_JetEC2Up", &in->vbfMass_JetEC2Up);
  tree->SetBranchAddress("vbfMass_JetEta0to3Down", &in->vbfMass_JetEta0to3Down);
  tree->SetBranchAddress("vbfMass_JetEta0to3Up", &in->vbfMass_JetEta0to3Up);
  tree->SetBranchAddress("vbfMass_JetEta0to5Down", &in->vbfMass_JetEta0to5Down);
  tree->SetBranchAddress("vbfMass_JetEta0to5Up", &in->vbfMass_JetEta0to5Up);
  tree->SetBranchAddress("vbfMass_JetEta3to5Down", &in->vbfMass_JetEta3to5Down);
  tree->SetBranchAddress("vbfMass_JetEta3to5Up", &in->vbfMass_JetEta3to5Up);
  tree->SetBranchAddress("vbfMass_JetRelativeBalDown", &in->vbfMass_JetRelativeBalDown);
  tree->SetBranchAddress("vbfMass_JetRelativeBalUp", &in->vbfMass_JetRelativeBalUp);
  tree->SetBranchAddress("vbfMass_JetRelativeSampleDown", &in->vbfMass_JetRelativeSampleDown);
  tree->SetBranchAddress("vbfMass_JetRelativeSampleUp", &in->vbfMass_JetRelativeSampleUp);
  tree->SetBranchAddress("vbfMass_JetTotalDown", &in->vbfMass_JetTotalDown);
  tree->SetBranchAddress("vbfMass_JetTotalUp", &in->vbfMass_JetTotalUp);
  tree->SetBranchAddress("vbfNJets20", &in->vbfNJets20);
  tree->SetBranchAddress("vbfNJets30", &in->vbfNJets30);
  tree->SetBranchAddress("vbfj1eta", &in->vbfj1eta);
  tree->SetBranchAddress("vbfj1pt", &in->vbfj1pt);
  tree->SetBranchAddress("vbfj2eta", &in->vbfj2eta);
  tree->SetBranchAddress("vbfj2pt", &in->vbfj2pt);
  tree->SetBranchAddress("vispX", &in->vispX);
  tree->SetBranchAddress("vispY", &in->vispY);
  tree->SetBranchAddress("idx", &in->idx);
}

#endif  // ROOT_SRC_ETAU_TREE_2017_H_
