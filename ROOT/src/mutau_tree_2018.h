// Copyright 2019 Tyler Mitchell

#ifndef ROOT_SRC_MUTAU_TREE_2018_H_
#define ROOT_SRC_MUTAU_TREE_2018_H_

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "../interface/mutau_input_branches.h"
#include "./base_tree.h"
#include "RecoilCorrector.h"
#include "TLorentzVector.h"
#include "TTree.h"

class mutau_tree2018 : public virtual base_tree {
 private:
  TTree *tree, *original;
  mutau_input_branches* in;
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
  mutau_tree2018(TTree* orig, TTree* itree, bool isMC, bool isEmbed, Int_t rec);
  virtual ~mutau_tree2018() {}
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
mutau_tree2018::mutau_tree2018(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec) : tree(itree),
                                                                                                    original(Original),
                                                                                                    in(new mutau_input_branches(Original)),
                                                                                                    isMC(IsMC),
                                                                                                    isEmbed(IsEmbed),
                                                                                                    recoil(rec) {}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void mutau_tree2018::do_skimming(TH1F* cutflow) {
  // declare variables for sorting
  ULong64_t evt_now(0);
  ULong64_t evt_before(1);
  int best_evt(-1);
  std::pair<float, float> muCandidate, tauCandidate;

  std::cout << "Starting the skim..." << std::endl;

  bool isData = !isEmbed && !isMC;
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
          tau *= 0.987;
        } else if (in->tDecayMode == 1) {
          tau *= 0.995;
        } else if (in->tDecayMode == 10) {
          tau *= 0.988;
        }
      } else if (in->tZTTGenMatching == 1 || in->tZTTGenMatching == 3) {
        if (in->tDecayMode == 0) {
          tau *= 0.968;
        } else if (in->tDecayMode == 1) {
          tau *= 1.026;
        }
      } else if (in->tZTTGenMatching == 2 || in->tZTTGenMatching == 4) {
        if (in->tDecayMode == 0) {
          tau *= 0.998;
        } else if (in->tDecayMode == 1) {
          tau *= 0.990;
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

    cutflow->Fill(1., 1.);
    // apply event selection

    auto Mu24 = in->IsoMu24Pass && in->mMatchesIsoMu24Path && in->mMatchesIsoMu24Filter;
    auto Mu27 = in->IsoMu27Pass && in->mMatchesIsoMu27Path && in->mMatchesIsoMu27Filter;
    auto Cross_base = in->mMatchesIsoMu20HPSTau27Filter && in->mMatchesIsoMu20HPSTau27Path
                      && in->tMatchesIsoMu20HPSTau27Filter && in->tMatchesIsoMu20HPSTau27Path;
    auto Cross_v1 = Cross_base && in->Mu20LooseHPSTau27Pass && in->run < 317509 && isData;
    auto Cross_v2 = Cross_base && in->Mu20LooseHPSTau27TightIDPass && (isMC || (in->run > 317509 && isData));
    auto Mu24_emb = in->mMatchEmbeddedFilterMu24;
    auto Mu27_emb = in->mMatchEmbeddedFilterMu27;
    auto Cross_emb = in->mMatchEmbeddedFilterMu20Tau27_2018 && in->tMatchEmbeddedFilterMu20HPSTau27;

    if (isEmbed) {
      if (Mu27_emb && in->mPt > 28) {
        cutflow->Fill(2., 1.);
      } else if (Mu24_emb && in->mPt > 25) {
        cutflow->Fill(2., 1.);
      } else if (Cross_emb && in->mPt > 21 && in->mPt < 25 && tau.Pt() > 31 && fabs(in->mEta) < 2.1 && fabs(tau.Eta()) < 2.1) {
        cutflow->Fill(2., 1.);
      } else {
        continue;
      }
    } else {
      if (Mu27 && in->mPt > 28) {
        cutflow->Fill(2., 1.);
      } else if (Mu24 && in->mPt > 25) {
        cutflow->Fill(2., 1.);
      } else if ((Cross_v1 || Cross_v2) && in->mPt > 21 && in->mPt < 25 && tau.Pt() > 31 && fabs(in->mEta) < 2.1 && fabs(tau.Eta()) < 2.1) {
        cutflow->Fill(2., 1.);
      } else {
        continue;
      }
    }

    if (in->mPt > 21. && fabs(in->mEta) < 2.4 && fabs(in->mPVDZ) < 0.2 && fabs(in->mPVDXY) < 0.045)
      cutflow->Fill(3., 1.);  // electron kinematic selection
    else
      continue;

    if (in->mPFIDMedium)
      cutflow->Fill(4., 1.);  // muon quality selection
    else
      continue;

    if (tau.Pt() > 30. && fabs(tau.Eta()) < 2.3 && fabs(in->tPVDZ) < 0.2)
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

    if (in->m_t_DR > 0.5) {
      cutflow->Fill(10., 1.);
    } else {
      continue;
    }

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

Float_t mutau_tree2018::do_tes_met_corr(Float_t decayMode, Float_t sf1, Float_t sf2, Float_t sf3, TLorentzVector& met, TLorentzVector tau) {
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
TTree* mutau_tree2018::fill_tree(RecoilCorrector recoilPFMetCorrector) {
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
    njets = in->jetVeto30;
    nbtag = in->bjetDeepCSVVeto20Medium_2018_DR0p5;
    njetspt20 = in->jetVeto20;

    // TLorentzVector ele, tau;
    mu.SetPtEtaPhiM(in->mPt, in->mEta, in->mPhi, in->mMass);
    tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);

    met_px = in->type1_pfMetEt * cos(in->type1_pfMetPhi);
    met_py = in->type1_pfMetEt * sin(in->type1_pfMetPhi);

    extraelec_veto = in->eVetoZTTp001dxyzR0 > 0;
    extramuon_veto = in->muVetoZTTp001dxyzR0 > 1;
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
          jetVeto30 + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,          // number of jets (hadronic jet multiplicity) (int)
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
          jetVeto30,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          jetVeto30,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,            // generator Z/W/Higgs px (float)
          in->genpY,            // generator Z/W/Higgs py (float)
          in->vispX,            // generator visible Z/W/Higgs px (float)
          in->vispY,            // generator visible Z/W/Higgs py (float)
          jetVeto30,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          jetVeto30,              // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,              // generator Z/W/Higgs px (float)
          in->genpY,              // generator Z/W/Higgs py (float)
          in->vispX,              // generator visible Z/W/Higgs px (float)
          in->vispY,              // generator visible Z/W/Higgs py (float)
          jetVeto30,              // number of jets (hadronic jet multiplicity) (int)
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
        auto sf = do_tes_met_corr(in->tDecayMode, 0.987, 0.995, 0.988, MET, tau);
        do_tes_met_corr(in->tDecayMode, 0.987, 0.995, 0.988, MET_JESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.987, 0.995, 0.988, MET_JESDown, tau);
        do_tes_met_corr(in->tDecayMode, 0.987, 0.995, 0.988, MET_UESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.987, 0.995, 0.988, MET_UESDown, tau);
        tau *= sf;
      } else if (in->tZTTGenMatching == 1 || in->tZTTGenMatching == 3) {
        auto sf = do_tes_met_corr(in->tDecayMode, 0.968, 1.026, 1.00, MET, tau);
        do_tes_met_corr(in->tDecayMode, 0.968, 1.026, 1.000, MET_JESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.968, 1.026, 1.000, MET_JESDown, tau);
        do_tes_met_corr(in->tDecayMode, 0.968, 1.026, 1.000, MET_UESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.968, 1.026, 1.000, MET_UESDown, tau);
        tau *= sf;
      } else if (in->tZTTGenMatching == 2 || in->tZTTGenMatching == 4) {
        auto sf = do_tes_met_corr(in->tDecayMode, 0.998, 0.990, 1.00, MET, tau);
        do_tes_met_corr(in->tDecayMode, 0.998, 0.990, 1.000, MET_JESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.998, 0.990, 1.000, MET_JESDown, tau);
        do_tes_met_corr(in->tDecayMode, 0.998, 0.990, 1.000, MET_UESUp, tau);
        do_tes_met_corr(in->tDecayMode, 0.998, 0.990, 1.000, MET_UESDown, tau);
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
void mutau_tree2018::set_branches() {
  // new branches
  tree->Branch("run", &Run);
  tree->Branch("lumi", &Lumi);
  tree->Branch("gen_match_1", &gen_match_1);
  tree->Branch("gen_match_2", &gen_match_2);
  tree->Branch("met_px", &met_px);
  tree->Branch("met_py", &met_py);
  tree->Branch("extraelec_veto", &extraelec_veto);
  tree->Branch("extramuon_veto", &extramuon_veto);
  tree->Branch("dilepton_veto", &dilepton_veto);
  tree->Branch("pfmetcorr_ex", &pfmetcorr_ex);
  tree->Branch("pfmetcorr_ey", &pfmetcorr_ey);
  tree->Branch("pfmetcorr_ex_UESUp", &pfmetcorr_ex_UESUp);
  tree->Branch("pfmetcorr_ey_UESUp", &pfmetcorr_ey_UESUp);
  tree->Branch("pfmetcorr_ex_UESDown", &pfmetcorr_ex_UESDown);
  tree->Branch("pfmetcorr_ey_UESDown", &pfmetcorr_ey_UESDown);
  tree->Branch("pfmetcorr_ex_JESUp", &pfmetcorr_ex_JESUp);
  tree->Branch("pfmetcorr_ey_JESUp", &pfmetcorr_ey_JESUp);
  tree->Branch("pfmetcorr_ex_JESDown", &pfmetcorr_ex_JESDown);
  tree->Branch("pfmetcorr_ey_JESDown", &pfmetcorr_ey_JESDown);
  tree->Branch("met", &met);
  tree->Branch("metphi", &metphi);
  tree->Branch("met_px", &met_px);
  tree->Branch("met_py", &met_py);
  tree->Branch("met_JESUp", &met_JESUp);
  tree->Branch("met_JESDown", &met_JESDown);
  tree->Branch("met_UESUp", &met_UESUp);
  tree->Branch("met_UESDown", &met_UESDown);
  tree->Branch("metphi_JESUp", &metphi_JESUp);
  tree->Branch("metphi_JESDown", &metphi_JESDown);
  tree->Branch("metphi_UESUp", &metphi_UESUp);
  tree->Branch("metphi_UESDown", &metphi_UESDown);
  tree->Branch("m_1", &m_1);
  tree->Branch("px_1", &px_1);
  tree->Branch("py_1", &py_1);
  tree->Branch("pz_1", &pz_1);
  tree->Branch("e_1", &e_1);
  tree->Branch("pt_1", &pt_1);
  tree->Branch("phi_1", &phi_1);
  tree->Branch("eta_1", &eta_1);
  tree->Branch("m_2", &m_2);
  tree->Branch("px_2", &px_2);
  tree->Branch("py_2", &py_2);
  tree->Branch("pz_2", &pz_2);
  tree->Branch("e_2", &e_2);
  tree->Branch("pt_2", &pt_2);
  tree->Branch("phi_2", &phi_2);
  tree->Branch("eta_2", &eta_2);

  // SVFit and MELA branches
  tree->Branch("q_1", &in->mCharge);
  tree->Branch("q_2", &in->tCharge);
  tree->Branch("jpt_1", &in->j1pt);
  tree->Branch("jeta_1", &in->j1eta);
  tree->Branch("jphi_1", &in->j1phi);
  tree->Branch("jpt_2", &in->j2pt);
  tree->Branch("jeta_2", &in->j2eta);
  tree->Branch("jphi_2", &in->j2phi);
  tree->Branch("l2_decayMode", &in->tDecayMode);
  tree->Branch("njets", &njets);
  tree->Branch("nbtag", &nbtag);
  tree->Branch("njetspt20", &njetspt20);

  // copy the rest
  tree->Branch("EmbPtWeight", &in->EmbPtWeight);
  tree->Branch("Eta", &in->Eta);
  tree->Branch("Flag_BadChargedCandidateFilter", &in->Flag_BadChargedCandidateFilter);
  tree->Branch("Flag_BadPFMuonFilter", &in->Flag_BadPFMuonFilter);
  tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &in->Flag_EcalDeadCellTriggerPrimitiveFilter);
  tree->Branch("Flag_HBHENoiseFilter", &in->Flag_HBHENoiseFilter);
  tree->Branch("Flag_HBHENoiseIsoFilter", &in->Flag_HBHENoiseIsoFilter);
  tree->Branch("Flag_badMuons", &in->Flag_badMuons);
  tree->Branch("Flag_duplicateMuons", &in->Flag_duplicateMuons);
  tree->Branch("Flag_ecalBadCalibFilter", &in->Flag_ecalBadCalibFilter);
  tree->Branch("Flag_eeBadScFilter", &in->Flag_eeBadScFilter);
  tree->Branch("Flag_globalSuperTightHalo2016Filter", &in->Flag_globalSuperTightHalo2016Filter);
  tree->Branch("Flag_globalTightHalo2016Filter", &in->Flag_globalTightHalo2016Filter);
  tree->Branch("Flag_goodVertices", &in->Flag_goodVertices);
  tree->Branch("GenWeight", &in->GenWeight);
  tree->Branch("Ht", &in->Ht);
  tree->Branch("IsoMu24Pass", &in->IsoMu24Pass);
  tree->Branch("IsoMu27Pass", &in->IsoMu27Pass);
  tree->Branch("LT", &in->LT);
  tree->Branch("Mass", &in->Mass);
  tree->Branch("MassError", &in->MassError);
  tree->Branch("MassErrord1", &in->MassErrord1);
  tree->Branch("MassErrord2", &in->MassErrord2);
  tree->Branch("MassErrord3", &in->MassErrord3);
  tree->Branch("MassErrord4", &in->MassErrord4);
  tree->Branch("Mt", &in->Mt);
  tree->Branch("Mu20LooseHPSTau27Pass", &in->Mu20LooseHPSTau27Pass);
  tree->Branch("Mu20LooseHPSTau27TightIDPass", &in->Mu20LooseHPSTau27TightIDPass);
  tree->Branch("Mu20LooseTau27Pass", &in->Mu20LooseTau27Pass);
  tree->Branch("Mu20LooseTau27TightIDPass", &in->Mu20LooseTau27TightIDPass);
  tree->Branch("Mu50Pass", &in->Mu50Pass);
  tree->Branch("NUP", &in->NUP);
  tree->Branch("Phi", &in->Phi);
  tree->Branch("Pt", &in->Pt);
  tree->Branch("Rivet_VEta", &in->Rivet_VEta);
  tree->Branch("Rivet_VPt", &in->Rivet_VPt);
  tree->Branch("Rivet_errorCode", &in->Rivet_errorCode);
  tree->Branch("Rivet_higgsEta", &in->Rivet_higgsEta);
  tree->Branch("Rivet_higgsPt", &in->Rivet_higgsPt);
  tree->Branch("Rivet_nJets25", &in->Rivet_nJets25);
  tree->Branch("Rivet_nJets30", &in->Rivet_nJets30);
  tree->Branch("Rivet_p4decay_VEta", &in->Rivet_p4decay_VEta);
  tree->Branch("Rivet_p4decay_VPt", &in->Rivet_p4decay_VPt);
  tree->Branch("Rivet_prodMode", &in->Rivet_prodMode);
  tree->Branch("Rivet_stage0_cat", &in->Rivet_stage0_cat);
  tree->Branch("Rivet_stage1_cat_pTjet25GeV", &in->Rivet_stage1_cat_pTjet25GeV);
  tree->Branch("Rivet_stage1_cat_pTjet30GeV", &in->Rivet_stage1_cat_pTjet30GeV);
  tree->Branch("Rivet_stage1p1_cat", &in->Rivet_stage1p1_cat);
  tree->Branch("bjetDeepCSVVeto20Loose_2016_DR0p5", &in->bjetDeepCSVVeto20Loose_2016_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Loose_2017_DR0p5", &in->bjetDeepCSVVeto20Loose_2017_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Loose_2018_DR0p5", &in->bjetDeepCSVVeto20Loose_2018_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Medium_2016_DR0", &in->bjetDeepCSVVeto20Medium_2016_DR0);
  tree->Branch("bjetDeepCSVVeto20Medium_2016_DR0p5", &in->bjetDeepCSVVeto20Medium_2016_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Medium_2017_DR0", &in->bjetDeepCSVVeto20Medium_2017_DR0);
  tree->Branch("bjetDeepCSVVeto20Medium_2017_DR0p5", &in->bjetDeepCSVVeto20Medium_2017_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Medium_2018_DR0", &in->bjetDeepCSVVeto20Medium_2018_DR0);
  tree->Branch("bjetDeepCSVVeto20Medium_2018_DR0p5", &in->bjetDeepCSVVeto20Medium_2018_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Tight_2016_DR0p5", &in->bjetDeepCSVVeto20Tight_2016_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Tight_2017_DR0p5", &in->bjetDeepCSVVeto20Tight_2017_DR0p5);
  tree->Branch("bjetDeepCSVVeto20Tight_2018_DR0p5", &in->bjetDeepCSVVeto20Tight_2018_DR0p5);
  tree->Branch("charge", &in->charge);
  tree->Branch("dielectronVeto", &in->dielectronVeto);
  tree->Branch("dimu9ele9Pass", &in->dimu9ele9Pass);
  tree->Branch("dimuonVeto", &in->dimuonVeto);
  tree->Branch("doubleMuDZminMass3p8Pass", &in->doubleMuDZminMass3p8Pass);
  tree->Branch("doubleMuDZminMass8Pass", &in->doubleMuDZminMass8Pass);
  tree->Branch("doubleMuSingleEPass", &in->doubleMuSingleEPass);
  tree->Branch("eVetoHZZPt5", &in->eVetoHZZPt5);
  tree->Branch("eVetoZTTp001dxyz", &in->eVetoZTTp001dxyz);
  tree->Branch("eVetoZTTp001dxyzR0", &in->eVetoZTTp001dxyzR0);
  tree->Branch("evt", &in->evt);
  tree->Branch("genEta", &in->genEta);
  tree->Branch("genHTT", &in->genHTT);
  tree->Branch("genM", &in->genM);
  tree->Branch("genMass", &in->genMass);
  tree->Branch("genPhi", &in->genPhi);
  tree->Branch("genpT", &in->genpT);
  tree->Branch("genpX", &in->genpX);
  tree->Branch("genpY", &in->genpY);
  tree->Branch("isdata", &in->isdata);
  tree->Branch("isembed", &in->isembed);
  tree->Branch("j1csv", &in->j1csv);
  tree->Branch("j1csvWoNoisyJets", &in->j1csvWoNoisyJets);
  tree->Branch("j1eta", &in->j1eta);
  tree->Branch("j1etaWoNoisyJets", &in->j1etaWoNoisyJets);
  tree->Branch("j1hadronflavor", &in->j1hadronflavor);
  tree->Branch("j1hadronflavorWoNoisyJets", &in->j1hadronflavorWoNoisyJets);
  tree->Branch("j1phi", &in->j1phi);
  tree->Branch("j1phiWoNoisyJets", &in->j1phiWoNoisyJets);
  tree->Branch("j1pt", &in->j1pt);
  tree->Branch("j1ptWoNoisyJets", &in->j1ptWoNoisyJets);
  tree->Branch("j1ptWoNoisyJets_JetEC2Down", &in->j1ptWoNoisyJets_JetEC2Down);
  tree->Branch("j1ptWoNoisyJets_JetEC2Up", &in->j1ptWoNoisyJets_JetEC2Up);
  tree->Branch("j1ptWoNoisyJets_JetEta0to3Down", &in->j1ptWoNoisyJets_JetEta0to3Down);
  tree->Branch("j1ptWoNoisyJets_JetEta0to3Up", &in->j1ptWoNoisyJets_JetEta0to3Up);
  tree->Branch("j1ptWoNoisyJets_JetEta0to5Down", &in->j1ptWoNoisyJets_JetEta0to5Down);
  tree->Branch("j1ptWoNoisyJets_JetEta0to5Up", &in->j1ptWoNoisyJets_JetEta0to5Up);
  tree->Branch("j1ptWoNoisyJets_JetEta3to5Down", &in->j1ptWoNoisyJets_JetEta3to5Down);
  tree->Branch("j1ptWoNoisyJets_JetEta3to5Up", &in->j1ptWoNoisyJets_JetEta3to5Up);
  tree->Branch("j1ptWoNoisyJets_JetRelativeBalDown", &in->j1ptWoNoisyJets_JetRelativeBalDown);
  tree->Branch("j1ptWoNoisyJets_JetRelativeBalUp", &in->j1ptWoNoisyJets_JetRelativeBalUp);
  tree->Branch("j1ptWoNoisyJets_JetRelativeSampleDown", &in->j1ptWoNoisyJets_JetRelativeSampleDown);
  tree->Branch("j1ptWoNoisyJets_JetRelativeSampleUp", &in->j1ptWoNoisyJets_JetRelativeSampleUp);
  tree->Branch("j1pt_JetEC2Down", &in->j1pt_JetEC2Down);
  tree->Branch("j1pt_JetEC2Up", &in->j1pt_JetEC2Up);
  tree->Branch("j1pt_JetEta0to3Down", &in->j1pt_JetEta0to3Down);
  tree->Branch("j1pt_JetEta0to3Up", &in->j1pt_JetEta0to3Up);
  tree->Branch("j1pt_JetEta0to5Down", &in->j1pt_JetEta0to5Down);
  tree->Branch("j1pt_JetEta0to5Up", &in->j1pt_JetEta0to5Up);
  tree->Branch("j1pt_JetEta3to5Down", &in->j1pt_JetEta3to5Down);
  tree->Branch("j1pt_JetEta3to5Up", &in->j1pt_JetEta3to5Up);
  tree->Branch("j1pt_JetRelativeBalDown", &in->j1pt_JetRelativeBalDown);
  tree->Branch("j1pt_JetRelativeBalUp", &in->j1pt_JetRelativeBalUp);
  tree->Branch("j1pt_JetRelativeSampleDown", &in->j1pt_JetRelativeSampleDown);
  tree->Branch("j1pt_JetRelativeSampleUp", &in->j1pt_JetRelativeSampleUp);
  tree->Branch("j2csv", &in->j2csv);
  tree->Branch("j2csvWoNoisyJets", &in->j2csvWoNoisyJets);
  tree->Branch("j2eta", &in->j2eta);
  tree->Branch("j2etaWoNoisyJets", &in->j2etaWoNoisyJets);
  tree->Branch("j2hadronflavor", &in->j2hadronflavor);
  tree->Branch("j2hadronflavorWoNoisyJets", &in->j2hadronflavorWoNoisyJets);
  tree->Branch("j2phi", &in->j2phi);
  tree->Branch("j2phiWoNoisyJets", &in->j2phiWoNoisyJets);
  tree->Branch("j2pt", &in->j2pt);
  tree->Branch("j2ptWoNoisyJets", &in->j2ptWoNoisyJets);
  tree->Branch("j2ptWoNoisyJets_JetEC2Down", &in->j2ptWoNoisyJets_JetEC2Down);
  tree->Branch("j2ptWoNoisyJets_JetEC2Up", &in->j2ptWoNoisyJets_JetEC2Up);
  tree->Branch("j2ptWoNoisyJets_JetEta0to3Down", &in->j2ptWoNoisyJets_JetEta0to3Down);
  tree->Branch("j2ptWoNoisyJets_JetEta0to3Up", &in->j2ptWoNoisyJets_JetEta0to3Up);
  tree->Branch("j2ptWoNoisyJets_JetEta0to5Down", &in->j2ptWoNoisyJets_JetEta0to5Down);
  tree->Branch("j2ptWoNoisyJets_JetEta0to5Up", &in->j2ptWoNoisyJets_JetEta0to5Up);
  tree->Branch("j2ptWoNoisyJets_JetEta3to5Down", &in->j2ptWoNoisyJets_JetEta3to5Down);
  tree->Branch("j2ptWoNoisyJets_JetEta3to5Up", &in->j2ptWoNoisyJets_JetEta3to5Up);
  tree->Branch("j2ptWoNoisyJets_JetRelativeBalDown", &in->j2ptWoNoisyJets_JetRelativeBalDown);
  tree->Branch("j2ptWoNoisyJets_JetRelativeBalUp", &in->j2ptWoNoisyJets_JetRelativeBalUp);
  tree->Branch("j2ptWoNoisyJets_JetRelativeSampleDown", &in->j2ptWoNoisyJets_JetRelativeSampleDown);
  tree->Branch("j2ptWoNoisyJets_JetRelativeSampleUp", &in->j2ptWoNoisyJets_JetRelativeSampleUp);
  tree->Branch("j2pt_JetEC2Down", &in->j2pt_JetEC2Down);
  tree->Branch("j2pt_JetEC2Up", &in->j2pt_JetEC2Up);
  tree->Branch("j2pt_JetEta0to3Down", &in->j2pt_JetEta0to3Down);
  tree->Branch("j2pt_JetEta0to3Up", &in->j2pt_JetEta0to3Up);
  tree->Branch("j2pt_JetEta0to5Down", &in->j2pt_JetEta0to5Down);
  tree->Branch("j2pt_JetEta0to5Up", &in->j2pt_JetEta0to5Up);
  tree->Branch("j2pt_JetEta3to5Down", &in->j2pt_JetEta3to5Down);
  tree->Branch("j2pt_JetEta3to5Up", &in->j2pt_JetEta3to5Up);
  tree->Branch("j2pt_JetRelativeBalDown", &in->j2pt_JetRelativeBalDown);
  tree->Branch("j2pt_JetRelativeBalUp", &in->j2pt_JetRelativeBalUp);
  tree->Branch("j2pt_JetRelativeSampleDown", &in->j2pt_JetRelativeSampleDown);
  tree->Branch("j2pt_JetRelativeSampleUp", &in->j2pt_JetRelativeSampleUp);
  tree->Branch("jb1eta_2016", &in->jb1eta_2016);
  tree->Branch("jb1eta_2017", &in->jb1eta_2017);
  tree->Branch("jb1eta_2018", &in->jb1eta_2018);
  tree->Branch("jb1hadronflavor_2016", &in->jb1hadronflavor_2016);
  tree->Branch("jb1hadronflavor_2017", &in->jb1hadronflavor_2017);
  tree->Branch("jb1hadronflavor_2018", &in->jb1hadronflavor_2018);
  tree->Branch("jb1phi_2016", &in->jb1phi_2016);
  tree->Branch("jb1phi_2017", &in->jb1phi_2017);
  tree->Branch("jb1phi_2018", &in->jb1phi_2018);
  tree->Branch("jb1pt_2016", &in->jb1pt_2016);
  tree->Branch("jb1pt_2017", &in->jb1pt_2017);
  tree->Branch("jb1pt_2018", &in->jb1pt_2018);
  tree->Branch("jb2eta_2016", &in->jb2eta_2016);
  tree->Branch("jb2eta_2017", &in->jb2eta_2017);
  tree->Branch("jb2eta_2018", &in->jb2eta_2018);
  tree->Branch("jb2hadronflavor_2016", &in->jb2hadronflavor_2016);
  tree->Branch("jb2hadronflavor_2017", &in->jb2hadronflavor_2017);
  tree->Branch("jb2hadronflavor_2018", &in->jb2hadronflavor_2018);
  tree->Branch("jb2phi_2016", &in->jb2phi_2016);
  tree->Branch("jb2phi_2017", &in->jb2phi_2017);
  tree->Branch("jb2phi_2018", &in->jb2phi_2018);
  tree->Branch("jb2pt_2016", &in->jb2pt_2016);
  tree->Branch("jb2pt_2017", &in->jb2pt_2017);
  tree->Branch("jb2pt_2018", &in->jb2pt_2018);
  tree->Branch("jetVeto20", &in->jetVeto20);
  tree->Branch("jetVeto20_JetEnDown", &in->jetVeto20_JetEnDown);
  tree->Branch("jetVeto20_JetEnUp", &in->jetVeto20_JetEnUp);
  tree->Branch("jetVeto30", &in->jetVeto30);
  tree->Branch("jetVeto30WoNoisyJets_JetEC2Down", &in->jetVeto30WoNoisyJets_JetEC2Down);
  tree->Branch("jetVeto30WoNoisyJets_JetEC2Up", &in->jetVeto30WoNoisyJets_JetEC2Up);
  tree->Branch("jetVeto30WoNoisyJets_JetEta0to3Down", &in->jetVeto30WoNoisyJets_JetEta0to3Down);
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
  tree->Branch("jetVeto30_JetEC2Down", &in->jetVeto30_JetEC2Down);
  tree->Branch("jetVeto30_JetEC2Up", &in->jetVeto30_JetEC2Up);
  tree->Branch("jetVeto30_JetEnDown", &in->jetVeto30_JetEnDown);
  tree->Branch("jetVeto30_JetEnUp", &in->jetVeto30_JetEnUp);
  tree->Branch("jetVeto30_JetEta0to3Down", &in->jetVeto30_JetEta0to3Down);
  tree->Branch("jetVeto30_JetEta0to3Up", &in->jetVeto30_JetEta0to3Up);
  tree->Branch("jetVeto30_JetEta0to5Down", &in->jetVeto30_JetEta0to5Down);
  tree->Branch("jetVeto30_JetEta0to5Up", &in->jetVeto30_JetEta0to5Up);
  tree->Branch("jetVeto30_JetEta3to5Down", &in->jetVeto30_JetEta3to5Down);
  tree->Branch("jetVeto30_JetEta3to5Up", &in->jetVeto30_JetEta3to5Up);
  tree->Branch("jetVeto30_JetRelativeBalDown", &in->jetVeto30_JetRelativeBalDown);
  tree->Branch("jetVeto30_JetRelativeBalUp", &in->jetVeto30_JetRelativeBalUp);
  tree->Branch("jetVeto30_JetRelativeSampleDown", &in->jetVeto30_JetRelativeSampleDown);
  tree->Branch("jetVeto30_JetRelativeSampleUp", &in->jetVeto30_JetRelativeSampleUp);
  tree->Branch("jetVeto30_JetTotalDown", &in->jetVeto30_JetTotalDown);
  tree->Branch("jetVeto30_JetTotalUp", &in->jetVeto30_JetTotalUp);
  tree->Branch("lumi", &in->lumi);
  tree->Branch("mBestTrackType", &in->mBestTrackType);
  tree->Branch("mCharge", &in->mCharge);
  tree->Branch("mChi2LocalPosition", &in->mChi2LocalPosition);
  tree->Branch("mComesFromHiggs", &in->mComesFromHiggs);
  tree->Branch("mCutBasedIdGlobalHighPt", &in->mCutBasedIdGlobalHighPt);
  tree->Branch("mCutBasedIdLoose", &in->mCutBasedIdLoose);
  tree->Branch("mCutBasedIdMedium", &in->mCutBasedIdMedium);
  tree->Branch("mCutBasedIdMediumPrompt", &in->mCutBasedIdMediumPrompt);
  tree->Branch("mCutBasedIdTight", &in->mCutBasedIdTight);
  tree->Branch("mCutBasedIdTrkHighPt", &in->mCutBasedIdTrkHighPt);
  tree->Branch("mEcalIsoDR03", &in->mEcalIsoDR03);
  tree->Branch("mEffectiveArea2011", &in->mEffectiveArea2011);
  tree->Branch("mEffectiveArea2012", &in->mEffectiveArea2012);
  tree->Branch("mEta", &in->mEta);
  tree->Branch("mEta_MuonEnDown", &in->mEta_MuonEnDown);
  tree->Branch("mEta_MuonEnUp", &in->mEta_MuonEnUp);
  tree->Branch("mGenCharge", &in->mGenCharge);
  tree->Branch("mGenDirectPromptTauDecayFinalState", &in->mGenDirectPromptTauDecayFinalState);
  tree->Branch("mGenEnergy", &in->mGenEnergy);
  tree->Branch("mGenEta", &in->mGenEta);
  tree->Branch("mGenIsPrompt", &in->mGenIsPrompt);
  tree->Branch("mGenMotherPdgId", &in->mGenMotherPdgId);
  tree->Branch("mGenParticle", &in->mGenParticle);
  tree->Branch("mGenPdgId", &in->mGenPdgId);
  tree->Branch("mGenPhi", &in->mGenPhi);
  tree->Branch("mGenPrompt", &in->mGenPrompt);
  tree->Branch("mGenPromptFinalState", &in->mGenPromptFinalState);
  tree->Branch("mGenPromptTauDecay", &in->mGenPromptTauDecay);
  tree->Branch("mGenPt", &in->mGenPt);
  tree->Branch("mGenTauDecay", &in->mGenTauDecay);
  tree->Branch("mGenVZ", &in->mGenVZ);
  tree->Branch("mGenVtxPVMatch", &in->mGenVtxPVMatch);
  tree->Branch("mHcalIsoDR03", &in->mHcalIsoDR03);
  tree->Branch("mIP3D", &in->mIP3D);
  tree->Branch("mIP3DErr", &in->mIP3DErr);
  tree->Branch("mIsGlobal", &in->mIsGlobal);
  tree->Branch("mIsPFMuon", &in->mIsPFMuon);
  tree->Branch("mIsTracker", &in->mIsTracker);
  tree->Branch("mIsoDB03", &in->mIsoDB03);
  tree->Branch("mIsoDB04", &in->mIsoDB04);
  tree->Branch("mJetArea", &in->mJetArea);
  tree->Branch("mJetBtag", &in->mJetBtag);
  tree->Branch("mJetDR", &in->mJetDR);
  tree->Branch("mJetEtaEtaMoment", &in->mJetEtaEtaMoment);
  tree->Branch("mJetEtaPhiMoment", &in->mJetEtaPhiMoment);
  tree->Branch("mJetEtaPhiSpread", &in->mJetEtaPhiSpread);
  tree->Branch("mJetHadronFlavour", &in->mJetHadronFlavour);
  tree->Branch("mJetPFCISVBtag", &in->mJetPFCISVBtag);
  tree->Branch("mJetPartonFlavour", &in->mJetPartonFlavour);
  tree->Branch("mJetPhiPhiMoment", &in->mJetPhiPhiMoment);
  tree->Branch("mJetPt", &in->mJetPt);
  tree->Branch("mLowestMll", &in->mLowestMll);
  tree->Branch("mMass", &in->mMass);
  tree->Branch("mMatchedStations", &in->mMatchedStations);
  tree->Branch("mMatchesIsoMu19Tau20Filter", &in->mMatchesIsoMu19Tau20Filter);
  tree->Branch("mMatchesIsoMu19Tau20Path", &in->mMatchesIsoMu19Tau20Path);
  tree->Branch("mMatchesIsoMu19Tau20SingleL1Filter", &in->mMatchesIsoMu19Tau20SingleL1Filter);
  tree->Branch("mMatchesIsoMu19Tau20SingleL1Path", &in->mMatchesIsoMu19Tau20SingleL1Path);
  tree->Branch("mMatchesIsoMu20HPSTau27Filter", &in->mMatchesIsoMu20HPSTau27Filter);
  tree->Branch("mMatchesIsoMu20HPSTau27Path", &in->mMatchesIsoMu20HPSTau27Path);
  tree->Branch("mMatchesIsoMu20Tau27Filter", &in->mMatchesIsoMu20Tau27Filter);
  tree->Branch("mMatchesIsoMu20Tau27Path", &in->mMatchesIsoMu20Tau27Path);
  tree->Branch("mMatchesIsoMu22Filter", &in->mMatchesIsoMu22Filter);
  tree->Branch("mMatchesIsoMu22Path", &in->mMatchesIsoMu22Path);
  tree->Branch("mMatchesIsoMu22eta2p1Filter", &in->mMatchesIsoMu22eta2p1Filter);
  tree->Branch("mMatchesIsoMu22eta2p1Path", &in->mMatchesIsoMu22eta2p1Path);
  tree->Branch("mMatchesIsoMu24Filter", &in->mMatchesIsoMu24Filter);
  tree->Branch("mMatchesIsoMu24Path", &in->mMatchesIsoMu24Path);
  tree->Branch("mMatchesIsoMu27Filter", &in->mMatchesIsoMu27Filter);
  tree->Branch("mMatchesIsoMu27Path", &in->mMatchesIsoMu27Path);
  tree->Branch("mMatchesIsoTkMu22Filter", &in->mMatchesIsoTkMu22Filter);
  tree->Branch("mMatchesIsoTkMu22Path", &in->mMatchesIsoTkMu22Path);
  tree->Branch("mMatchesIsoTkMu22eta2p1Filter", &in->mMatchesIsoTkMu22eta2p1Filter);
  tree->Branch("mMatchesIsoTkMu22eta2p1Path", &in->mMatchesIsoTkMu22eta2p1Path);
  tree->Branch("mMiniIsoLoose", &in->mMiniIsoLoose);
  tree->Branch("mMiniIsoMedium", &in->mMiniIsoMedium);
  tree->Branch("mMiniIsoTight", &in->mMiniIsoTight);
  tree->Branch("mMiniIsoVeryTight", &in->mMiniIsoVeryTight);
  tree->Branch("mMuonHits", &in->mMuonHits);
  tree->Branch("mMvaLoose", &in->mMvaLoose);
  tree->Branch("mMvaMedium", &in->mMvaMedium);
  tree->Branch("mMvaTight", &in->mMvaTight);
  tree->Branch("mNearestZMass", &in->mNearestZMass);
  tree->Branch("mNormTrkChi2", &in->mNormTrkChi2);
  tree->Branch("mNormalizedChi2", &in->mNormalizedChi2);
  tree->Branch("mPFChargedHadronIsoR04", &in->mPFChargedHadronIsoR04);
  tree->Branch("mPFChargedIso", &in->mPFChargedIso);
  tree->Branch("mPFIDLoose", &in->mPFIDLoose);
  tree->Branch("mPFIDMedium", &in->mPFIDMedium);
  tree->Branch("mPFIDTight", &in->mPFIDTight);
  tree->Branch("mPFIsoLoose", &in->mPFIsoLoose);
  tree->Branch("mPFIsoMedium", &in->mPFIsoMedium);
  tree->Branch("mPFIsoTight", &in->mPFIsoTight);
  tree->Branch("mPFIsoVeryLoose", &in->mPFIsoVeryLoose);
  tree->Branch("mPFIsoVeryTight", &in->mPFIsoVeryTight);
  tree->Branch("mPFNeutralHadronIsoR04", &in->mPFNeutralHadronIsoR04);
  tree->Branch("mPFNeutralIso", &in->mPFNeutralIso);
  tree->Branch("mPFPUChargedIso", &in->mPFPUChargedIso);
  tree->Branch("mPFPhotonIso", &in->mPFPhotonIso);
  tree->Branch("mPFPhotonIsoR04", &in->mPFPhotonIsoR04);
  tree->Branch("mPFPileupIsoR04", &in->mPFPileupIsoR04);
  tree->Branch("mPVDXY", &in->mPVDXY);
  tree->Branch("mPVDZ", &in->mPVDZ);
  tree->Branch("mPhi", &in->mPhi);
  tree->Branch("mPhi_MuonEnDown", &in->mPhi_MuonEnDown);
  tree->Branch("mPhi_MuonEnUp", &in->mPhi_MuonEnUp);
  tree->Branch("mPixHits", &in->mPixHits);
  tree->Branch("mPt", &in->mPt);
  tree->Branch("mPt_MuonEnDown", &in->mPt_MuonEnDown);
  tree->Branch("mPt_MuonEnUp", &in->mPt_MuonEnUp);
  tree->Branch("mRelPFIsoDBDefault", &in->mRelPFIsoDBDefault);
  tree->Branch("mRelPFIsoDBDefaultR04", &in->mRelPFIsoDBDefaultR04);
  tree->Branch("mRelPFIsoRho", &in->mRelPFIsoRho);
  tree->Branch("mRho", &in->mRho);
  tree->Branch("mSIP2D", &in->mSIP2D);
  tree->Branch("mSIP3D", &in->mSIP3D);
  tree->Branch("mSegmentCompatibility", &in->mSegmentCompatibility);
  tree->Branch("mSoftCutBasedId", &in->mSoftCutBasedId);
  tree->Branch("mTkIsoLoose", &in->mTkIsoLoose);
  tree->Branch("mTkIsoTight", &in->mTkIsoTight);
  tree->Branch("mTkLayersWithMeasurement", &in->mTkLayersWithMeasurement);
  tree->Branch("mTrkIsoDR03", &in->mTrkIsoDR03);
  tree->Branch("mTrkKink", &in->mTrkKink);
  tree->Branch("mTypeCode", &in->mTypeCode);
  tree->Branch("mVZ", &in->mVZ);
  tree->Branch("mValidFraction", &in->mValidFraction);
  tree->Branch("mZTTGenMatching", &in->mZTTGenMatching);
  tree->Branch("m_t_DR", &in->m_t_DR);
  tree->Branch("m_t_Mass", &in->m_t_Mass);
  tree->Branch("m_t_doubleL1IsoTauMatch", &in->m_t_doubleL1IsoTauMatch);
  tree->Branch("metSig", &in->metSig);
  tree->Branch("metcov00", &in->metcov00);
  tree->Branch("metcov01", &in->metcov01);
  tree->Branch("metcov10", &in->metcov10);
  tree->Branch("metcov11", &in->metcov11);
  tree->Branch("mu12e23DZPass", &in->mu12e23DZPass);
  tree->Branch("mu12e23Pass", &in->mu12e23Pass);
  tree->Branch("mu23e12DZPass", &in->mu23e12DZPass);
  tree->Branch("mu23e12Pass", &in->mu23e12Pass);
  tree->Branch("mu8diele12DZPass", &in->mu8diele12DZPass);
  tree->Branch("mu8diele12Pass", &in->mu8diele12Pass);
  tree->Branch("mu8e23DZPass", &in->mu8e23DZPass);
  tree->Branch("mu8e23Pass", &in->mu8e23Pass);
  tree->Branch("muGlbIsoVetoPt10", &in->muGlbIsoVetoPt10);
  tree->Branch("muVeto5", &in->muVeto5);
  tree->Branch("muVetoZTTp001dxyz", &in->muVetoZTTp001dxyz);
  tree->Branch("muVetoZTTp001dxyzR0", &in->muVetoZTTp001dxyzR0);
  tree->Branch("nTruePU", &in->nTruePU);
  tree->Branch("npNLO", &in->npNLO);
  tree->Branch("numGenJets", &in->numGenJets);
  tree->Branch("nvtx", &in->nvtx);
  tree->Branch("processID", &in->processID);
  tree->Branch("puppiMetEt", &in->puppiMetEt);
  tree->Branch("puppiMetPhi", &in->puppiMetPhi);
  tree->Branch("pvChi2", &in->pvChi2);
  tree->Branch("pvDX", &in->pvDX);
  tree->Branch("pvDY", &in->pvDY);
  tree->Branch("pvDZ", &in->pvDZ);
  tree->Branch("pvIsFake", &in->pvIsFake);
  tree->Branch("pvIsValid", &in->pvIsValid);
  tree->Branch("pvNormChi2", &in->pvNormChi2);
  tree->Branch("pvRho", &in->pvRho);
  tree->Branch("pvX", &in->pvX);
  tree->Branch("pvY", &in->pvY);
  tree->Branch("pvZ", &in->pvZ);
  tree->Branch("pvndof", &in->pvndof);
  tree->Branch("raw_pfMetEt", &in->raw_pfMetEt);
  tree->Branch("raw_pfMetPhi", &in->raw_pfMetPhi);
  tree->Branch("recoilDaught", &in->recoilDaught);
  tree->Branch("recoilWithMet", &in->recoilWithMet);
  tree->Branch("rho", &in->rho);
  tree->Branch("run", &in->run);
  tree->Branch("singleIsoMu22Pass", &in->singleIsoMu22Pass);
  tree->Branch("singleIsoMu22eta2p1Pass", &in->singleIsoMu22eta2p1Pass);
  tree->Branch("singleIsoTkMu22Pass", &in->singleIsoTkMu22Pass);
  tree->Branch("singleIsoTkMu22eta2p1Pass", &in->singleIsoTkMu22eta2p1Pass);
  tree->Branch("singleMu19eta2p1LooseTau20Pass", &in->singleMu19eta2p1LooseTau20Pass);
  tree->Branch("singleMu19eta2p1LooseTau20singleL1Pass", &in->singleMu19eta2p1LooseTau20singleL1Pass);
  tree->Branch("tAgainstElectronLooseMVA6", &in->tAgainstElectronLooseMVA6);
  tree->Branch("tAgainstElectronLooseMVA62018", &in->tAgainstElectronLooseMVA62018);
  tree->Branch("tAgainstElectronMVA6Raw", &in->tAgainstElectronMVA6Raw);
  tree->Branch("tAgainstElectronMVA6Raw2018", &in->tAgainstElectronMVA6Raw2018);
  tree->Branch("tAgainstElectronMVA6category", &in->tAgainstElectronMVA6category);
  tree->Branch("tAgainstElectronMVA6category2018", &in->tAgainstElectronMVA6category2018);
  tree->Branch("tAgainstElectronMediumMVA6", &in->tAgainstElectronMediumMVA6);
  tree->Branch("tAgainstElectronMediumMVA62018", &in->tAgainstElectronMediumMVA62018);
  tree->Branch("tAgainstElectronTightMVA6", &in->tAgainstElectronTightMVA6);
  tree->Branch("tAgainstElectronTightMVA62018", &in->tAgainstElectronTightMVA62018);
  tree->Branch("tAgainstElectronVLooseMVA6", &in->tAgainstElectronVLooseMVA6);
  tree->Branch("tAgainstElectronVLooseMVA62018", &in->tAgainstElectronVLooseMVA62018);
  tree->Branch("tAgainstElectronVTightMVA6", &in->tAgainstElectronVTightMVA6);
  tree->Branch("tAgainstElectronVTightMVA62018", &in->tAgainstElectronVTightMVA62018);
  tree->Branch("tAgainstMuonLoose3", &in->tAgainstMuonLoose3);
  tree->Branch("tAgainstMuonTight3", &in->tAgainstMuonTight3);
  tree->Branch("tByCombinedIsolationDeltaBetaCorrRaw3Hits", &in->tByCombinedIsolationDeltaBetaCorrRaw3Hits);
  tree->Branch("tByIsolationMVArun2v1DBdR03oldDMwLTraw", &in->tByIsolationMVArun2v1DBdR03oldDMwLTraw);
  tree->Branch("tByIsolationMVArun2v1DBnewDMwLTraw", &in->tByIsolationMVArun2v1DBnewDMwLTraw);
  tree->Branch("tByIsolationMVArun2v1DBoldDMwLTraw", &in->tByIsolationMVArun2v1DBoldDMwLTraw);
  tree->Branch("tByLooseCombinedIsolationDeltaBetaCorr3Hits", &in->tByLooseCombinedIsolationDeltaBetaCorr3Hits);
  tree->Branch("tByLooseIsolationMVArun2v1DBdR03oldDMwLT", &in->tByLooseIsolationMVArun2v1DBdR03oldDMwLT);
  tree->Branch("tByLooseIsolationMVArun2v1DBnewDMwLT", &in->tByLooseIsolationMVArun2v1DBnewDMwLT);
  tree->Branch("tByLooseIsolationMVArun2v1DBoldDMwLT", &in->tByLooseIsolationMVArun2v1DBoldDMwLT);
  tree->Branch("tByMediumCombinedIsolationDeltaBetaCorr3Hits", &in->tByMediumCombinedIsolationDeltaBetaCorr3Hits);
  tree->Branch("tByMediumIsolationMVArun2v1DBdR03oldDMwLT", &in->tByMediumIsolationMVArun2v1DBdR03oldDMwLT);
  tree->Branch("tByMediumIsolationMVArun2v1DBnewDMwLT", &in->tByMediumIsolationMVArun2v1DBnewDMwLT);
  tree->Branch("tByMediumIsolationMVArun2v1DBoldDMwLT", &in->tByMediumIsolationMVArun2v1DBoldDMwLT);
  tree->Branch("tByPhotonPtSumOutsideSignalCone", &in->tByPhotonPtSumOutsideSignalCone);
  tree->Branch("tByTightCombinedIsolationDeltaBetaCorr3Hits", &in->tByTightCombinedIsolationDeltaBetaCorr3Hits);
  tree->Branch("tByTightIsolationMVArun2v1DBdR03oldDMwLT", &in->tByTightIsolationMVArun2v1DBdR03oldDMwLT);
  tree->Branch("tByTightIsolationMVArun2v1DBnewDMwLT", &in->tByTightIsolationMVArun2v1DBnewDMwLT);
  tree->Branch("tByTightIsolationMVArun2v1DBoldDMwLT", &in->tByTightIsolationMVArun2v1DBoldDMwLT);
  tree->Branch("tByVLooseIsolationMVArun2v1DBdR03oldDMwLT", &in->tByVLooseIsolationMVArun2v1DBdR03oldDMwLT);
  tree->Branch("tByVLooseIsolationMVArun2v1DBnewDMwLT", &in->tByVLooseIsolationMVArun2v1DBnewDMwLT);
  tree->Branch("tByVLooseIsolationMVArun2v1DBoldDMwLT", &in->tByVLooseIsolationMVArun2v1DBoldDMwLT);
  tree->Branch("tByVTightIsolationMVArun2v1DBdR03oldDMwLT", &in->tByVTightIsolationMVArun2v1DBdR03oldDMwLT);
  tree->Branch("tByVTightIsolationMVArun2v1DBnewDMwLT", &in->tByVTightIsolationMVArun2v1DBnewDMwLT);
  tree->Branch("tByVTightIsolationMVArun2v1DBoldDMwLT", &in->tByVTightIsolationMVArun2v1DBoldDMwLT);
  tree->Branch("tByVVTightIsolationMVArun2v1DBdR03oldDMwLT", &in->tByVVTightIsolationMVArun2v1DBdR03oldDMwLT);
  tree->Branch("tByVVTightIsolationMVArun2v1DBnewDMwLT", &in->tByVVTightIsolationMVArun2v1DBnewDMwLT);
  tree->Branch("tByVVTightIsolationMVArun2v1DBoldDMwLT", &in->tByVVTightIsolationMVArun2v1DBoldDMwLT);
  tree->Branch("tCharge", &in->tCharge);
  tree->Branch("tChargedIsoPtSum", &in->tChargedIsoPtSum);
  tree->Branch("tChargedIsoPtSumdR03", &in->tChargedIsoPtSumdR03);
  tree->Branch("tComesFromHiggs", &in->tComesFromHiggs);
  tree->Branch("tDecayMode", &in->tDecayMode);
  tree->Branch("tDecayModeFinding", &in->tDecayModeFinding);
  tree->Branch("tDecayModeFindingNewDMs", &in->tDecayModeFindingNewDMs);
  tree->Branch("tDeepTau2017v1VSeraw", &in->tDeepTau2017v1VSeraw);
  tree->Branch("tDeepTau2017v1VSjetraw", &in->tDeepTau2017v1VSjetraw);
  tree->Branch("tDeepTau2017v1VSmuraw", &in->tDeepTau2017v1VSmuraw);
  tree->Branch("tDpfTau2016v0VSallraw", &in->tDpfTau2016v0VSallraw);
  tree->Branch("tDpfTau2016v1VSallraw", &in->tDpfTau2016v1VSallraw);
  tree->Branch("tEta", &in->tEta);
  tree->Branch("tFootprintCorrection", &in->tFootprintCorrection);
  tree->Branch("tFootprintCorrectiondR03", &in->tFootprintCorrectiondR03);
  tree->Branch("tGenCharge", &in->tGenCharge);
  tree->Branch("tGenDecayMode", &in->tGenDecayMode);
  tree->Branch("tGenEnergy", &in->tGenEnergy);
  tree->Branch("tGenEta", &in->tGenEta);
  tree->Branch("tGenJetEta", &in->tGenJetEta);
  tree->Branch("tGenJetPt", &in->tGenJetPt);
  tree->Branch("tGenMotherEnergy", &in->tGenMotherEnergy);
  tree->Branch("tGenMotherEta", &in->tGenMotherEta);
  tree->Branch("tGenMotherPdgId", &in->tGenMotherPdgId);
  tree->Branch("tGenMotherPhi", &in->tGenMotherPhi);
  tree->Branch("tGenMotherPt", &in->tGenMotherPt);
  tree->Branch("tGenPdgId", &in->tGenPdgId);
  tree->Branch("tGenPhi", &in->tGenPhi);
  tree->Branch("tGenPt", &in->tGenPt);
  tree->Branch("tGenStatus", &in->tGenStatus);
  tree->Branch("tJetArea", &in->tJetArea);
  tree->Branch("tJetBtag", &in->tJetBtag);
  tree->Branch("tJetDR", &in->tJetDR);
  tree->Branch("tJetEtaEtaMoment", &in->tJetEtaEtaMoment);
  tree->Branch("tJetEtaPhiMoment", &in->tJetEtaPhiMoment);
  tree->Branch("tJetEtaPhiSpread", &in->tJetEtaPhiSpread);
  tree->Branch("tJetHadronFlavour", &in->tJetHadronFlavour);
  tree->Branch("tJetPFCISVBtag", &in->tJetPFCISVBtag);
  tree->Branch("tJetPartonFlavour", &in->tJetPartonFlavour);
  tree->Branch("tJetPhiPhiMoment", &in->tJetPhiPhiMoment);
  tree->Branch("tJetPt", &in->tJetPt);
  tree->Branch("tL1IsoTauMatch", &in->tL1IsoTauMatch);
  tree->Branch("tL1IsoTauPt", &in->tL1IsoTauPt);
  tree->Branch("tLeadTrackPt", &in->tLeadTrackPt);
  tree->Branch("tLooseDeepTau2017v1VSe", &in->tLooseDeepTau2017v1VSe);
  tree->Branch("tLooseDeepTau2017v1VSjet", &in->tLooseDeepTau2017v1VSjet);
  tree->Branch("tLooseDeepTau2017v1VSmu", &in->tLooseDeepTau2017v1VSmu);
  tree->Branch("tLowestMll", &in->tLowestMll);
  tree->Branch("tMass", &in->tMass);
  tree->Branch("tMatchesIsoMu19Tau20Filter", &in->tMatchesIsoMu19Tau20Filter);
  tree->Branch("tMatchesIsoMu19Tau20Path", &in->tMatchesIsoMu19Tau20Path);
  tree->Branch("tMatchesIsoMu19Tau20SingleL1Filter", &in->tMatchesIsoMu19Tau20SingleL1Filter);
  tree->Branch("tMatchesIsoMu19Tau20SingleL1Path", &in->tMatchesIsoMu19Tau20SingleL1Path);
  tree->Branch("tMatchesIsoMu20HPSTau27Filter", &in->tMatchesIsoMu20HPSTau27Filter);
  tree->Branch("tMatchesIsoMu20HPSTau27Path", &in->tMatchesIsoMu20HPSTau27Path);
  tree->Branch("tMatchesIsoMu20Tau27Filter", &in->tMatchesIsoMu20Tau27Filter);
  tree->Branch("tMatchesIsoMu20Tau27Path", &in->tMatchesIsoMu20Tau27Path);
  tree->Branch("tMediumDeepTau2017v1VSe", &in->tMediumDeepTau2017v1VSe);
  tree->Branch("tMediumDeepTau2017v1VSjet", &in->tMediumDeepTau2017v1VSjet);
  tree->Branch("tMediumDeepTau2017v1VSmu", &in->tMediumDeepTau2017v1VSmu);
  tree->Branch("tNChrgHadrIsolationCands", &in->tNChrgHadrIsolationCands);
  tree->Branch("tNChrgHadrSignalCands", &in->tNChrgHadrSignalCands);
  tree->Branch("tNGammaSignalCands", &in->tNGammaSignalCands);
  tree->Branch("tNNeutralHadrSignalCands", &in->tNNeutralHadrSignalCands);
  tree->Branch("tNSignalCands", &in->tNSignalCands);
  tree->Branch("tNearestZMass", &in->tNearestZMass);
  tree->Branch("tNeutralIsoPtSum", &in->tNeutralIsoPtSum);
  tree->Branch("tNeutralIsoPtSumWeight", &in->tNeutralIsoPtSumWeight);
  tree->Branch("tNeutralIsoPtSumWeightdR03", &in->tNeutralIsoPtSumWeightdR03);
  tree->Branch("tNeutralIsoPtSumdR03", &in->tNeutralIsoPtSumdR03);
  tree->Branch("tPVDXY", &in->tPVDXY);
  tree->Branch("tPVDZ", &in->tPVDZ);
  tree->Branch("tPhi", &in->tPhi);
  tree->Branch("tPhotonPtSumOutsideSignalCone", &in->tPhotonPtSumOutsideSignalCone);
  tree->Branch("tPhotonPtSumOutsideSignalConedR03", &in->tPhotonPtSumOutsideSignalConedR03);
  tree->Branch("tPt", &in->tPt);
  tree->Branch("tPuCorrPtSum", &in->tPuCorrPtSum);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTLoose", &in->tRerunMVArun2v2DBoldDMwLTLoose);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTMedium", &in->tRerunMVArun2v2DBoldDMwLTMedium);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTTight", &in->tRerunMVArun2v2DBoldDMwLTTight);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVLoose", &in->tRerunMVArun2v2DBoldDMwLTVLoose);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVTight", &in->tRerunMVArun2v2DBoldDMwLTVTight);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVVLoose", &in->tRerunMVArun2v2DBoldDMwLTVVLoose);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVVTight", &in->tRerunMVArun2v2DBoldDMwLTVVTight);
  tree->Branch("tRerunMVArun2v2DBoldDMwLTraw", &in->tRerunMVArun2v2DBoldDMwLTraw);
  tree->Branch("tTightDeepTau2017v1VSe", &in->tTightDeepTau2017v1VSe);
  tree->Branch("tTightDeepTau2017v1VSjet", &in->tTightDeepTau2017v1VSjet);
  tree->Branch("tTightDeepTau2017v1VSmu", &in->tTightDeepTau2017v1VSmu);
  tree->Branch("tTightDpfTau2016v0VSall", &in->tTightDpfTau2016v0VSall);
  tree->Branch("tTightDpfTau2016v1VSall", &in->tTightDpfTau2016v1VSall);
  tree->Branch("tVLooseDeepTau2017v1VSe", &in->tVLooseDeepTau2017v1VSe);
  tree->Branch("tVLooseDeepTau2017v1VSjet", &in->tVLooseDeepTau2017v1VSjet);
  tree->Branch("tVLooseDeepTau2017v1VSmu", &in->tVLooseDeepTau2017v1VSmu);
  tree->Branch("tVTightDeepTau2017v1VSe", &in->tVTightDeepTau2017v1VSe);
  tree->Branch("tVTightDeepTau2017v1VSjet", &in->tVTightDeepTau2017v1VSjet);
  tree->Branch("tVTightDeepTau2017v1VSmu", &in->tVTightDeepTau2017v1VSmu);
  tree->Branch("tVVLooseDeepTau2017v1VSe", &in->tVVLooseDeepTau2017v1VSe);
  tree->Branch("tVVLooseDeepTau2017v1VSjet", &in->tVVLooseDeepTau2017v1VSjet);
  tree->Branch("tVVLooseDeepTau2017v1VSmu", &in->tVVLooseDeepTau2017v1VSmu);
  tree->Branch("tVVTightDeepTau2017v1VSe", &in->tVVTightDeepTau2017v1VSe);
  tree->Branch("tVVTightDeepTau2017v1VSjet", &in->tVVTightDeepTau2017v1VSjet);
  tree->Branch("tVVTightDeepTau2017v1VSmu", &in->tVVTightDeepTau2017v1VSmu);
  tree->Branch("tVVVLooseDeepTau2017v1VSe", &in->tVVVLooseDeepTau2017v1VSe);
  tree->Branch("tVVVLooseDeepTau2017v1VSjet", &in->tVVVLooseDeepTau2017v1VSjet);
  tree->Branch("tVVVLooseDeepTau2017v1VSmu", &in->tVVVLooseDeepTau2017v1VSmu);
  tree->Branch("tVZ", &in->tVZ);
  tree->Branch("tZTTGenDR", &in->tZTTGenDR);
  tree->Branch("tZTTGenEta", &in->tZTTGenEta);
  tree->Branch("tZTTGenMatching", &in->tZTTGenMatching);
  tree->Branch("tZTTGenPhi", &in->tZTTGenPhi);
  tree->Branch("tZTTGenPt", &in->tZTTGenPt);
  tree->Branch("tauVetoPt20Loose3HitsVtx", &in->tauVetoPt20Loose3HitsVtx);
  tree->Branch("tauVetoPt20TightMVALTVtx", &in->tauVetoPt20TightMVALTVtx);
  tree->Branch("topQuarkPt1", &in->topQuarkPt1);
  tree->Branch("topQuarkPt2", &in->topQuarkPt2);
  tree->Branch("tripleEPass", &in->tripleEPass);
  tree->Branch("tripleMu10_5_5Pass", &in->tripleMu10_5_5Pass);
  tree->Branch("tripleMu12_10_5Pass", &in->tripleMu12_10_5Pass);
  tree->Branch("type1_pfMetEt", &in->type1_pfMetEt);
  tree->Branch("type1_pfMetPhi", &in->type1_pfMetPhi);
  tree->Branch("type1_pfMet_shiftedPhi_JetEnDown", &in->type1_pfMet_shiftedPhi_JetEnDown);
  tree->Branch("type1_pfMet_shiftedPhi_JetEnUp", &in->type1_pfMet_shiftedPhi_JetEnUp);
  tree->Branch("type1_pfMet_shiftedPhi_JetEta0to3Down", &in->type1_pfMet_shiftedPhi_JetEta0to3Down);
  tree->Branch("type1_pfMet_shiftedPhi_JetEta0to3Up", &in->type1_pfMet_shiftedPhi_JetEta0to3Up);
  tree->Branch("type1_pfMet_shiftedPhi_JetEta0to5Down", &in->type1_pfMet_shiftedPhi_JetEta0to5Down);
  tree->Branch("type1_pfMet_shiftedPhi_JetEta0to5Up", &in->type1_pfMet_shiftedPhi_JetEta0to5Up);
  tree->Branch("type1_pfMet_shiftedPhi_JetEta3to5Down", &in->type1_pfMet_shiftedPhi_JetEta3to5Down);
  tree->Branch("type1_pfMet_shiftedPhi_JetEta3to5Up", &in->type1_pfMet_shiftedPhi_JetEta3to5Up);
  tree->Branch("type1_pfMet_shiftedPhi_JetRelativeBalDown", &in->type1_pfMet_shiftedPhi_JetRelativeBalDown);
  tree->Branch("type1_pfMet_shiftedPhi_JetRelativeBalUp", &in->type1_pfMet_shiftedPhi_JetRelativeBalUp);
  tree->Branch("type1_pfMet_shiftedPhi_JetRelativeSampleDown", &in->type1_pfMet_shiftedPhi_JetRelativeSampleDown);
  tree->Branch("type1_pfMet_shiftedPhi_JetRelativeSampleUp", &in->type1_pfMet_shiftedPhi_JetRelativeSampleUp);
  tree->Branch("type1_pfMet_shiftedPhi_JetResDown", &in->type1_pfMet_shiftedPhi_JetResDown);
  tree->Branch("type1_pfMet_shiftedPhi_JetResUp", &in->type1_pfMet_shiftedPhi_JetResUp);
  tree->Branch("type1_pfMet_shiftedPhi_JetTotalDown", &in->type1_pfMet_shiftedPhi_JetTotalDown);
  tree->Branch("type1_pfMet_shiftedPhi_JetTotalUp", &in->type1_pfMet_shiftedPhi_JetTotalUp);
  tree->Branch("type1_pfMet_shiftedPhi_UnclusteredEnDown", &in->type1_pfMet_shiftedPhi_UnclusteredEnDown);
  tree->Branch("type1_pfMet_shiftedPhi_UnclusteredEnUp", &in->type1_pfMet_shiftedPhi_UnclusteredEnUp);
  tree->Branch("type1_pfMet_shiftedPt_JetEnDown", &in->type1_pfMet_shiftedPt_JetEnDown);
  tree->Branch("type1_pfMet_shiftedPt_JetEnUp", &in->type1_pfMet_shiftedPt_JetEnUp);
  tree->Branch("type1_pfMet_shiftedPt_JetEta0to3Down", &in->type1_pfMet_shiftedPt_JetEta0to3Down);
  tree->Branch("type1_pfMet_shiftedPt_JetEta0to3Up", &in->type1_pfMet_shiftedPt_JetEta0to3Up);
  tree->Branch("type1_pfMet_shiftedPt_JetEta0to5Down", &in->type1_pfMet_shiftedPt_JetEta0to5Down);
  tree->Branch("type1_pfMet_shiftedPt_JetEta0to5Up", &in->type1_pfMet_shiftedPt_JetEta0to5Up);
  tree->Branch("type1_pfMet_shiftedPt_JetEta3to5Down", &in->type1_pfMet_shiftedPt_JetEta3to5Down);
  tree->Branch("type1_pfMet_shiftedPt_JetEta3to5Up", &in->type1_pfMet_shiftedPt_JetEta3to5Up);
  tree->Branch("type1_pfMet_shiftedPt_JetRelativeBalDown", &in->type1_pfMet_shiftedPt_JetRelativeBalDown);
  tree->Branch("type1_pfMet_shiftedPt_JetRelativeBalUp", &in->type1_pfMet_shiftedPt_JetRelativeBalUp);
  tree->Branch("type1_pfMet_shiftedPt_JetRelativeSampleDown", &in->type1_pfMet_shiftedPt_JetRelativeSampleDown);
  tree->Branch("type1_pfMet_shiftedPt_JetRelativeSampleUp", &in->type1_pfMet_shiftedPt_JetRelativeSampleUp);
  tree->Branch("type1_pfMet_shiftedPt_JetResDown", &in->type1_pfMet_shiftedPt_JetResDown);
  tree->Branch("type1_pfMet_shiftedPt_JetResUp", &in->type1_pfMet_shiftedPt_JetResUp);
  tree->Branch("type1_pfMet_shiftedPt_JetTotalDown", &in->type1_pfMet_shiftedPt_JetTotalDown);
  tree->Branch("type1_pfMet_shiftedPt_JetTotalUp", &in->type1_pfMet_shiftedPt_JetTotalUp);
  tree->Branch("type1_pfMet_shiftedPt_UnclusteredEnDown", &in->type1_pfMet_shiftedPt_UnclusteredEnDown);
  tree->Branch("type1_pfMet_shiftedPt_UnclusteredEnUp", &in->type1_pfMet_shiftedPt_UnclusteredEnUp);
  tree->Branch("vbfDeta", &in->vbfDeta);
  tree->Branch("vbfJetVeto20", &in->vbfJetVeto20);
  tree->Branch("vbfJetVeto30", &in->vbfJetVeto30);
  tree->Branch("vbfMass", &in->vbfMass);
  tree->Branch("vbfMassWoNoisyJets", &in->vbfMassWoNoisyJets);
  tree->Branch("vbfMassWoNoisyJets_JetEC2Down", &in->vbfMassWoNoisyJets_JetEC2Down);
  tree->Branch("vbfMassWoNoisyJets_JetEC2Up", &in->vbfMassWoNoisyJets_JetEC2Up);
  tree->Branch("vbfMassWoNoisyJets_JetEta0to3Down", &in->vbfMassWoNoisyJets_JetEta0to3Down);
  tree->Branch("vbfMassWoNoisyJets_JetEta0to3Up", &in->vbfMassWoNoisyJets_JetEta0to3Up);
  tree->Branch("vbfMassWoNoisyJets_JetEta0to5Down", &in->vbfMassWoNoisyJets_JetEta0to5Down);
  tree->Branch("vbfMassWoNoisyJets_JetEta0to5Up", &in->vbfMassWoNoisyJets_JetEta0to5Up);
  tree->Branch("vbfMassWoNoisyJets_JetEta3to5Down", &in->vbfMassWoNoisyJets_JetEta3to5Down);
  tree->Branch("vbfMassWoNoisyJets_JetEta3to5Up", &in->vbfMassWoNoisyJets_JetEta3to5Up);
  tree->Branch("vbfMassWoNoisyJets_JetRelativeBalDown", &in->vbfMassWoNoisyJets_JetRelativeBalDown);
  tree->Branch("vbfMassWoNoisyJets_JetRelativeBalUp", &in->vbfMassWoNoisyJets_JetRelativeBalUp);
  tree->Branch("vbfMassWoNoisyJets_JetRelativeSampleDown", &in->vbfMassWoNoisyJets_JetRelativeSampleDown);
  tree->Branch("vbfMassWoNoisyJets_JetRelativeSampleUp", &in->vbfMassWoNoisyJets_JetRelativeSampleUp);
  tree->Branch("vbfMassWoNoisyJets_JetTotalDown", &in->vbfMassWoNoisyJets_JetTotalDown);
  tree->Branch("vbfMassWoNoisyJets_JetTotalUp", &in->vbfMassWoNoisyJets_JetTotalUp);
  tree->Branch("vbfMass_JetEC2Down", &in->vbfMass_JetEC2Down);
  tree->Branch("vbfMass_JetEC2Up", &in->vbfMass_JetEC2Up);
  tree->Branch("vbfMass_JetEta0to3Down", &in->vbfMass_JetEta0to3Down);
  tree->Branch("vbfMass_JetEta0to3Up", &in->vbfMass_JetEta0to3Up);
  tree->Branch("vbfMass_JetEta0to5Down", &in->vbfMass_JetEta0to5Down);
  tree->Branch("vbfMass_JetEta0to5Up", &in->vbfMass_JetEta0to5Up);
  tree->Branch("vbfMass_JetEta3to5Down", &in->vbfMass_JetEta3to5Down);
  tree->Branch("vbfMass_JetEta3to5Up", &in->vbfMass_JetEta3to5Up);
  tree->Branch("vbfMass_JetRelativeBalDown", &in->vbfMass_JetRelativeBalDown);
  tree->Branch("vbfMass_JetRelativeBalUp", &in->vbfMass_JetRelativeBalUp);
  tree->Branch("vbfMass_JetRelativeSampleDown", &in->vbfMass_JetRelativeSampleDown);
  tree->Branch("vbfMass_JetRelativeSampleUp", &in->vbfMass_JetRelativeSampleUp);
  tree->Branch("vbfMass_JetTotalDown", &in->vbfMass_JetTotalDown);
  tree->Branch("vbfMass_JetTotalUp", &in->vbfMass_JetTotalUp);
  tree->Branch("vbfNJets20", &in->vbfNJets20);
  tree->Branch("vbfNJets30", &in->vbfNJets30);
  tree->Branch("vbfj1eta", &in->vbfj1eta);
  tree->Branch("vbfj1pt", &in->vbfj1pt);
  tree->Branch("vbfj2eta", &in->vbfj2eta);
  tree->Branch("vbfj2pt", &in->vbfj2pt);
  tree->Branch("vispX", &in->vispX);
  tree->Branch("vispY", &in->vispY);
  tree->Branch("idx", &in->idx);
}

#endif  // ROOT_SRC_MUTAU_TREE_2018_H_
