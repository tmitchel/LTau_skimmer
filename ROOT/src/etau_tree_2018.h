// Copyright 2019 Tyler Mitchell

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "./base_tree.h"
#include "./et_2018_input_branches.h"
#include "RecoilCorrector.h"
#include "TLorentzVector.h"
#include "TTree.h"

class etau_tree2018 : public virtual base_tree {
 private:
  TTree *tree, *original;
  et_2018_input_branches* in;
  bool isMC, isEmbed;
  std::vector<Int_t> good_events;
  TLorentzVector ele, tau, MET, MET_UESUp, MET_UESDown, MET_JESUp, MET_JESDown;

 public:
  // Member variables

  Int_t Run, Lumi, recoil;
  Float_t placeholder;  // for all branches not present in 2018

  // // Constructed while running
  Int_t gen_match_1, gen_match_2, njets, nbtag, njetspt20;
  Float_t jetVeto20, jetVeto30, met, metphi, met_px, met_py, extraelec_veto, extramuon_veto, dilepton_veto, pfmetcorr_ex, pfmetcorr_ey;
  Float_t pfmetcorr_ex_UESUp, pfmetcorr_ey_UESUp, pfmetcorr_ex_UESDown, pfmetcorr_ey_UESDown, pfmetcorr_ex_JESUp, pfmetcorr_ey_JESUp,
      pfmetcorr_ex_JESDown, pfmetcorr_ey_JESDown;
  Float_t met_UESUp, met_UESDown, met_JESUp, met_JESDown, metphi_UESUp, metphi_UESDown, metphi_JESUp, metphi_JESDown;

  Float_t pt_1, eta_1, phi_1, m_1, e_1, px_1, py_1, pz_1, pt_2, eta_2, phi_2, m_2, e_2, px_2, py_2, pz_2;

  // Member functions
  etau_tree2018(TTree* orig, TTree* itree, bool isMC, bool isEmbed, Int_t rec);
  virtual ~etau_tree2018() {}
  void do_skimming(TH1F*);
  void set_branches();
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
etau_tree2018::etau_tree2018(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec) : tree(itree),
                                                                                                  original(Original),
                                                                                                  in(new et_2018_input_branches(Original)),
                                                                                                  isMC(IsMC),
                                                                                                  isEmbed(IsEmbed),
                                                                                                  recoil(rec) {
  original->SetBranchStatus("Double*", 0);
  original->SetBranchStatus("Mu*", 0);
}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void etau_tree2018::do_skimming(TH1F* cutflow) {
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
    auto Cross = in->eMatchesEle24Tau30Filter && in->eMatchesEle24Tau30Path && in->Ele24Tau30Pass && in->tMatchesEle24Tau30Path && in->tMatchesEle24Tau30Filter;

    if (isEmbed || (Ele27 || Ele32 || Ele35 || Cross))
      cutflow->Fill(2., 1.);
    else
      continue;

    if (!isEmbed || (in->Ele27WPTightPass || in->Ele32WPTightPass || in->Ele35WPTightPass || in->Ele24Tau30Pass))
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

//////////////////////////////////////////////////////////////////
// Purpose: Fill tree with variables from original and new      //
//          variables. Only events that have been skimmed,      //
//          sorted, and stored in the good_events vector will   //
//          be stored in the tree.                              //
//////////////////////////////////////////////////////////////////
// Return: The same TTree passed to the constructor and stored  //
//         in original, but now it is filled with good events   //
//////////////////////////////////////////////////////////////////
TTree* etau_tree2018::fill_tree(RecoilCorrector recoilPFMetCorrector) {
  set_branches();  // get all the branches set up

  // loop through all events pasing skimming/sorting
  for (auto& ievt : good_events) {
    original->GetEntry(ievt);

    // convert from Float_t in FSA to Int_t for analyzer
    gen_match_1 = in->eZTTGenMatching;
    gen_match_2 = in->tZTTGenMatching;
    njets = in->jetVeto30;
    nbtag = in->bjetDeepCSVVeto20Medium;
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
    // not in ntuples yet
    // MET_UESUp.SetPtEtaPhiM(type1_pfMet_shiftedPt_UnclusteredEnUp, 0, type1_pfMet_shiftedPhi_UnclusteredEnUp, 0);
    // MET_UESDown.SetPtEtaPhiM(type1_pfMet_shiftedPt_UnclusteredEnDown, 0, type1_pfMet_shiftedPhi_UnclusteredEnDown, 0);
    // MET_JESUp.SetPtEtaPhiM(type1_pfMet_shiftedPt_JetEnUp, 0, type1_pfMet_shiftedPhi_JetEnUp, 0);
    // MET_JESDown.SetPtEtaPhiM(type1_pfMet_shiftedPt_JetEnDown, 0, type1_pfMet_shiftedPhi_JetEnDown, 0);

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
          in->genpX,          // generator Z/W/Higgs px (float)
          in->genpY,          // generator Z/W/Higgs py (float)
          in->vispX,          // generator visible Z/W/Higgs px (float)
          in->vispY,          // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                // generator Z/W/Higgs px (float)
          in->genpY,                // generator Z/W/Higgs py (float)
          in->vispX,                // generator visible Z/W/Higgs px (float)
          in->vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                // generator Z/W/Higgs px (float)
          in->genpY,                // generator Z/W/Higgs py (float)
          in->vispX,                // generator visible Z/W/Higgs px (float)
          in->vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                  // generator Z/W/Higgs px (float)
          in->genpY,                  // generator Z/W/Higgs py (float)
          in->vispX,                  // generator visible Z/W/Higgs px (float)
          in->vispY,                  // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                  // generator Z/W/Higgs px (float)
          in->genpY,                  // generator Z/W/Higgs py (float)
          in->vispX,                  // generator visible Z/W/Higgs px (float)
          in->vispY,                  // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESDown);  // corrected type I pf met py (float)

    } else if (recoil == 2) {
      recoilPFMetCorrector.CorrectByMeanResolution(
          MET.Px(),       // uncorrected type I pf met px (float)
          MET.Py(),       // uncorrected type I pf met py (float)
          in->genpX,          // generator Z/W/Higgs px (float)
          in->genpY,          // generator Z/W/Higgs py (float)
          in->vispX,          // generator visible Z/W/Higgs px (float)
          in->vispY,          // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                // generator Z/W/Higgs px (float)
          in->genpY,                // generator Z/W/Higgs py (float)
          in->vispX,                // generator visible Z/W/Higgs px (float)
          in->vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                // generator Z/W/Higgs px (float)
          in->genpY,                // generator Z/W/Higgs py (float)
          in->vispX,                // generator visible Z/W/Higgs px (float)
          in->vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                  // generator Z/W/Higgs px (float)
          in->genpY,                  // generator Z/W/Higgs py (float)
          in->vispX,                  // generator visible Z/W/Higgs px (float)
          in->vispY,                  // generator visible Z/W/Higgs py (float)
          jetVeto30,              // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          in->genpX,                  // generator Z/W/Higgs px (float)
          in->genpY,                  // generator Z/W/Higgs py (float)
          in->vispX,                  // generator visible Z/W/Higgs px (float)
          in->vispY,                  // generator visible Z/W/Higgs py (float)
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
        if (in->tDecayMode == 0) {
          MET = MET + tau - 1.007 * tau;
          MET_JESUp = MET_JESUp + tau - 1.007 * tau;
          MET_JESDown = MET_JESDown + tau - 1.007 * tau;
          MET_UESUp = MET_UESUp + tau - 1.007 * tau;
          MET_UESDown = MET_UESDown + tau - 1.007 * tau;
          tau *= 1.007;
        } else if (in->tDecayMode == 1) {
          MET = MET + tau - 0.998 * tau;
          MET_JESUp = MET_JESUp + tau - 0.998 * tau;
          MET_JESDown = MET_JESDown + tau - 0.998 * tau;
          MET_UESUp = MET_UESUp + tau - 0.998 * tau;
          MET_UESDown = MET_UESDown + tau - 0.998 * tau;
          tau *= 0.998;
        } else if (in->tDecayMode == 10) {
          MET = MET + tau - 1.001 * tau;
          MET_JESUp = MET_JESUp + tau - 1.001 * tau;
          MET_JESDown = MET_JESDown + tau - 1.001 * tau;
          MET_UESUp = MET_UESUp + tau - 1.001 * tau;
          MET_UESDown = MET_UESDown + tau - 1.001 * tau;
          tau *= 1.001;
        }
      } else if (in->tZTTGenMatching == 1 || in->tZTTGenMatching == 3) {
        if (in->tDecayMode == 0) {
          MET = MET + tau - 1.003 * tau;
          MET_JESUp = MET_JESUp + tau - 1.003 * tau;
          MET_JESDown = MET_JESDown + tau - 1.003 * tau;
          MET_UESUp = MET_UESUp + tau - 1.003 * tau;
          MET_UESDown = MET_UESDown + tau - 1.003 * tau;
          tau *= 1.003;
        } else if (in->tDecayMode == 1) {
          MET = MET + tau - 1.036 * tau;
          MET_JESUp = MET_JESUp + tau - 1.036 * tau;
          MET_JESDown = MET_JESDown + tau - 1.036 * tau;
          MET_UESUp = MET_UESUp + tau - 1.036 * tau;
          MET_UESDown = MET_UESDown + tau - 1.036 * tau;
          tau *= 1.036;
        }
      }
    } else if (isEmbed) {
      if (in->tZTTGenMatching == 5) {
        if (in->tDecayMode == 0) {
          MET = MET + tau - 0.975 * tau;
          MET_JESUp = MET_JESUp + tau - 0.975 * tau;
          MET_JESDown = MET_JESDown + tau - 0.975 * tau;
          MET_UESUp = MET_UESUp + tau - 0.975 * tau;
          MET_UESDown = MET_UESDown + tau - 0.975 * tau;
          tau *= 0.975;
        } else if (in->tDecayMode == 1) {
          MET = MET + tau - 0.975 * 1.051 * tau;
          MET_JESUp = MET_JESUp + tau - 0.975 * 1.051 * tau;
          MET_JESDown = MET_JESDown + tau - 0.975 * 1.051 * tau;
          MET_UESUp = MET_UESUp + tau - 0.975 * 1.051 * tau;
          MET_UESDown = MET_UESDown + tau - 0.975 * 1.051 * tau;
          tau *= 0.975 * 1.051;
        } else if (in->tDecayMode == 10) {
          MET = MET + tau - 0.975 * 0.975 * 0.975 * tau;
          MET_JESUp = MET_JESUp + tau - 0.975 * 0.975 * 0.975 * tau;
          MET_JESDown = MET_JESDown + tau - 0.975 * 0.975 * 0.975 * tau;
          MET_UESUp = MET_UESUp + tau - 0.975 * 0.975 * 0.975 * tau;
          MET_UESDown = MET_UESDown + tau - 0.975 * 0.975 * 0.975 * tau;
          tau *= 0.975 * 0.975 * 0.975;
        }
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
void etau_tree2018::set_branches() {
  // output file branches
  tree->Branch("evt", &in->evt);
  tree->Branch("run", &Run);
  tree->Branch("lumi", &Lumi);
  tree->Branch("gen_match_1", &gen_match_1, "gen_match_1/I");
  tree->Branch("gen_match_2", &gen_match_2, "gen_match_2/I");
  tree->Branch("njets", &njets, "njets/I");
  tree->Branch("nbtag", &nbtag, "nbtag/I");
  tree->Branch("njetspt20", &njetspt20, "njetspt20/I");
  tree->Branch("vbfMassWoNoisyJets", &placeholder, "vbfMassWoNoisyJets/F");

  tree->Branch("eMatchesEle27Filter", &placeholder, "eMatchesEle27Filter/F");
  tree->Branch("eMatchesEle32Filter", &placeholder, "eMatchesEle32Filter/F");
  tree->Branch("eMatchesEle35Filter", &placeholder, "eMatchesEle35Filter/F");
  tree->Branch("eMatchesEle24Tau30Filter", &placeholder, "eMatchesEle24Tau30Filter/F");
  tree->Branch("tMatchesEle24Tau30Filter", &placeholder, "tMatchesEle24Tau30Filter/F");
  tree->Branch("eMatchesEle27Path", &placeholder, "eMatchesEle27Path/F");
  tree->Branch("eMatchesEle32Path", &placeholder, "eMatchesEle32Path/F");
  tree->Branch("eMatchesEle35Path", &placeholder, "eMatchesEle35Path/F");
  tree->Branch("eMatchesEle24Tau30Path", &placeholder, "eMatchesEle24Tau30Path/F");
  tree->Branch("tMatchesEle24Tau30Path", &placeholder, "tMatchesEle24Tau30Path/F");
  tree->Branch("Ele24Tau30Pass", &placeholder, "Ele24Tau30Pass/F");
  tree->Branch("Ele27WPTightPass", &placeholder, "Ele27WPTightPass/F");
  tree->Branch("Ele32WPTightPass", &placeholder, "Ele32WPTightPass/F");
  tree->Branch("Ele35WPTightPass", &placeholder, "Ele35WPTightPass/F");

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
  tree->Branch("dZ_1", &in->ePVDZ, "dZ_1/F");
  tree->Branch("d0_1", &in->ePVDXY, "d0_1/F");
  tree->Branch("q_1", &in->eCharge, "q_1/F");
  tree->Branch("iso_1", &in->eIsoDB03, "iso_1/F");
  tree->Branch("NoisoID80_1", &in->eMVANoisoWP80, "NoisoID80_1/F");
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
  tree->Branch("l2_decayMode", &in->tDecayMode, "l2_decayMode/F");

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

  tree->Branch("rho", &in->rho, "rho/F");
  tree->Branch("metcov00", &in->metcov00, "metcov00/F");
  tree->Branch("metcov01", &in->metcov01, "metcov01/F");
  tree->Branch("metcov10", &in->metcov10, "metcov10/F");
  tree->Branch("metcov11", &in->metcov11, "metcov11/F");
  tree->Branch("metcov00_v2", &placeholder, "metcov00_v2/F");
  tree->Branch("metcov01_v2", &placeholder, "metcov01_v2/F");
  tree->Branch("metcov10_v2", &placeholder, "metcov10_v2/F");
  tree->Branch("metcov11_v2", &placeholder, "metcov11_v2/F");
  tree->Branch("NUP", &in->NUP, "NUP/F");
  tree->Branch("genM", &in->genM, "genM/F");
  tree->Branch("genpT", &in->genpT, "genpT/F");
  tree->Branch("genEta", &in->genEta, "genEta/F");
  tree->Branch("numGenJets", &in->numGenJets, "numGenJets/F");
  tree->Branch("npu", &in->nTruePU, "npu/F");
  tree->Branch("npv", &in->nvtx, "npv/F");
  tree->Branch("eMVANoisoWP80", &in->eMVANoisoWP80, "eMVANoisoWP80/F");
  tree->Branch("genweight", &in->GenWeight, "genweight/F");
  tree->Branch("metSig", &in->metSig, "metSig/F");
  tree->Branch("Rivet_higgsPt", &in->Rivet_higgsPt, "Rivet_higgsPt/F");
  tree->Branch("Rivet_nJets30", &in->Rivet_nJets30, "Rivet_nJets30/F");

  tree->Branch("jpt_1", &in->j1pt, "jpt_1/F");
  tree->Branch("jeta_1", &in->j1eta, "jeta_1/F");
  tree->Branch("jphi_1", &in->j1phi, "jphi_1/F");
  tree->Branch("jcsv_1", &in->j1csv, "jcsv_1/F");
  tree->Branch("jpt_2", &in->j2pt, "jpt_2/F");
  tree->Branch("jeta_2", &in->j2eta, "jeta_2/F");
  tree->Branch("jphi_2", &in->j2phi, "jphi_2/F");
  tree->Branch("jcsv_2", &in->j2csv, "jcsv_2/F");
  tree->Branch("bpt_1", &in->jb1pt, "bpt_1/F");
  tree->Branch("beta_1", &in->jb1eta, "beta_1/F");
  tree->Branch("bphi_1", &in->jb1phi, "bphi_1/F");
  tree->Branch("bcsv_1", &in->jb1csv, "bcsv_1/F");
  tree->Branch("bflavor_1", &in->jb1hadronflavor, "bflavor_1/F");
  tree->Branch("bpt_2", &in->jb2pt, "bpt_2/F");
  tree->Branch("beta_2", &in->jb2eta, "beta_2/F");
  tree->Branch("bphi_2", &in->jb2phi, "bphi_2/F");
  tree->Branch("bcsv_2", &in->jb2csv, "bcsv_2/F");
  tree->Branch("bflavor_2", &in->jb2hadronflavor, "bflavor_2/F");

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
  tree->Branch("eGenCharge", &in->eGenCharge, "eGenCharge/F");
  tree->Branch("eGenDirectPromptTauDecay", &in->eGenDirectPromptTauDecay, "eGenDirectPromptTauDecay/F");
  tree->Branch("eGenEnergy", &in->eGenEnergy, "eGenEnergy/F");
  tree->Branch("eGenEta", &in->eGenEta, "eGenEta/F");
  tree->Branch("eGenIsPrompt", &in->eGenIsPrompt, "eGenIsPrompt/F");
  tree->Branch("eGenMotherPdgId", &in->eGenMotherPdgId, "eGenMotherPdgId/F");
  tree->Branch("eGenParticle", &in->eGenParticle, "eGenParticle/F");
  tree->Branch("eGenPdgId", &in->eGenPdgId, "eGenPdgId/F");
  tree->Branch("eGenPhi", &in->eGenPhi, "eGenPhi/F");
  tree->Branch("eGenPrompt", &in->eGenPrompt, "eGenPrompt/F");
  tree->Branch("eGenPromptTauDecay", &in->eGenPromptTauDecay, "eGenPromptTauDecay/F");
  tree->Branch("eGenPt", &in->eGenPt, "eGenPt/F");
  tree->Branch("eGenTauDecay", &in->eGenTauDecay, "eGenTauDecay/F");
  tree->Branch("eGenVZ", &in->eGenVZ, "eGenVZ/F");
  tree->Branch("eGenVtxPVMatch", &in->eGenVtxPVMatch, "eGenVtxPVMatch/F");

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
  // tree->Branch("type1_pfMet_shiftedPt_UnclusteredEnUp", &type1_pfMet_shiftedPt_UnclusteredEnUp, "type1_pfMet_shiftedPt_UnclusteredEnUp/F");
  // tree->Branch("type1_pfMet_shiftedPhi_UnclusteredEnUp", &type1_pfMet_shiftedPhi_UnclusteredEnUp, "type1_pfMet_shiftedPhi_UnclusteredEnUp/F");
  // tree->Branch("type1_pfMet_shiftedPt_UnclusteredEnDown", &type1_pfMet_shiftedPt_UnclusteredEnDown, "type1_pfMet_shiftedPt_UnclusteredEnDown/F");
  // tree->Branch("type1_pfMet_shiftedPhi_UnclusteredEnDown", &type1_pfMet_shiftedPhi_UnclusteredEnDown, "type1_pfMet_shiftedPhi_UnclusteredEnDown/F");
  // tree->Branch("type1_pfMet_shiftedPt_JetEnUp", &type1_pfMet_shiftedPt_JetEnUp, "type1_pfMet_shiftedPt_JetEnUp/F");
  // tree->Branch("type1_pfMet_shiftedPhi_JetEnUp", &type1_pfMet_shiftedPhi_JetEnUp, "type1_pfMet_shiftedPhi_JetEnUp/F");
  // tree->Branch("type1_pfMet_shiftedPt_JetEnDown", &type1_pfMet_shiftedPt_JetEnDown, "type1_pfMet_shiftedPt_JetEnDown/F");
  // tree->Branch("type1_pfMet_shiftedPhi_JetEnDown", &type1_pfMet_shiftedPhi_JetEnDown, "type1_pfMet_shiftedPhi_JetEnDown/F");

  // tree->Branch("jetVeto30WoNoisyJets_JetEta0to3Down", &jetVeto30WoNoisyJets_JetEta0to3Down, "jetVeto30WoNoisyJets_JetEta0to3Down/F");
  // tree->Branch("jetVeto30WoNoisyJets_JetEta0to3Up", &jetVeto30WoNoisyJets_JetEta0to3Up);
  // tree->Branch("jetVeto30WoNoisyJets_JetEta0to5Down", &jetVeto30WoNoisyJets_JetEta0to5Down);
  // tree->Branch("jetVeto30WoNoisyJets_JetEta0to5Up", &jetVeto30WoNoisyJets_JetEta0to5Up);
  // tree->Branch("jetVeto30WoNoisyJets_JetEta3to5Down", &jetVeto30WoNoisyJets_JetEta3to5Down);
  // tree->Branch("jetVeto30WoNoisyJets_JetEta3to5Up", &jetVeto30WoNoisyJets_JetEta3to5Up);
  // tree->Branch("jetVeto30WoNoisyJets_JetRelativeBalDownWoNoisyJets", &jetVeto30WoNoisyJets_JetRelativeBalDownWoNoisyJets);
  // tree->Branch("jetVeto30WoNoisyJets_JetRelativeBalUpWoNoisyJets", &jetVeto30WoNoisyJets_JetRelativeBalUpWoNoisyJets);
  // tree->Branch("jetVeto30WoNoisyJets_JetRelativeSampleDown", &jetVeto30WoNoisyJets_JetRelativeSampleDown);
  // tree->Branch("jetVeto30WoNoisyJets_JetRelativeSampleUp", &jetVeto30WoNoisyJets_JetRelativeSampleUp);
  // tree->Branch("jetVeto30WoNoisyJets_JetTotalDown", &jetVeto30WoNoisyJets_JetTotalDown);
  // tree->Branch("jetVeto30WoNoisyJets_JetTotalUp", &jetVeto30WoNoisyJets_JetTotalUp);
  // tree->Branch("jetVeto30_JetAbsoluteFlavMapDown", &jetVeto30_JetAbsoluteFlavMapDown);
  // tree->Branch("jetVeto30_JetAbsoluteFlavMapUp", &jetVeto30_JetAbsoluteFlavMapUp);
  // tree->Branch("jetVeto30_JetAbsoluteMPFBiasDown", &jetVeto30_JetAbsoluteMPFBiasDown);
  // tree->Branch("jetVeto30_JetAbsoluteMPFBiasUp", &jetVeto30_JetAbsoluteMPFBiasUp);
  // tree->Branch("jetVeto30_JetAbsoluteScaleDown", &jetVeto30_JetAbsoluteScaleDown);
  // tree->Branch("jetVeto30_JetAbsoluteScaleUp", &jetVeto30_JetAbsoluteScaleUp);
  // tree->Branch("jetVeto30_JetAbsoluteStatDown", &jetVeto30_JetAbsoluteStatDown);
  // tree->Branch("jetVeto30_JetAbsoluteStatUp", &jetVeto30_JetAbsoluteStatUp);
  // tree->Branch("jetVeto30_JetClosureDown", &jetVeto30_JetClosureDown);
  // tree->Branch("jetVeto30_JetClosureUp", &jetVeto30_JetClosureUp);
  // tree->Branch("jetVeto30_JetEnDown", &jetVeto30_JetEnDown);
  // tree->Branch("jetVeto30_JetFlavorQCDDown", &jetVeto30_JetFlavorQCDDown);
  // tree->Branch("jetVeto30_JetFlavorQCDUp", &jetVeto30_JetFlavorQCDUp);
  // tree->Branch("jetVeto30_JetFragmentationDown", &jetVeto30_JetFragmentationDown);
  // tree->Branch("jetVeto30_JetFragmentationUp", &jetVeto30_JetFragmentationUp);
  // tree->Branch("jetVeto30_JetPileUpDataMCDown", &jetVeto30_JetPileUpDataMCDown);
  // tree->Branch("jetVeto30_JetPileUpDataMCUp", &jetVeto30_JetPileUpDataMCUp);
  // tree->Branch("jetVeto30_JetPileUpPtBBDown", &jetVeto30_JetPileUpPtBBDown);
  // tree->Branch("jetVeto30_JetPileUpPtBBUp", &jetVeto30_JetPileUpPtBBUp);
  // tree->Branch("jetVeto30_JetPileUpPtEC1Down", &jetVeto30_JetPileUpPtEC1Down);
  // tree->Branch("jetVeto30_JetPileUpPtEC1Up", &jetVeto30_JetPileUpPtEC1Up);
  // tree->Branch("jetVeto30_JetPileUpPtEC2Down", &jetVeto30_JetPileUpPtEC2Down);
  // tree->Branch("jetVeto30_JetPileUpPtEC2Up", &jetVeto30_JetPileUpPtEC2Up);
  // tree->Branch("jetVeto30_JetPileUpPtHFDown", &jetVeto30_JetPileUpPtHFDown);
  // tree->Branch("jetVeto30_JetPileUpPtHFUp", &jetVeto30_JetPileUpPtHFUp);
  // tree->Branch("jetVeto30_JetPileUpPtRefDown", &jetVeto30_JetPileUpPtRefDown);
  // tree->Branch("jetVeto30_JetPileUpPtRefUp", &jetVeto30_JetPileUpPtRefUp);
  // tree->Branch("jetVeto30_JetRelativeBalDown", &jetVeto30_JetRelativeBalDown);
  // tree->Branch("jetVeto30_JetRelativeBalUp", &jetVeto30_JetRelativeBalUp);
  // tree->Branch("jetVeto30_JetRelativeFSRDown", &jetVeto30_JetRelativeFSRDown);
  // tree->Branch("jetVeto30_JetRelativeFSRUp", &jetVeto30_JetRelativeFSRUp);
  // tree->Branch("jetVeto30_JetRelativeJEREC1Down", &jetVeto30_JetRelativeJEREC1Down);
  // tree->Branch("jetVeto30_JetRelativeJEREC1Up", &jetVeto30_JetRelativeJEREC1Up);
  // tree->Branch("jetVeto30_JetRelativeJEREC2Down", &jetVeto30_JetRelativeJEREC2Down);
  // tree->Branch("jetVeto30_JetRelativeJEREC2Up", &jetVeto30_JetRelativeJEREC2Up);
  // tree->Branch("jetVeto30_JetRelativeJERHFDown", &jetVeto30_JetRelativeJERHFDown);
  // tree->Branch("jetVeto30_JetRelativeJERHFUp", &jetVeto30_JetRelativeJERHFUp);
  // tree->Branch("jetVeto30_JetRelativePtBBDown", &jetVeto30_JetRelativePtBBDown);
  // tree->Branch("jetVeto30_JetRelativePtBBUp", &jetVeto30_JetRelativePtBBUp);
  // tree->Branch("jetVeto30_JetRelativePtEC1Down", &jetVeto30_JetRelativePtEC1Down);
  // tree->Branch("jetVeto30_JetRelativePtEC1Up", &jetVeto30_JetRelativePtEC1Up);
  // tree->Branch("jetVeto30_JetRelativePtEC2Down", &jetVeto30_JetRelativePtEC2Down);
  // tree->Branch("jetVeto30_JetRelativePtEC2Up", &jetVeto30_JetRelativePtEC2Up);
  // tree->Branch("jetVeto30_JetRelativePtHFDown", &jetVeto30_JetRelativePtHFDown);
  // tree->Branch("jetVeto30_JetRelativePtHFUp", &jetVeto30_JetRelativePtHFUp);
  // tree->Branch("jetVeto30_JetRelativeSampleDown", &jetVeto30_JetRelativeSampleDown);
  // tree->Branch("jetVeto30_JetRelativeSampleUp", &jetVeto30_JetRelativeSampleUp);
  // tree->Branch("jetVeto30_JetRelativeStatECDown", &jetVeto30_JetRelativeStatECDown);
  // tree->Branch("jetVeto30_JetRelativeStatECUp", &jetVeto30_JetRelativeStatECUp);
  // tree->Branch("jetVeto30_JetRelativeStatFSRDown", &jetVeto30_JetRelativeStatFSRDown);
  // tree->Branch("jetVeto30_JetRelativeStatFSRUp", &jetVeto30_JetRelativeStatFSRUp);
  // tree->Branch("jetVeto30_JetRelativeStatHFDown", &jetVeto30_JetRelativeStatHFDown);
  // tree->Branch("jetVeto30_JetRelativeStatHFUp", &jetVeto30_JetRelativeStatHFUp);
  // tree->Branch("jetVeto30_JetSinglePionECALDown", &jetVeto30_JetSinglePionECALDown);
  // tree->Branch("jetVeto30_JetSinglePionECALUp", &jetVeto30_JetSinglePionECALUp);
  // tree->Branch("jetVeto30_JetSinglePionHCALDown", &jetVeto30_JetSinglePionHCALDown);
  // tree->Branch("jetVeto30_JetSinglePionHCALUp", &jetVeto30_JetSinglePionHCALUp);
  // tree->Branch("jetVeto30_JetTimePtEtaDown", &jetVeto30_JetTimePtEtaDown);
  // tree->Branch("jetVeto30_JetTimePtEtaUp", &jetVeto30_JetTimePtEtaUp);
  // tree->Branch("jetVeto30_JetTotalDown", &jetVeto30_JetTotalDown);
  // tree->Branch("jetVeto30_JetTotalUp", &jetVeto30_JetTotalUp);

  // tree->Branch("vbfMassWoNoisyJets_JetEta0to3Down", &vbfMassWoNoisyJets_JetEta0to3Down);
  // tree->Branch("vbfMassWoNoisyJets_JetEta0to3Up", &vbfMassWoNoisyJets_JetEta0to3Up);
  // tree->Branch("vbfMassWoNoisyJets_JetEta0to5Down", &vbfMassWoNoisyJets_JetEta0to5Down);
  // tree->Branch("vbfMassWoNoisyJets_JetEta0to5Up", &vbfMassWoNoisyJets_JetEta0to5Up);
  // tree->Branch("vbfMassWoNoisyJets_JetEta3to5Down", &vbfMassWoNoisyJets_JetEta3to5Down);
  // tree->Branch("vbfMassWoNoisyJets_JetEta3to5Up", &vbfMassWoNoisyJets_JetEta3to5Up);
  // tree->Branch("vbfMassWoNoisyJets_JetRelativeSampleDown", &vbfMassWoNoisyJets_JetRelativeSampleDown);
  // tree->Branch("vbfMassWoNoisyJets_JetRelativeSampleUp", &vbfMassWoNoisyJets_JetRelativeSampleUp);
  // tree->Branch("vbfMassWoNoisyJets_JetTotalDown", &vbfMassWoNoisyJets_JetTotalDown);
  // tree->Branch("vbfMassWoNoisyJets_JetTotalUp", &vbfMassWoNoisyJets_JetTotalUp);
  // tree->Branch("vbfMass_JetAbsoluteFlavMapDown", &vbfMass_JetAbsoluteFlavMapDown);
  // tree->Branch("vbfMass_JetAbsoluteFlavMapUp", &vbfMass_JetAbsoluteFlavMapUp);
  // tree->Branch("vbfMass_JetAbsoluteMPFBiasDown", &vbfMass_JetAbsoluteMPFBiasDown);
  // tree->Branch("vbfMass_JetAbsoluteMPFBiasUp", &vbfMass_JetAbsoluteMPFBiasUp);
  // tree->Branch("vbfMass_JetAbsoluteScaleDown", &vbfMass_JetAbsoluteScaleDown);
  // tree->Branch("vbfMass_JetAbsoluteScaleUp", &vbfMass_JetAbsoluteScaleUp);
  // tree->Branch("vbfMass_JetAbsoluteStatDown", &vbfMass_JetAbsoluteStatDown);
  // tree->Branch("vbfMass_JetAbsoluteStatUp", &vbfMass_JetAbsoluteStatUp);
  // tree->Branch("vbfMass_JetClosureDown", &vbfMass_JetClosureDown);
  // tree->Branch("vbfMass_JetClosureUp", &vbfMass_JetClosureUp);
  // tree->Branch("vbfMass_JetFlavorQCDDown", &vbfMass_JetFlavorQCDDown);
  // tree->Branch("vbfMass_JetFlavorQCDUp", &vbfMass_JetFlavorQCDUp);
  // tree->Branch("vbfMass_JetFragmentationDown", &vbfMass_JetFragmentationDown);
  // tree->Branch("vbfMass_JetFragmentationUp", &vbfMass_JetFragmentationUp);
  // tree->Branch("vbfMass_JetPileUpDataMCDown", &vbfMass_JetPileUpDataMCDown);
  // tree->Branch("vbfMass_JetPileUpDataMCUp", &vbfMass_JetPileUpDataMCUp);
  // tree->Branch("vbfMass_JetPileUpPtBBDown", &vbfMass_JetPileUpPtBBDown);
  // tree->Branch("vbfMass_JetPileUpPtBBUp", &vbfMass_JetPileUpPtBBUp);
  // tree->Branch("vbfMass_JetPileUpPtEC1Down", &vbfMass_JetPileUpPtEC1Down);
  // tree->Branch("vbfMass_JetPileUpPtEC1Up", &vbfMass_JetPileUpPtEC1Up);
  // tree->Branch("vbfMass_JetPileUpPtEC2Down", &vbfMass_JetPileUpPtEC2Down);
  // tree->Branch("vbfMass_JetPileUpPtEC2Up", &vbfMass_JetPileUpPtEC2Up);
  // tree->Branch("vbfMass_JetPileUpPtHFDown", &vbfMass_JetPileUpPtHFDown);
  // tree->Branch("vbfMass_JetPileUpPtHFUp", &vbfMass_JetPileUpPtHFUp);
  // tree->Branch("vbfMass_JetPileUpPtRefDown", &vbfMass_JetPileUpPtRefDown);
  // tree->Branch("vbfMass_JetPileUpPtRefUp", &vbfMass_JetPileUpPtRefUp);
  // tree->Branch("vbfMass_JetRelativeBalDown", &vbfMass_JetRelativeBalDown);
  // tree->Branch("vbfMass_JetRelativeBalUp", &vbfMass_JetRelativeBalUp);
  // tree->Branch("vbfMass_JetRelativeFSRDown", &vbfMass_JetRelativeFSRDown);
  // tree->Branch("vbfMass_JetRelativeFSRUp", &vbfMass_JetRelativeFSRUp);
  // tree->Branch("vbfMass_JetRelativeJEREC1Down", &vbfMass_JetRelativeJEREC1Down);
  // tree->Branch("vbfMass_JetRelativeJEREC1Up", &vbfMass_JetRelativeJEREC1Up);
  // tree->Branch("vbfMass_JetRelativeJEREC2Down", &vbfMass_JetRelativeJEREC2Down);
  // tree->Branch("vbfMass_JetRelativeJEREC2Up", &vbfMass_JetRelativeJEREC2Up);
  // tree->Branch("vbfMass_JetRelativeJERHFDown", &vbfMass_JetRelativeJERHFDown);
  // tree->Branch("vbfMass_JetRelativeJERHFUp", &vbfMass_JetRelativeJERHFUp);
  // tree->Branch("vbfMass_JetRelativePtBBDown", &vbfMass_JetRelativePtBBDown);
  // tree->Branch("vbfMass_JetRelativePtBBUp", &vbfMass_JetRelativePtBBUp);
  // tree->Branch("vbfMass_JetRelativePtEC1Down", &vbfMass_JetRelativePtEC1Down);
  // tree->Branch("vbfMass_JetRelativePtEC1Up", &vbfMass_JetRelativePtEC1Up);
  // tree->Branch("vbfMass_JetRelativePtEC2Down", &vbfMass_JetRelativePtEC2Down);
  // tree->Branch("vbfMass_JetRelativePtEC2Up", &vbfMass_JetRelativePtEC2Up);
  // tree->Branch("vbfMass_JetRelativePtHFDown", &vbfMass_JetRelativePtHFDown);
  // tree->Branch("vbfMass_JetRelativePtHFUp", &vbfMass_JetRelativePtHFUp);
  // tree->Branch("vbfMass_JetRelativeSampleDown", &vbfMass_JetRelativeSampleDown);
  // tree->Branch("vbfMass_JetRelativeSampleUp", &vbfMass_JetRelativeSampleUp);
  // tree->Branch("vbfMass_JetRelativeStatECDown", &vbfMass_JetRelativeStatECDown);
  // tree->Branch("vbfMass_JetRelativeStatECUp", &vbfMass_JetRelativeStatECUp);
  // tree->Branch("vbfMass_JetRelativeStatFSRDown", &vbfMass_JetRelativeStatFSRDown);
  // tree->Branch("vbfMass_JetRelativeStatFSRUp", &vbfMass_JetRelativeStatFSRUp);
  // tree->Branch("vbfMass_JetRelativeStatHFDown", &vbfMass_JetRelativeStatHFDown);
  // tree->Branch("vbfMass_JetRelativeStatHFUp", &vbfMass_JetRelativeStatHFUp);
  // tree->Branch("vbfMass_JetSinglePionECALDown", &vbfMass_JetSinglePionECALDown);
  // tree->Branch("vbfMass_JetSinglePionECALUp", &vbfMass_JetSinglePionECALUp);
  // tree->Branch("vbfMass_JetSinglePionHCALDown", &vbfMass_JetSinglePionHCALDown);
  // tree->Branch("vbfMass_JetSinglePionHCALUp", &vbfMass_JetSinglePionHCALUp);
  // tree->Branch("vbfMass_JetTimePtEtaDown", &vbfMass_JetTimePtEtaDown);
  // tree->Branch("vbfMass_JetTimePtEtaUp", &vbfMass_JetTimePtEtaUp);
  // tree->Branch("vbfMass_JetTotalDown", &vbfMass_JetTotalDown);
  // tree->Branch("vbfMass_JetTotalUp", &vbfMass_JetTotalUp);

  // 2016 placeholders
  tree->Branch("amcatNLO_weight", &placeholder, "amcatNLO_weight/F");
  tree->Branch("eMatchesSingleE25Tight", &placeholder, "eMatchesSingleE25Tight/F");
  tree->Branch("eMatchesEle25TightFilter", &placeholder, "eMatchesEle25TightFilter/F");
  tree->Branch("singleE25eta2p1TightPass", &placeholder, "singleE25eta2p1TightPass/F");
  tree->Branch("againstElectronTightMVA6_2", &placeholder, "againstElectronTightMVA6_2/F");
  tree->Branch("againstElectronVLooseMVA6_2", &placeholder, "againstElectronVLooseMVA6_2/F");
  tree->Branch("againstMuonTight3_2", &placeholder, "againstMuonTight3_2/F");
  tree->Branch("againstMuonLoose3_2", &placeholder, "againstMuonLoose3_2/F");
  tree->Branch("byVLooseIsolationMVArun2v1DBoldDMwLT_2", &placeholder, "byVLooseIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("decayModeFindingNewDMs_2", &placeholder, "decayModeFindingNewDMs_2/F");
}
