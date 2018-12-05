// Copyright 2018 Tyler Mitchell

#ifndef ROOT_SRC_MUTAU_TREE_2017_H_
#define ROOT_SRC_MUTAU_TREE_2017_H_

#include <cmath>
#include <iostream>
#include <vector>
#include <utility>
#include "RecoilCorrector.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TMath.h"

class mutau_tree {
 private:
  TTree *tree, *original;
  bool isMC, isEmbed;
  std::vector<Int_t> good_events;
  TLorentzVector mu, tau, pfMET, mvaMET, puppiMET, rawPfMET;

 public:
  // Member variables
  Int_t Run, Lumi, recoil;

  // Selection Variables.
  Float_t mIsGlobal, mNormalizedChi2, mChi2LocalPosition, mTrkKink, mValidFraction, mSegmentCompatibility, tRerunMVArun2v2DBoldDMwLTVVLoose, tDecayModeFinding;

  // Intermediate Variables. (Missing trigweight_2, idisoweight_2)
  Float_t muVetoZTTp001dxyzR0, eVetoZTTp001dxyzR0, dielectronVeto, dimuonVeto, mMatchesIsoMu20Tau27Filter, mMatchesIsoMu20Tau27Path, tMatchesIsoMu20Tau27Filter,
      tMatchesIsoMu20Tau27Path, mMatchesIsoMu24Filter, mMatchesIsoMu24Path, mMatchesIsoMu27Filter, mMatchesIsoMu27Path, Mu20Tau27Pass, IsoMu27Pass,
      IsoMu24Pass, Flag_BadChargedCandidateFilter, Flag_BadPFMuonFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_HBHENoiseFilter, Flag_HBHENoiseIsoFilter,
      Flag_badMuons, Flag_duplicateMuons, Flag_ecalBadCalibFilter, Flag_eeBadScFilter, Flag_globalSuperTightHalo2016Filter, Flag_globalTightHalo2016Filter,
      Flag_goodVertices, mPt, mPhi, mEta, mMass, mCharge, mPVDXY, mPVDZ, mRelPFIsoDBDefaultR04, mZTTGenMatching, tPt, tEta, tPhi, tMass, tCharge, tPVDXY, tPVDZ, tZTTGenMatching,
      tDecayMode, type1_pfMetEt, type1_pfMetPhi, raw_pfMetEt, raw_pfMetPhi, puppiMetEt,
      puppiMetPhi, tAgainstElectronVLooseMVA6, tAgainstElectronLooseMVA6, tAgainstElectronMediumMVA6, tAgainstElectronTightMVA6, tAgainstElectronVTightMVA6,
      tAgainstMuonLoose3, tAgainstMuonTight3, tRerunMVArun2v2DBoldDMwLTraw, vbfMassWoNoisyJets, jetVeto20WoNoisyJets, jetVeto30WoNoisyJets, bjetDeepCSVVeto20MediumWoNoisyJets,
      j1ptWoNoisyJets, j1etaWoNoisyJets, j1phiWoNoisyJets, j1csvWoNoisyJets, j2ptWoNoisyJets, j2etaWoNoisyJets, j2phiWoNoisyJets, j2csvWoNoisyJets,
      jb1ptWoNoisyJets, jb1etaWoNoisyJets, jb1phiWoNoisyJets, jb1csvWoNoisyJets, jb1hadronflavorWoNoisyJets, jb2ptWoNoisyJets, jb2etaWoNoisyJets,
      jb2phiWoNoisyJets, jb2csvWoNoisyJets, jb2hadronflavorWoNoisyJets, mPFIDLoose, mPFIDMedium, mPFIDTight, tRerunMVArun2v2DBoldDMwLTLoose,
      tRerunMVArun2v2DBoldDMwLTMedium, tRerunMVArun2v2DBoldDMwLTTight, nvtx, nTruePU, genpX, genpY, vispX, vispY, pfmetcorr_ex, pfmetcorr_ey;

  // Output Variables.
  ULong64_t evt;
  UInt_t run, lumi;
  Float_t dilepton_veto, extraelec_veto, extramuon_veto, trg_singleelectron, trg_singlemuon, trg_singletau, trg_muonelectron, trg_mutaucross, trg_doubletau,
      flagFilter, pt_1, phi_1, eta_1, m_1, q_1, d0_1, dZ_1, mt_1, pfmt_1, puppimt_1, iso_1, gen_match_1, againstElectronLooseMVA6_1, againstElectronMediumMVA6_1,
      againstElectronTightMVA6_1, againstElectronVLooseMVA6_1, againstElectronVTightMVA6_1, againstMuonLoose3_1, againstMuonTight3_1, byIsolationMVA3oldDMwLTraw_1,
      trigweight_1, idisoweight_1, pt_2, phi_2, eta_2, m_2, q_2, d0_2, dZ_2, mt_2, iso_2, gen_match_2, againstElectronLooseMVA6_2, againstElectronMediumMVA6_2,
      againstElectronTightMVA6_2, againstElectronVLooseMVA6_2, againstElectronVTightMVA6_2, againstMuonLoose3_2, againstMuonTight3_2, byIsolationMVA3oldDMwLTraw_2,
      trigweight_2, idisoweight_2, pt_tt, mt_tot, m_vis, m_sv, mt_sv, met, metphi, puppimet, puppimetphi, pzetavis, pzetamiss, pfpzetamiss, puppipzetamiss,
      metcov00, metcov01, metcov10, metcov11, mjj, jdeta, njetingap, njetingap20, jdphi, dijetpt, dijetphi, ptvis, nbtag, njets, njetspt20, jpt_1, jeta_1,
      jphi_1, jcsv_1, jpt_2, jeta_2, jphi_2, jcsv_2, bpt_1, beta_1, bphi_1, bcsv_1, bpt_2, beta_2, bphi_2, bcsv_2, puweight, NUP, weight, id_m_loose_1,
      id_m_medium_1, id_m_tight_1, id_m_loose_2, id_m_medium_2, id_m_tight_2, pt_sv, eta_sv, phi_sv, met_sv, jpfid_1, jpuid_1, jpfid_2, jpuid_2, bpfid_1,
      bpuid_1, bpfid_2, bpuid_2, npv, npu, rho;



  // Member functions
  mutau_tree(TTree *orig, TTree *itree, bool isMC, bool isEmbed, Int_t rec);
  virtual ~mutau_tree() {}
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
mutau_tree::mutau_tree(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec) :
tree(itree),
original(Original),
isMC(IsMC),
isEmbed(IsEmbed),
recoil(rec) {
  original->SetBranchAddress("mPt", &mPt);
  original->SetBranchAddress("mEta", &mEta);
  original->SetBranchAddress("mPhi", &mPhi);
  original->SetBranchAddress("mMass", &mMass);
  original->SetBranchAddress("mPVDZ", &mPVDZ);
  original->SetBranchAddress("mPVDXY", &mPVDXY);
  original->SetBranchAddress("tPt", &tPt);
  original->SetBranchAddress("tEta", &tEta);
  original->SetBranchAddress("tPhi", &tPhi);
  original->SetBranchAddress("tMass", &tMass);
  original->SetBranchAddress("tCharge", &tCharge);
  original->SetBranchAddress("tPVDZ", &tPVDZ);
  original->SetBranchAddress("tPVDXY", &tPVDXY);
  original->SetBranchAddress("mIsGlobal", &mIsGlobal);
  original->SetBranchAddress("mNormalizedChi2", &mNormalizedChi2);
  original->SetBranchAddress("mChi2LocalPosition", &mChi2LocalPosition);
  original->SetBranchAddress("mTrkKink", &mTrkKink);
  original->SetBranchAddress("mPFIDLoose", &mPFIDLoose);
  original->SetBranchAddress("mPFIDMedium", &mPFIDMedium);
  original->SetBranchAddress("mValidFraction", &mValidFraction);
  original->SetBranchAddress("mSegmentCompatibility", &mSegmentCompatibility);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVVLoose", &tRerunMVArun2v2DBoldDMwLTVVLoose);

  original->SetBranchAddress("tDecayMode", &tDecayMode);
  original->SetBranchAddress("tDecayModeFinding", &tDecayModeFinding);
  original->SetBranchAddress("tZTTGenMatching", &tZTTGenMatching);
  original->SetBranchAddress("tMatchesIsoMu20Tau27Filter", &tMatchesIsoMu20Tau27Filter);
  original->SetBranchAddress("tMatchesIsoMu20Tau27Path", &tMatchesIsoMu20Tau27Path);
  original->SetBranchAddress("mMatchesIsoMu20Tau27Filter", &mMatchesIsoMu20Tau27Filter);
  original->SetBranchAddress("mMatchesIsoMu20Tau27Path", &mMatchesIsoMu20Tau27Path);
  original->SetBranchAddress("mMatchesIsoMu24Filter", &mMatchesIsoMu24Filter);
  original->SetBranchAddress("mMatchesIsoMu24Path", &mMatchesIsoMu24Path);
  original->SetBranchAddress("mMatchesIsoMu27Filter", &mMatchesIsoMu27Filter);
  original->SetBranchAddress("mMatchesIsoMu27Path", &mMatchesIsoMu27Path);
  original->SetBranchAddress("Mu20Tau27Pass", &Mu20Tau27Pass);
  original->SetBranchAddress("IsoMu27Pass", &IsoMu27Pass);
  original->SetBranchAddress("IsoMu24Pass", &IsoMu24Pass);
}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void mutau_tree::do_skimming(TH1F* cutflow) {
  // declare variables for sorting
  ULong64_t evt_now(0);
  ULong64_t evt_before(1);
  int best_evt(-1);
  std::pair<float, float> muCandidate, tauCandidate;

  Int_t nevt = (Int_t)original->GetEntries();
  for (auto ievt = 0; ievt < nevt; ievt++) {
    original->GetEntry(ievt);
    evt_now = evt;

    // TLorentzVector ele, tau;
    mu.SetPtEtaPhiM(mPt, mEta, mPhi, mMass);
    tau.SetPtEtaPhiM(tPt, tEta, tPhi, tMass);

    // apply TES
    if (isMC && !isEmbed) {
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          tau *= 1.007;
        } else if (tDecayMode == 1) {
          tau *= 0.998;
        } else if (tDecayMode == 10) {
          tau *= 1.001;
        }
      } else if (tZTTGenMatching == 1 || tZTTGenMatching == 3) {
        if (tDecayMode == 0) {
          tau *= 1.003;
        } else if (tDecayMode == 1) {
          tau *= 1.036;
        }
      }
    } else if (isEmbed) {
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          tau *= 0.975;
        } else if (tDecayMode == 1) {
          tau *= (0.975*1.051);
        } else if (tDecayMode == 10) {
          tau *= (0.975*0.975*0.975);
        }
      }
    }

    cutflow->Fill(1., 1.);
    // apply event selection

    auto Mu24 = IsoMu24Pass && mMatchesIsoMu24Path && mMatchesIsoMu24Filter;
    auto Mu27 = IsoMu27Pass && mMatchesIsoMu27Path && mMatchesIsoMu27Filter;
    auto Cross = Mu20Tau27Pass && mMatchesIsoMu20Tau27Filter && mMatchesIsoMu20Tau27Path && tMatchesIsoMu20Tau27Filter && tMatchesIsoMu20Tau27Path;

    if (Mu24 || Mu27 || Cross) cutflow->Fill(2., 1.);
    else  continue;

    if (mPt > 21. && fabs(mEta) < 2.1 && fabs(mPVDZ) < 0.2 && fabs(mPVDXY) < 0.045) cutflow->Fill(3., 1.);  // electron kinematic selection
    else  continue;

    bool goodglob = mIsGlobal  && mNormalizedChi2 < 3  && mChi2LocalPosition < 12 && mTrkKink < 20;
    bool isMedium = mPFIDLoose && mValidFraction> 0.49 && mSegmentCompatibility > (goodglob ? 0.303 : 0.451);

    if (mPFIDMedium) cutflow->Fill(4., 1.);  // muon quality selection
    else  continue;

    if (isMC || isEmbed || isMedium) cutflow->Fill(5., 1.);  // muon quality selection
    else  continue;

    if (tau.Pt() > 20. && fabs(tau.Eta()) < 2.3 && fabs(tPVDZ) < 0.2) cutflow->Fill(6., 1.);  // tau kinematic selection
    else  continue;

    if (tRerunMVArun2v2DBoldDMwLTVVLoose && tDecayModeFinding > 0 && fabs(tCharge) < 2) cutflow->Fill(7., 1.);  // tau quality selection
    else  continue;

    if (mu.DeltaR(tau) > 0.5) cutflow->Fill(8., 1.);
    else  continue;

    // implement new sorting per
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
    if (evt_now != evt_before) {  // new event, save the tau candidates
      // since it is new event, do we have the best entry to save? If yes, save it!
      if ( best_evt > -1  )
        good_events.push_back(best_evt);

      //  this is a new event, so the first tau pair is the best! :)
      best_evt = ievt;
      muCandidate = std::make_pair(mPt, mRelPFIsoDBDefaultR04);
      tauCandidate  = std::make_pair(tPt,  tRerunMVArun2v2DBoldDMwLTraw);
    } else {  // not a new event
      std::pair<float, float> currEleCandidate(mPt, mRelPFIsoDBDefaultR04);
      std::pair<float, float> currTauCandidate(tPt, tRerunMVArun2v2DBoldDMwLTraw);

      // clause 1, select the pair that has most isolated tau lepton 1
      if (currEleCandidate.second - muCandidate.second  > 0.0001 ) best_evt = ievt;

      // check if the first tau is the same, and if so - move to clause 2
      if ( fabs(currEleCandidate.second - muCandidate.second)  <  0.0001 ) {
        // pick up  the pair with the highest pT of the first candidate
        if (currEleCandidate.first - muCandidate.first > 0.0001 ) best_evt = ievt;
        if ( fabs(currEleCandidate.first -muCandidate.first) < 0.0001 ) {
          // same pT, same iso, move to clause 3
          if (currTauCandidate.second - tauCandidate.second > 0.0001 ) best_evt = ievt;
          if ( fabs(currTauCandidate.second - tauCandidate.second) < 0.0001 ) {
            // same iso - pick the pair with the highest pT
            if ( currTauCandidate.first - tauCandidate.first  > 0.0001 ) best_evt = ievt;
          }  // tau2 has the same isolation
        }  // tau1 has the same pT
      }  // tau1 has the same isolation
    }  // not a new event
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
TTree* mutau_tree::fill_tree(RecoilCorrector recoilPFMetCorrector) {
  set_branches();  // get all the branches set up

  // loop through all events pasing skimming/sorting
  for (auto& ievt : good_events) {
    original->GetEntry(ievt);

    // convert from Float_t in FSA to Int_t for analyzer
    // evt
    run = Run;
    lumi = Lumi;
    extraelec_veto = eVetoZTTp001dxyzR0 > 1;
    extramuon_veto = muVetoZTTp001dxyzR0 > 0;
    dilepton_veto = dimuonVeto > 0;
    // triggers
    flagFilter = (Flag_BadChargedCandidateFilter == 0 && Flag_BadPFMuonFilter == 0 && Flag_EcalDeadCellTriggerPrimitiveFilter == 0
                  && Flag_HBHENoiseFilter == 0 && Flag_HBHENoiseIsoFilter == 0 && Flag_ecalBadCalibFilter == 0 && Flag_eeBadScFilter == 0
                  && Flag_globalSuperTightHalo2016Filter == 0 && Flag_goodVertices == 0);

    // TLorentzVector ele, tau;
    mu.SetPtEtaPhiM(mPt, mEta, mPhi, mMass);
    tau.SetPtEtaPhiM(tPt, tEta, tPhi, tMass);
    rawPfMET.SetPtEtaPhi(raw_pfMetEt, 0, raw_pfMetPhi, 0);
    pfMET.SetPtEtaPhiM(type1_pfMetEt, 0, type1_pfMetPhi, 0);
    pfmetcorr_ex = pfMET.Px();
    pfmetcorr_ey = pfMET.Py();
    mvaMET.SetPtEtaPhiM(0, 0, 0, 0);
    puppiMET.SetPtEtaPhiM(puppiMetEt, 0, puppiMetPhi, 0);

    if (recoil == 1) {
      recoilPFMetCorrector.CorrectByMeanResolution(
          rawPfMET.Px(),             // uncorrected type I pf met px (float)
          rawPfMET.Py(),             // uncorrected type I pf met py (float)
          genpX,                     // generator Z/W/Higgs px (float)
          genpY,                     // generator Z/W/Higgs py (float)
          vispX,                     // generator visible Z/W/Higgs px (float)
          vispY,                     // generator visible Z/W/Higgs py (float)
          jetVeto30WoNoisyJets + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,              // corrected type I pf met px (float)
          pfmetcorr_ey);             // corrected type I pf met py (float)
    } else if (recoil == 2) {
      recoilPFMetCorrector.CorrectByMeanResolution(
          rawPfMET.Px(),             // uncorrected type I pf met px (float)
          rawPfMET.Py(),             // uncorrected type I pf met py (float)
          genpX,                     // generator Z/W/Higgs px (float)
          genpY,                     // generator Z/W/Higgs py (float)
          vispX,                     // generator visible Z/W/Higgs px (float)
          vispY,                     // generator visible Z/W/Higgs py (float)
          jetVeto30WoNoisyJets + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,              // corrected type I pf met px (float)
          pfmetcorr_ey);             // corrected type I pf met py (float)
    }

    pfMET.SetPxPyPzE(pfMET.Px(), pfMET.Py(), 0, sqrt(pfMET.Px() * pfMET.Px() + pfMET.Py() * pfMET.Py()));
    mvaMET.SetPxPyPzE(mvaMET.Px(), mvaMET.Py(), 0, sqrt(mvaMET.Px() * mvaMET.Px() + mvaMET.Py() * mvaMET.Py()));
    puppiMET.SetPxPyPzE(puppiMET.Px(), puppiMET.Py(), 0, sqrt(puppiMET.Px() * puppiMET.Px() + puppiMET.Py() * puppiMET.Py()));

    if (isMC && !isEmbed) {
      // met correction due to tau energy scale
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          pfMET = pfMET + tau - 1.007*tau;
          tau *= 1.007;
        } else if (tDecayMode == 1) {
          pfMET = pfMET + tau - 0.998*tau;
          tau *= 0.998;
        } else if (tDecayMode == 10) {
          pfMET = pfMET + tau - 1.001*tau;
          tau *= 1.001;
        }
      } else if (tZTTGenMatching == 1 || tZTTGenMatching == 3) {
        if (tDecayMode == 0) {
          pfMET = pfMET+tau-1.003*tau;
          tau *= 1.003;
        } else if (tDecayMode == 1) {
          pfMET = pfMET+tau-1.036*tau;
          tau *= 1.036;
        }
      }
    } else if (isEmbed) {
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          pfMET = pfMET + tau - 0.975*tau;
          tau *= 0.975;
        } else if (tDecayMode == 1) {
          pfMET = pfMET + tau - 0.975*1.051*tau;
          tau *= 0.975*1.051;
        } else if (tDecayMode == 10) {
          pfMET = pfMET + tau - 0.975*0.975*0.975*tau;
          tau *= 0.975*0.975*0.975;
        }
      }
    }

    pt_1 = mu.Pt();
    eta_1 = mu.Eta();
    phi_1 = mu.Phi();
    m_1 = mu.M();
    q_1 = mCharge;
    d0_1 = mPVDXY;
    dZ_1 = mPVDZ;
    mt_1 = TMath::Sqrt(2 * mu.Pt() * mvaMET.Pt() * (1 - TMath::Cos(TMath::Abs(mu.Phi() - mvaMET.Phi()))));
    pfmt_1 = TMath::Sqrt(2 * mu.Pt() * pfMET.Pt() * (1 - TMath::Cos(TMath::Abs(mu.Phi() - pfMET.Phi()))));
    puppimt_1 = TMath::Sqrt(2 * mu.Pt() * puppiMET.Pt() * (1 - TMath::Cos(TMath::Abs(mu.Phi() - puppiMET.Phi()))));
    iso_1 = mRelPFIsoDBDefaultR04;
    gen_match_1 = mZTTGenMatching;
    againstElectronVLooseMVA6_1 = 0;
    againstElectronLooseMVA6_1 = 0;
    againstElectronMediumMVA6_1 = 0;
    againstElectronTightMVA6_1 = 0;
    againstElectronVTightMVA6_1 = 0;
    againstMuonLoose3_1 = 0;
    againstMuonTight3_1 = 0;
    byIsolationMVA3oldDMwLTraw_1 = 0;  // tau isolation
    // trigweight_1
    // idisoweight_1
    pt_2 = tau.Pt();
    eta_2 = tau.Eta();
    phi_2 = tau.Phi();
    m_2 = tau.M();
    q_2 = tCharge;
    d0_2 = tPVDXY;
    dZ_2 = tPVDZ;
    mt_2 = TMath::Sqrt(2 * tau.Pt() * mvaMET.Pt() * (1 - TMath::Cos(TMath::Abs(tau.Phi() - mvaMET.Phi()))));
    iso_2 = tRerunMVArun2v2DBoldDMwLTraw;
    gen_match_2 = tZTTGenMatching;
    againstElectronVLooseMVA6_2 = tAgainstElectronVLooseMVA6;
    againstElectronLooseMVA6_2 = tAgainstElectronLooseMVA6;
    againstElectronMediumMVA6_2 = tAgainstElectronMediumMVA6;
    againstElectronTightMVA6_2 = tAgainstElectronTightMVA6;
    againstElectronVTightMVA6_2 = tAgainstElectronVTightMVA6;
    againstMuonLoose3_2 = tAgainstMuonLoose3;
    againstMuonTight3_2 = tAgainstMuonTight3;
    byIsolationMVA3oldDMwLTraw_2 = tRerunMVArun2v2DBoldDMwLTraw;  // tau isolation
    pt_tt = (mu + tau + pfMET).Pt();
    mt_tot = TMath::Sqrt(2 * mu.Pt()  * mvaMET.Pt() * (1 - TMath::Cos(TMath::Abs(mu.Phi()  - mvaMET.Phi())))) +
             TMath::Sqrt(2 * tau.Pt() * mvaMET.Pt() * (1 - TMath::Cos(TMath::Abs(tau.Phi() - mvaMET.Phi())))) +
             TMath::Sqrt(2 * mu.Pt()  * tau.Pt()    * (1 - TMath::Cos(TMath::Abs(mu.Phi()  - tau.Phi()))));
    m_vis = (mu + tau).M();
    m_sv = 0;   // calculated later
    mt_sv = 0;  // calculated later
    met = type1_pfMetEt;
    metphi = type1_pfMetPhi;
    puppimet = puppiMetEt;
    puppimetphi = puppiMetPhi;
    pzetavis = ((mu.Px() + tau.Px()) * TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()) + TMath::Sin(mu.Phi()) + TMath::Cos(tau.Phi()) * (mu.Py() + tau.Py()))
              / (TMath::Sqrt(TMath::Power(TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()), 2) +
                            TMath::Power(TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()), 2)));
    pzetamiss = (pfMET.Px() * (TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi())) + pfMET.Py() * (TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi())))
              / (TMath::Sqrt(TMath::Power(TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()), 2) +
                            TMath::Power(TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()), 2)));
    pfpzetamiss = 0;
    puppipzetamiss = (puppiMET.Px() * (TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi())) + puppiMET.Py() * (TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi())))
              / (TMath::Sqrt(TMath::Power(TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()), 2) +
                            TMath::Power(TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()), 2)));
    // metcov00 (straight from tree)
    // metcov10 (straight from tree)
    // metcov10 (straight from tree)
    // metcov11 (straight from tree)

    TLorentzVector jet1, jet2;
    if (jetVeto20WoNoisyJets > 0 && j1ptWoNoisyJets > 0)
      jet1.SetPtEtaPhiM(j1ptWoNoisyJets, j1etaWoNoisyJets, j1phiWoNoisyJets, 0);
    if (jetVeto20WoNoisyJets > 1 && j2ptWoNoisyJets > 0)
      jet2.SetPtEtaPhiM(j2ptWoNoisyJets, j2etaWoNoisyJets, j2phiWoNoisyJets, 0);
    TLorentzVector dijet = jet1 + jet2;

    // need to add VBF system variables
    mjj = (jet1 + jet2).M();
    jdeta = TMath::Abs(jet1.Eta() - jet2.Eta());
    // njetingap
    // njetingap20
    jdphi = TMath::Abs(jet1.Phi() - jet2.Phi());
    dijetpt = (jet1 + jet2).Pt();
    dijetphi = (jet1 + jet2).Phi();
    ptvis = (mu + tau).Pt();
    nbtag = bjetDeepCSVVeto20MediumWoNoisyJets;
    njets = jetVeto30WoNoisyJets;
    njetspt20 = jetVeto20WoNoisyJets;
    jpt_1 = j1ptWoNoisyJets;
    jeta_1 = j1etaWoNoisyJets;
    jphi_1 = j1phiWoNoisyJets;
    jcsv_1 = j1csvWoNoisyJets;
    jpt_2 = j2ptWoNoisyJets;
    jeta_2 = j2etaWoNoisyJets;
    jphi_2 = j2phiWoNoisyJets;
    jcsv_2 = j2csvWoNoisyJets;
    bpt_1 = jb1ptWoNoisyJets;
    beta_1 = jb1etaWoNoisyJets;
    bphi_1 = jb1phiWoNoisyJets;
    bcsv_1 = jb1csvWoNoisyJets;
    bpt_2 = jb2ptWoNoisyJets;
    beta_2 = jb2etaWoNoisyJets;
    bphi_2 = jb2phiWoNoisyJets;
    bcsv_2 = jb2csvWoNoisyJets;
    // NUP (straight from tree)
    id_m_loose_1 = mPFIDLoose;
    id_m_medium_1 = mPFIDMedium;
    id_m_tight_1 = mPFIDTight;
    id_m_loose_2 = tRerunMVArun2v2DBoldDMwLTLoose;
    id_m_medium_2 = tRerunMVArun2v2DBoldDMwLTMedium;
    id_m_tight_2 = tRerunMVArun2v2DBoldDMwLTTight;
    pt_sv = 0;    // calculated later
    eta_sv = 0;   // calculated later
    phi_sv = 0;   // calculated later
    met_sv = 0;   // calculated later
    puweight = 0;

    // need to implement
    // jpfid_1
    // jpuid_1
    // jpfid_2
    // jpuid_2
    // bpfid_1
    // bpuid_1
    // bpfid_2
    // bpuid_2
    npv = nvtx;
    npu = nTruePU;
    // rho (straight from tree)

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
void mutau_tree::set_branches() {
  tree->Branch("evt", &evt);
  tree->Branch("run", &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("dilepton_veto", &dilepton_veto);
  tree->Branch("extraelec_veto", &extraelec_veto);
  tree->Branch("extramuon_veto", &extramuon_veto);
  tree->Branch("trg_singleelectron", &trg_singleelectron);
  tree->Branch("trg_singlemuon", &trg_singlemuon);
  tree->Branch("trg_singletau", &trg_singletau);
  tree->Branch("trg_muonelectron", &trg_muonelectron);
  tree->Branch("trg_mutaucross", &trg_mutaucross);
  tree->Branch("trg_doubletau", &trg_doubletau);
  tree->Branch("flagFilter", &flagFilter);
  tree->Branch("pt_1", &pt_1);
  tree->Branch("phi_1", &phi_1);
  tree->Branch("eta_1", &eta_1);
  tree->Branch("m_1", &m_1);
  tree->Branch("q_1", &q_1);
  tree->Branch("d0_1", &d0_1);
  tree->Branch("dZ_1", &dZ_1);
  tree->Branch("mt_1", &mt_1);
  tree->Branch("pfmt_1", &pfmt_1);
  tree->Branch("puppimt_1", &puppimt_1);
  tree->Branch("iso_1", &iso_1);
  tree->Branch("gen_match_1", &gen_match_1);
  tree->Branch("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1);
  tree->Branch("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1);
  tree->Branch("againstElectronTightMVA6_1", &againstElectronTightMVA6_1);
  tree->Branch("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1);
  tree->Branch("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1);
  tree->Branch("againstMuonLoose3_1", &againstMuonLoose3_1);
  tree->Branch("againstMuonTight3_1", &againstMuonTight3_1);
  tree->Branch("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1);
  tree->Branch("trigweight_1", &trigweight_1);
  tree->Branch("idisoweight_1", &idisoweight_1);
  tree->Branch("pt_2", &pt_2);
  tree->Branch("phi_2", &phi_2);
  tree->Branch("eta_2", &eta_2);
  tree->Branch("m_2", &m_2);
  tree->Branch("q_2", &q_2);
  tree->Branch("d0_2", &d0_2);
  tree->Branch("dZ_2", &dZ_2);
  tree->Branch("mt_2", &mt_2);
  tree->Branch("iso_2", &iso_2);
  tree->Branch("gen_match_2", &gen_match_2);
  tree->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2);
  tree->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2);
  tree->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2);
  tree->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2);
  tree->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2);
  tree->Branch("againstMuonLoose3_2", &againstMuonLoose3_2);
  tree->Branch("againstMuonTight3_2", &againstMuonTight3_2);
  tree->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2);
  tree->Branch("trigweight_2", &trigweight_2);
  tree->Branch("idisoweight_2", &idisoweight_2);
  tree->Branch("pt_tt", &pt_tt);
  tree->Branch("mt_tot", &mt_tot);
  tree->Branch("m_vis", &m_vis);
  tree->Branch("m_sv", &m_sv);
  tree->Branch("mt_sv", &mt_sv);
  tree->Branch("met", &met);
  tree->Branch("metphi", &metphi);
  tree->Branch("puppimet", &puppimet);
  tree->Branch("puppimetphi", &puppimetphi);
  tree->Branch("pzetavis", &pzetavis);
  tree->Branch("pzetamiss", &pzetamiss);
  tree->Branch("pfpzetamiss", &pfpzetamiss);
  tree->Branch("puppipzetamiss", &puppipzetamiss);
  tree->Branch("metcov00", &metcov00);
  tree->Branch("metcov01", &metcov01);
  tree->Branch("metcov10", &metcov10);
  tree->Branch("metcov11", &metcov11);
  tree->Branch("mjj", &mjj);
  tree->Branch("jdeta", &jdeta);
  tree->Branch("njetingap", &njetingap);
  tree->Branch("njetingap20", &njetingap20);
  tree->Branch("jdphi", &jdphi);
  tree->Branch("dijetpt", &dijetpt);
  tree->Branch("dijetphi", &dijetphi);
  tree->Branch("ptvis", &ptvis);
  tree->Branch("nbtag", &nbtag);
  tree->Branch("njets", &njets);
  tree->Branch("njetspt20", &njetspt20);
  tree->Branch("jpt_1", &jpt_1);
  tree->Branch("jeta_1", &jeta_1);
  tree->Branch("jphi_1", &jphi_1);
  tree->Branch("jcsv_1", &jcsv_1);
  tree->Branch("jpt_2", &jpt_2);
  tree->Branch("jeta_2", &jeta_2);
  tree->Branch("jphi_2", &jphi_2);
  tree->Branch("jcsv_2", &jcsv_2);
  tree->Branch("bpt_1", &bpt_1);
  tree->Branch("beta_1", &beta_1);
  tree->Branch("bphi_1", &bphi_1);
  tree->Branch("bcsv_1", &bcsv_1);
  tree->Branch("bpt_2", &bpt_2);
  tree->Branch("beta_2", &beta_2);
  tree->Branch("bphi_2", &bphi_2);
  tree->Branch("bcsv_2", &bcsv_2);
  tree->Branch("puweight", &puweight);
  tree->Branch("NUP", &NUP);
  tree->Branch("weight", &weight);
  tree->Branch("id_m_loose_1", &id_m_loose_1);
  tree->Branch("id_m_medium_1", &id_m_medium_1);
  tree->Branch("id_m_tight_1", &id_m_tight_1);
  tree->Branch("id_m_loose_2", &id_m_loose_2);
  tree->Branch("id_m_medium_2", &id_m_medium_2);
  tree->Branch("id_m_tight_2", &id_m_tight_2);
  tree->Branch("pt_sv", &pt_sv);
  tree->Branch("eta_sv", &eta_sv);
  tree->Branch("phi_sv", &phi_sv);
  tree->Branch("met_sv", &met_sv);
  tree->Branch("jpfid_1", &jpfid_1);
  tree->Branch("jpuid_1", &jpuid_1);
  tree->Branch("jpfid_2", &jpfid_2);
  tree->Branch("jpuid_2", &jpuid_2);
  tree->Branch("jpfid_2", &jpfid_2);
  tree->Branch("jpuid_2", &jpuid_2);
  tree->Branch("bpfid_1", &bpfid_1);
  tree->Branch("bpuid_1", &bpuid_1);
  tree->Branch("bpfid_2", &bpfid_2);
  tree->Branch("bpuid_2", &bpuid_2);
  tree->Branch("npv", &npv);
  tree->Branch("npu", &npu);
  tree->Branch("rho", &rho);

  // input branches
  original->SetBranchAddress("evt", &evt);
  original->SetBranchAddress("run", &Run);
  original->SetBranchAddress("lumi", &Lumi);
  original->SetBranchAddress("dimuonVeto", &dimuonVeto);
  original->SetBranchAddress("eVetoZTTp001dxyzR0", &eVetoZTTp001dxyzR0);
  original->SetBranchAddress("muVetoZTTp001dxyzR0", &muVetoZTTp001dxyzR0);
  original->SetBranchAddress("Flag_BadChargedCandidateFilter", &Flag_BadChargedCandidateFilter);
  original->SetBranchAddress("Flag_BadPFMuonFilter", &Flag_BadPFMuonFilter);
  original->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
  original->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter);
  original->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter);
  original->SetBranchAddress("Flag_badMuons", &Flag_badMuons);
  original->SetBranchAddress("Flag_duplicateMuons", &Flag_duplicateMuons);
  original->SetBranchAddress("Flag_ecalBadCalibFilter", &Flag_ecalBadCalibFilter);
  original->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter);
  original->SetBranchAddress("Flag_globalSuperTightHalo2016Filter", &Flag_globalSuperTightHalo2016Filter);
  original->SetBranchAddress("Flag_globalTightHalo2016Filter", &Flag_globalTightHalo2016Filter);
  original->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices);
  original->SetBranchAddress("mCharge", &mCharge);
  original->SetBranchAddress("mRelPFIsoDBDefaultR04", &mRelPFIsoDBDefaultR04);
  original->SetBranchAddress("mZTTGenMatching", &mZTTGenMatching);
  original->SetBranchAddress("raw_pfMetEt", &raw_pfMetEt);
  original->SetBranchAddress("raw_pfMetPhi", &raw_pfMetPhi);
  original->SetBranchAddress("type1_pfMetEt", &type1_pfMetEt);
  original->SetBranchAddress("type1_pfMetPhi", &type1_pfMetPhi);
  original->SetBranchAddress("puppiMetEt", &puppiMetEt);
  original->SetBranchAddress("puppiMetPhi", &puppiMetPhi);
  original->SetBranchAddress("tAgainstElectronVLooseMVA6", &tAgainstElectronVLooseMVA6);
  original->SetBranchAddress("tAgainstElectronLooseMVA6", &tAgainstElectronLooseMVA6);
  original->SetBranchAddress("tAgainstElectronMediumMVA6", &tAgainstElectronMediumMVA6);
  original->SetBranchAddress("tAgainstElectronTightMVA6", &tAgainstElectronTightMVA6);
  original->SetBranchAddress("tAgainstElectronVTightMVA6", &tAgainstElectronVTightMVA6);
  original->SetBranchAddress("tAgainstMuonLoose3", &tAgainstMuonLoose3);
  original->SetBranchAddress("tAgainstMuonTight3", &tAgainstMuonTight3);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTraw", &tRerunMVArun2v2DBoldDMwLTraw);
  original->SetBranchAddress("vbfMassWoNoisyJets", &vbfMassWoNoisyJets);
  original->SetBranchAddress("jetVeto20WoNoisyJets", &jetVeto20WoNoisyJets);
  original->SetBranchAddress("jetVeto30WoNoisyJets", &jetVeto30WoNoisyJets);
  original->SetBranchAddress("bjetDeepCSVVeto20MediumWoNoisyJets", &bjetDeepCSVVeto20MediumWoNoisyJets);
  original->SetBranchAddress("j1ptWoNoisyJets", &jpt_1);
  original->SetBranchAddress("j1etaWoNoisyJets", &jeta_1);
  original->SetBranchAddress("j1phiWoNoisyJets", &jphi_1);
  original->SetBranchAddress("j1csvWoNoisyJets", &jcsv_1);
  original->SetBranchAddress("j2ptWoNoisyJets", &jpt_2);
  original->SetBranchAddress("j2etaWoNoisyJets", &jeta_2);
  original->SetBranchAddress("j2phiWoNoisyJets", &jphi_2);
  original->SetBranchAddress("j2csvWoNoisyJets", &jcsv_2);
  original->SetBranchAddress("jb1ptWoNoisyJets", &bpt_1);
  original->SetBranchAddress("jb1etaWoNoisyJets", &beta_1);
  original->SetBranchAddress("jb1phiWoNoisyJets", &bphi_1);
  original->SetBranchAddress("jb1csvWoNoisyJets", &bcsv_1);
  original->SetBranchAddress("jb1hadronflavorWoNoisyJets", &jb1hadronflavorWoNoisyJets);
  original->SetBranchAddress("jb2ptWoNoisyJets", &bpt_2);
  original->SetBranchAddress("jb2etaWoNoisyJets", &beta_2);
  original->SetBranchAddress("jb2phiWoNoisyJets", &bphi_2);
  original->SetBranchAddress("jb2csvWoNoisyJets", &bcsv_2);
  original->SetBranchAddress("jb2hadronflavorWoNoisyJets", &jb2hadronflavorWoNoisyJets);
  original->SetBranchAddress("NUP", &NUP);
  original->SetBranchAddress("mPFIDTight", &mPFIDTight);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTLoose", &tRerunMVArun2v2DBoldDMwLTLoose);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTMedium", &tRerunMVArun2v2DBoldDMwLTMedium);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTTight", &tRerunMVArun2v2DBoldDMwLTTight);
  // potentially more
  original->SetBranchAddress("metcov00", &metcov00);
  original->SetBranchAddress("metcov01", &metcov01);
  original->SetBranchAddress("metcov10", &metcov10);
  original->SetBranchAddress("metcov11", &metcov11);
  original->SetBranchAddress("nvtx", &nvtx);
  original->SetBranchAddress("nTruePU", &nTruePU);
  original->SetBranchAddress("rho", &rho);

  original->SetBranchAddress("genpX", &genpX);
  original->SetBranchAddress("genpY", &genpY);
  original->SetBranchAddress("vispX", &vispX);
  original->SetBranchAddress("vispY", &vispY);
}

#endif  // ROOT_SRC_MUTAU_TREE_2017_H_
