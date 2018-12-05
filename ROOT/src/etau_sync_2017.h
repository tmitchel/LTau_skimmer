// Copyright 2018 Tyler Mitchell

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TLorentzVector.h"
#include "RecoilCorrector.h"

class etau_tree {
 private:
  TTree *tree, *original;
  bool isMC, isEmbed;
  std::vector<Int_t> good_events;
  TLorentzVector ele, tau, MET, MET_UESUp, MET_UESDown, MET_JESUp, MET_JESDown;

 public:
  // Member variables

  ULong64_t evt;
  Int_t Run, Lumi, recoil;

  // Selection variables
  Float_t ePt, eEta, ePhi, eMass, tPt, tEta, tPhi, tMass, tZTTGenMatching, tDecayMode, eMatchesEle27Filter, eMatchesEle27Path, Ele27WPTightPass, eMatchesEle32Filter, eMatchesEle32Path, Ele32WPTightPass, eMatchesEle35Filter, eMatchesEle35Path,
      Ele35WPTightPass, eMatchesEle24Tau30Filter, eMatchesEle24Tau30Path, Ele24Tau30Pass, tMatchesEle24Tau30Path, tMatchesEle24Tau30Filter, ePVDZ, ePVDXY, eMVANoisoWP90, ePassesConversionVeto, eMissingHits, tPVDZ,
      tByVLooseIsolationMVArun2v1DBoldDMwLT, tRerunMVArun2v2DBoldDMwLTVLoose, tRerunMVArun2v2DBoldDMwLTVVLoose, tDecayModeFinding, tCharge,
      tAgainstMuonLoose3, tAgainstElectronTightMVA6, muVetoZTTp001dxyzR0, eVetoZTTp001dxyzR0, dielectronVeto, eIsoDB03, tByIsolationMVArun2v1DBoldDMwLTraw, tRerunMVArun2v2DBoldDMwLTraw;

  // Constructed while running
  UInt_t run, lumi;
  Int_t gen_match_1, gen_match_2, njets, nbtag, njetspt20;
  Float_t jetVeto20, jetVeto30;
  Float_t met_px, met_py, extraelec_veto, extramuon_veto, dilepton_veto, pfmetcorr_ex, pfmetcorr_ey, genpX, genpY, vispX, vispY, met, metphi,
          pt_1, eta_1, phi_1, m_1, e_1, px_1, py_1, pz_1, pt_2, eta_2, phi_2, m_2, e_2, px_2, py_2, pz_2;

  // Tau Isolation
  Float_t tByLooseIsolationMVArun2v1DBoldDMwLT, tByMediumIsolationMVArun2v1DBoldDMwLT, tByTightIsolationMVArun2v1DBoldDMwLT, tByVTightIsolationMVArun2v1DBoldDMwLT, tByVVTightIsolationMVArun2v1DBoldDMwLT;
  Float_t tRerunMVArun2v2DBoldDMwLTLoose, tRerunMVArun2v2DBoldDMwLTMedium, tRerunMVArun2v2DBoldDMwLTTight, tRerunMVArun2v2DBoldDMwLTVTight, tRerunMVArun2v2DBoldDMwLTVVTight;

  // Jets
  Float_t vbfMassWoNoisyJets, jpt_1, jeta_1, jphi_1, jcsv_1, jpt_2, jeta_2, jphi_2, jcsv_2, bpt_1, beta_1, bphi_1, bcsv_1, bflavor_1, bpt_2, beta_2, bphi_2, bcsv_2, bflavor_2, bjetDeepCSVVeto20Tight, bjetDeepCSVVeto30Loose, bjetDeepCSVVeto30Medium, bjetDeepCSVVeto30Tight, topQuarkPt1, topQuarkPt2;

  // Gen Ino
  Float_t tZTTGenPt, tZTTGenPhi, tZTTGenEta, tZTTGenDR, tGenDecayMode, tGenEnergy, tGenEta, tGenJetEta, tGenJetPt, tGenMotherEnergy, tGenMotherEta, tGenMotherPdgId, tGenMotherPhi, tGenMotherPt, tGenPdgId, tGenPhi, tGenPt, tGenStatus,
          eGenCharge, eGenDirectPromptTauDecay, eGenEnergy, eGenEta, eGenIsPrompt, eGenMotherPdgId, eGenParticle, eGenPdgId, eGenPhi, eGenPrompt, eGenPromptTauDecay, eGenPt, eGenTauDecay, eGenVZ, eGenVtxPVMatch;

  // Others
  Float_t rho, metcov00, metcov01, metcov10, metcov11, NUP, genM, genpT, genEta, numGenJets, nTruePU, nvtx, eZTTGenMatching, eCharge, eMVANoisoWP80, GenWeight, metSig, metcov00_v2, metcov01_v2, metcov10_v2, metcov11_v2, Rivet_higgsPt, Rivet_nJets30,
  type1_pfMetEt, type1_pfMetPhi, bjetDeepCSVVeto20MediumWoNoisyJets;

  // Flags
  Float_t Flag_BadChargedCandidateFilter, Flag_BadPFMuonFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_HBHENoiseFilter, Flag_HBHENoiseIsoFilter, Flag_badMuons, Flag_duplicateMuons, Flag_ecalBadCalibFilter, Flag_eeBadScFilter,
          Flag_globalSuperTightHalo2016Filter, Flag_globalTightHalo2016Filter, Flag_goodVertices;

  // 2016 Placeholder
  Float_t amcatNLO_weight, eMatchesSingleE25Tight, eMatchesEle25TightFilter, singleE25eta2p1TightPass, eMatchesEle25eta2p1TightPath, tAgainstElectronVLooseMVA6, tAgainstMuonTight3, decayModeFindingNewDMs_2;

  // Member functions
  etau_tree(TTree *orig, TTree *itree, bool isMC, bool isEmbed, Int_t rec);
  virtual ~etau_tree() {}
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
etau_tree::etau_tree(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec) :
tree(itree),
original(Original),
isMC(IsMC),
isEmbed(IsEmbed),
recoil(rec) {
  // read only what is needed for skimming and sorting
  original->SetBranchAddress("evt", &evt);

  // electron variables
  original->SetBranchAddress("ePt", &ePt);
  original->SetBranchAddress("eEta", &eEta);
  original->SetBranchAddress("ePhi", &ePhi);
  original->SetBranchAddress("eMass", &eMass);
  original->SetBranchAddress("ePVDZ", &ePVDZ);
  original->SetBranchAddress("ePVDXY", &ePVDXY);
  original->SetBranchAddress("eIsoDB03", &eIsoDB03);
  original->SetBranchAddress("eMissingHits", &eMissingHits);
  original->SetBranchAddress("ePassesConversionVeto", &ePassesConversionVeto);
  original->SetBranchAddress("eMVANoisoWP90", &eMVANoisoWP90);

  // tau variables
  original->SetBranchAddress("tPt", &tPt);
  original->SetBranchAddress("tEta", &tEta);
  original->SetBranchAddress("tPhi", &tPhi);
  original->SetBranchAddress("tMass", &tMass);
  original->SetBranchAddress("tPVDZ", &tPVDZ);
  original->SetBranchAddress("tCharge", &tCharge);
  original->SetBranchAddress("tDecayMode", &tDecayMode);
  original->SetBranchAddress("tZTTGenMatching", &tZTTGenMatching);
  original->SetBranchAddress("tDecayModeFinding", &tDecayModeFinding);
  original->SetBranchAddress("tAgainstMuonLoose3", &tAgainstMuonLoose3);
  original->SetBranchAddress("tAgainstElectronTightMVA6", &tAgainstElectronTightMVA6);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTraw", &tRerunMVArun2v2DBoldDMwLTraw);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVLoose", &tRerunMVArun2v2DBoldDMwLTVLoose);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVVLoose", &tRerunMVArun2v2DBoldDMwLTVVLoose);
  original->SetBranchAddress("tByIsolationMVArun2v1DBoldDMwLTraw", &tByIsolationMVArun2v1DBoldDMwLTraw);
  original->SetBranchAddress("tByVLooseIsolationMVArun2v1DBoldDMwLT", &tByVLooseIsolationMVArun2v1DBoldDMwLT);

  // 2017 triggers
  original->SetBranchAddress("eMatchesEle27Filter", &eMatchesEle27Filter);
  original->SetBranchAddress("eMatchesEle32Filter", &eMatchesEle32Filter);
  original->SetBranchAddress("eMatchesEle35Filter", &eMatchesEle35Filter);
  original->SetBranchAddress("eMatchesEle24Tau30Filter", &eMatchesEle24Tau30Filter);
  original->SetBranchAddress("tMatchesEle24Tau30Filter", &tMatchesEle24Tau30Filter);
  original->SetBranchAddress("eMatchesEle27Path", &eMatchesEle27Path);
  original->SetBranchAddress("eMatchesEle32Path", &eMatchesEle32Path);
  original->SetBranchAddress("eMatchesEle35Path", &eMatchesEle35Path);
  original->SetBranchAddress("eMatchesEle24Tau30Path", &eMatchesEle24Tau30Path);
  original->SetBranchAddress("tMatchesEle24Tau30Path", &tMatchesEle24Tau30Path);
  original->SetBranchAddress("Ele24Tau30Pass", &Ele24Tau30Pass);
  original->SetBranchAddress("Ele27WPTightPass", &Ele27WPTightPass);
  original->SetBranchAddress("Ele32WPTightPass", &Ele32WPTightPass);
  original->SetBranchAddress("Ele35WPTightPass", &Ele35WPTightPass);

  // other
  original->SetBranchAddress("dielectronVeto", &dielectronVeto);
  original->SetBranchAddress("eVetoZTTp001dxyzR0", &eVetoZTTp001dxyzR0);
  original->SetBranchAddress("muVetoZTTp001dxyzR0", &muVetoZTTp001dxyzR0);
}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void etau_tree::do_skimming(TH1F* cutflow) {
  // declare variables for sorting
  ULong64_t evt_now(0);
  ULong64_t evt_before(1);
  int best_evt(-1);
  std::pair<float, float> eleCandidate, tauCandidate;

  Int_t nevt = (Int_t)original->GetEntries();
  for (auto ievt = 0; ievt < nevt; ievt++) {
    original->GetEntry(ievt);
    evt_now = evt;

    // TLorentzVector ele, tau;
    ele.SetPtEtaPhiM(ePt, eEta, ePhi, eMass);
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

    float el_pt_min(25.), tau_pt_min(23.);

    cutflow->Fill(1., 1.);
    // apply event selection

    auto Ele27 = eMatchesEle27Filter && eMatchesEle27Path && Ele27WPTightPass;
    auto Ele32 = eMatchesEle32Filter && eMatchesEle32Path && Ele32WPTightPass;
    auto Ele35 = eMatchesEle35Filter && eMatchesEle35Path && Ele35WPTightPass;
    auto Cross =  eMatchesEle24Tau30Filter && eMatchesEle24Tau30Path && Ele24Tau30Pass && tMatchesEle24Tau30Path && tMatchesEle24Tau30Filter;

    if (isEmbed || (Ele27 || Ele32 || Ele35 || Cross)) cutflow->Fill(2., 1.);
    else  continue;

    if (!isEmbed || (Ele27WPTightPass || Ele32WPTightPass || Ele35WPTightPass || Ele24Tau30Pass)) cutflow->Fill(3., 1.);
    else  continue;

    if (ePt > el_pt_min && fabs(eEta) < 2.1 && fabs(ePVDZ) < 0.2 && fabs(ePVDXY) < 0.045) cutflow->Fill(4., 1.);  // electron kinematic selection
    else  continue;

    if (eMVANoisoWP90 && ePassesConversionVeto && eMissingHits < 2) cutflow->Fill(5., 1.);  // electron quality selection
    else  continue;

    if (tau.Pt() > tau_pt_min && fabs(tau.Eta()) < 2.3 && fabs(tPVDZ) < 0.2) cutflow->Fill(6., 1.);  // tau kinematic selection
    else  continue;

    if (tRerunMVArun2v2DBoldDMwLTVLoose && tDecayModeFinding > 0 && fabs(tCharge) < 2) cutflow->Fill(7., 1.);  // tau quality selection
    else  continue;

    if (tAgainstMuonLoose3 > 0.5 && tAgainstElectronTightMVA6 > 0.5) cutflow->Fill(8., 1.);  // tau against leptons
    else  continue;

    if (muVetoZTTp001dxyzR0 == 0 && eVetoZTTp001dxyzR0 < 2 && dielectronVeto == 0) cutflow->Fill(9., 1.);  // vetos
    else  continue;

    // implement new sorting per
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
    if (evt_now != evt_before) {  // new event, save the tau candidates
      //   since it is new event, do we have the best entry to save? If yes, save it!
      if ( best_evt > -1  )
        good_events.push_back(best_evt);

      //  this is a new event, so the first tau pair is the best! :)
      best_evt = ievt;
      eleCandidate = std::make_pair(ePt, eIsoDB03);
      tauCandidate  = std::make_pair(tPt,  tRerunMVArun2v2DBoldDMwLTraw);
    } else {  // not a new event
      std::pair<float, float> currEleCandidate(ePt, eIsoDB03);
      std::pair<float, float> currTauCandidate(tPt, tRerunMVArun2v2DBoldDMwLTraw);

      // clause 1, select the pair that has most isolated tau lepton 1
      if (currEleCandidate.second - eleCandidate.second  > 0.0001 ) best_evt = ievt;

      // check if the first tau is the same, and if so - move to clause 2
      if ( fabs(currEleCandidate.second - eleCandidate.second)  <  0.0001 ) {
        // pick up  the pair with the highest pT of the first candidate
        if (currEleCandidate.first - eleCandidate.first > 0.0001 ) best_evt = ievt;
        if ( fabs(currEleCandidate.first -eleCandidate.first) < 0.0001 ) {
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
TTree* etau_tree::fill_tree(RecoilCorrector recoilPFMetCorrector) {
  set_branches();  // get all the branches set up

  // loop through all events pasing skimming/sorting
  for (auto& ievt : good_events) {
    original->GetEntry(ievt);

    // convert from Float_t in FSA to Int_t for analyzer
    run = Run;
    lumi = Lumi;
    gen_match_1 = eZTTGenMatching;
    gen_match_2 = tZTTGenMatching;
    njets = jetVeto30;
    nbtag = bjetDeepCSVVeto20MediumWoNoisyJets;
    njetspt20 = jetVeto20;

    // TLorentzVector ele, tau;
    ele.SetPtEtaPhiM(ePt, eEta, ePhi, eMass);
    tau.SetPtEtaPhiM(tPt, tEta, tPhi, tMass);

    met_px = type1_pfMetEt * cos(type1_pfMetPhi);
    met_py = type1_pfMetEt * sin(type1_pfMetPhi);

    extraelec_veto = eVetoZTTp001dxyzR0 > 1;
    extramuon_veto = muVetoZTTp001dxyzR0 > 0;
    dilepton_veto = dielectronVeto > 0;

    // TLorentzVector ele, tau;
    ele.SetPtEtaPhiM(ePt, eEta, ePhi, eMass);
    tau.SetPtEtaPhiM(tPt, tEta, tPhi, tMass);
    MET.SetPtEtaPhiM(type1_pfMetEt, 0, type1_pfMetPhi, 0);

    pfmetcorr_ex = MET.Px();
    pfmetcorr_ey = MET.Py();

    if (recoil == 1) {
      recoilPFMetCorrector.CorrectByMeanResolution(
          MET.Px(),       // uncorrected type I pf met px (float)
          MET.Py(),       // uncorrected type I pf met py (float)
          genpX,          // generator Z/W/Higgs px (float)
          genpY,          // generator Z/W/Higgs py (float)
          vispX,          // generator visible Z/W/Higgs px (float)
          vispY,          // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)
    } else if (recoil == 2) {
      recoilPFMetCorrector.CorrectByMeanResolution(
          MET.Px(),       // uncorrected type I pf met px (float)
          MET.Py(),       // uncorrected type I pf met py (float)
          genpX,          // generator Z/W/Higgs px (float)
          genpY,          // generator Z/W/Higgs py (float)
          vispX,          // generator visible Z/W/Higgs px (float)
          vispY,          // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)
    }

    MET.SetPxPyPzE(pfmetcorr_ex, pfmetcorr_ey, 0, sqrt(pfmetcorr_ex * pfmetcorr_ex + pfmetcorr_ey * pfmetcorr_ey));

    if (isMC && !isEmbed) {
      // met correction due to tau energy scale
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          MET = MET + tau - 1.007*tau;
          tau *= 1.007;
        } else if (tDecayMode == 1) {
          MET = MET + tau - 0.998*tau;
          tau *= 0.998;
        } else if (tDecayMode == 10) {
          MET = MET + tau - 1.001*tau;
          tau *= 1.001;
        }
      } else if (tZTTGenMatching == 1 || tZTTGenMatching == 3) {
        if (tDecayMode == 0) {
          MET = MET+tau-1.003*tau;
          tau *= 1.003;
        } else if (tDecayMode == 1) {
          MET = MET+tau-1.036*tau;
          tau *= 1.036;
        }
      }
    } else if (isEmbed) {
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          MET = MET + tau - 0.975*tau;
          tau *= 0.975;
        } else if (tDecayMode == 1) {
          MET = MET + tau - 0.975*1.051*tau;
          tau *= 0.975*1.051;
        } else if (tDecayMode == 10) {
          MET = MET + tau - 0.975*0.975*0.975*tau;
          tau *= 0.975*0.975*0.975;
        }
      }
    }

    met = MET.Pt();
    metphi = MET.Phi();
    met_px = MET.Px();
    met_py = MET.Py();

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
void etau_tree::set_branches() {
  // output file branches
  tree->Branch("evt", &evt);
  tree->Branch("run" , &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("gen_match_1", &gen_match_1, "gen_match_1/I");
  tree->Branch("gen_match_2", &gen_match_2, "gen_match_2/I");
  tree->Branch("njets", &njets, "njets/I");
  tree->Branch("nbtag", &nbtag, "nbtag/I");
  tree->Branch("njetspt20", &njetspt20, "njetspt20/I");
  tree->Branch("vbfMassWoNoisyJets", &vbfMassWoNoisyJets, "vbfMassWoNoisyJets/F");

  tree->Branch("eMatchesEle27Filter"     , &eMatchesEle27Filter     , "eMatchesEle27Filter/F");
  tree->Branch("eMatchesEle32Filter"     , &eMatchesEle32Filter     , "eMatchesEle32Filter/F");
  tree->Branch("eMatchesEle35Filter"     , &eMatchesEle35Filter     , "eMatchesEle35Filter/F");
  tree->Branch("eMatchesEle24Tau30Filter", &eMatchesEle24Tau30Filter, "eMatchesEle24Tau30Filter/F");
  tree->Branch("tMatchesEle24Tau30Filter", &tMatchesEle24Tau30Filter, "tMatchesEle24Tau30Filter/F");
  tree->Branch("eMatchesEle27Path"       , &eMatchesEle27Path       , "eMatchesEle27Path/F");
  tree->Branch("eMatchesEle32Path"       , &eMatchesEle32Path       , "eMatchesEle32Path/F");
  tree->Branch("eMatchesEle35Path"       , &eMatchesEle35Path       , "eMatchesEle35Path/F");
  tree->Branch("eMatchesEle24Tau30Path"  , &eMatchesEle24Tau30Path  , "eMatchesEle24Tau30Path/F");
  tree->Branch("tMatchesEle24Tau30Path"  , &tMatchesEle24Tau30Path  , "tMatchesEle24Tau30Path/F");
  tree->Branch("Ele24Tau30Pass"          , &Ele24Tau30Pass          , "Ele24Tau30Pass/F");
  tree->Branch("Ele27WPTightPass"        , &Ele27WPTightPass        , "Ele27WPTightPass/F");
  tree->Branch("Ele32WPTightPass"        , &Ele32WPTightPass        , "Ele32WPTightPass/F");
  tree->Branch("Ele35WPTightPass"        , &Ele35WPTightPass        , "Ele35WPTightPass/F");

  tree->Branch("met_px", &met_px, "met_px/F");
  tree->Branch("met_py", &met_py, "met_py/F");
  tree->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/F");
  tree->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/F");
  tree->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/F");
  tree->Branch("pfmetcorr_ex", &pfmetcorr_ex, "pfmetcorr_ex/F");
  tree->Branch("pfmetcorr_ey", &pfmetcorr_ey, "pfmetcorr_ey/F");
  tree->Branch("genpX", &genpX, "genpX/F");
  tree->Branch("genpY", &genpY, "genpY/F");
  tree->Branch("vispX", &vispX, "vispX/F");
  tree->Branch("vispY", &vispY, "vispY/F");
  tree->Branch("met", &met, "met/F");
  tree->Branch("metphi", &metphi, "metphi/F");
  tree->Branch("met_px", &met_px, "met_px/F");
  tree->Branch("met_py", &met_py, "met_py/F");
  tree->Branch("pt_1" , &pt_1 , "pt_1/F");
  tree->Branch("eta_1", &eta_1, "eta_1/F");
  tree->Branch("phi_1", &phi_1, "phi_1/F");
  tree->Branch("m_1"  , &m_1  , "m_1/F");
  tree->Branch("e_1"  , &e_1  , "e_1/F");
  tree->Branch("px_1" , &px_1 , "px_1/F");
  tree->Branch("py_1" , &py_1 , "py_1/F");
  tree->Branch("pz_1" , &pz_1 , "pz_1/F");
  tree->Branch("dZ_1" , &ePVDZ, "dZ_1/F");
  tree->Branch("d0_1" , &ePVDXY, "d0_1/F");
  tree->Branch("q_1"  , &eCharge, "q_1/F");
  tree->Branch("iso_1", &eIsoDB03, "iso_1/F");
  tree->Branch("NoisoID80_1", &eMVANoisoWP80, "NoisoID80_1/F");
  tree->Branch("pt_2" , &pt_2 , "pt_2/F");
  tree->Branch("eta_2", &eta_2, "eta_2/F");
  tree->Branch("phi_2", &phi_2, "phi_2/F");
  tree->Branch("m_2"  , &m_2  , "m_2/F");
  tree->Branch("e_2"  , &e_2  , "e_2/F");
  tree->Branch("px_2" , &px_2 , "px_2/F");
  tree->Branch("py_2" , &py_2 , "py_2/F");
  tree->Branch("pz_2" , &pz_2 , "pz_2/F");
  tree->Branch("dZ_2" , &tPVDZ , "dZ_2/F");
  tree->Branch("q_2"  , &tCharge, "q_2/F");
  tree->Branch("iso_2", &tRerunMVArun2v2DBoldDMwLTraw, "iso_2/F");
  tree->Branch("decayModeFinding_2", &tDecayModeFinding, "decayModeFinding_2/F");
  tree->Branch("l2_decayMode", &tDecayMode, "l2_decayMode/F");

  tree->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_2"  , &tByLooseIsolationMVArun2v1DBoldDMwLT     , "byLooseIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_2" , &tByMediumIsolationMVArun2v1DBoldDMwLT    , "byMediumIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2"  , &tByTightIsolationMVArun2v1DBoldDMwLT     , "byTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byVTightIsolationMVArun2v1DBoldDMwLT_2" , &tByVTightIsolationMVArun2v1DBoldDMwLT    , "byVTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byVVTightIsolationMVArun2v1DBoldDMwLT_2", &tByVVTightIsolationMVArun2v1DBoldDMwLT   , "byVVTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVLoose"        , &tRerunMVArun2v2DBoldDMwLTVLoose          , "tRerunMVArun2v2DBoldDMwLTVLoose/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTLoose"         , &tRerunMVArun2v2DBoldDMwLTLoose           , "tRerunMVArun2v2DBoldDMwLTLoose/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTMedium"        , &tRerunMVArun2v2DBoldDMwLTMedium          , "tRerunMVArun2v2DBoldDMwLTMedium/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTTight"         , &tRerunMVArun2v2DBoldDMwLTTight           , "tRerunMVArun2v2DBoldDMwLTTight/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVTight"        , &tRerunMVArun2v2DBoldDMwLTVTight          , "tRerunMVArun2v2DBoldDMwLTVTight/F");
  tree->Branch("tRerunMVArun2v2DBoldDMwLTVVTight"       , &tRerunMVArun2v2DBoldDMwLTVVTight         , "tRerunMVArun2v2DBoldDMwLTVVTight/F");

  tree->Branch("rho"        , &rho        , "rho/F");
  tree->Branch("metcov00"   , &metcov00   , "metcov00/F");
  tree->Branch("metcov01"   , &metcov01   , "metcov01/F");
  tree->Branch("metcov10"   , &metcov10   , "metcov10/F");
  tree->Branch("metcov11"   , &metcov11   , "metcov11/F");
  tree->Branch("metcov00_v2", &metcov00_v2, "metcov00_v2/F");
  tree->Branch("metcov01_v2", &metcov01_v2, "metcov01_v2/F");
  tree->Branch("metcov10_v2", &metcov10_v2, "metcov10_v2/F");
  tree->Branch("metcov11_v2", &metcov11_v2, "metcov11_v2/F");
  tree->Branch("NUP", &NUP, "NUP/F");
  tree->Branch("genM", &genM, "genM/F");
  tree->Branch("genpT", &genpT, "genpT/F");
  tree->Branch("genEta", &genEta, "genEta/F");
  tree->Branch("numGenJets", &numGenJets, "numGenJets/F");
  tree->Branch("npu", &nTruePU, "npu/F");
  tree->Branch("npv", &nvtx, "npv/F");
  tree->Branch("eMVANoisoWP80", &eMVANoisoWP80, "eMVANoisoWP80/F");
  tree->Branch("genweight", &GenWeight, "genweight/F");
  tree->Branch("metSig", &metSig, "metSig/F");
  tree->Branch("Rivet_higgsPt", &Rivet_higgsPt, "Rivet_higgsPt/F");
  tree->Branch("Rivet_nJets30", &Rivet_nJets30, "Rivet_nJets30/F");

  tree->Branch("jpt_1" , &jpt_1 , "jpt_1/F");
  tree->Branch("jeta_1", &jeta_1, "jeta_1/F");
  tree->Branch("jphi_1", &jphi_1, "jphi_1/F");
  tree->Branch("jcsv_1", &jcsv_1, "jcsv_1/F");
  tree->Branch("jpt_2" , &jpt_2 , "jpt_2/F");
  tree->Branch("jeta_2", &jeta_2, "jeta_2/F");
  tree->Branch("jphi_2", &jphi_2, "jphi_2/F");
  tree->Branch("jcsv_2", &jcsv_2, "jcsv_2/F");
  tree->Branch("bpt_1" , &bpt_1 , "bpt_1/F");
  tree->Branch("beta_1", &beta_1, "beta_1/F");
  tree->Branch("bphi_1", &bphi_1, "bphi_1/F");
  tree->Branch("bcsv_1", &bcsv_1, "bcsv_1/F");
  tree->Branch("bflavor_1", &bflavor_1, "bflavor_1/F");
  tree->Branch("bpt_2" , &bpt_2 , "bpt_2/F");
  tree->Branch("beta_2", &beta_2, "beta_2/F");
  tree->Branch("bphi_2", &bphi_2, "bphi_2/F");
  tree->Branch("bcsv_2", &bcsv_2, "bcsv_2/F");
  tree->Branch("bflavor_2", &bflavor_2, "bflavor_2/F");

  tree->Branch("bjetDeepCSVVeto20Tight", &bjetDeepCSVVeto20Tight, "bjetDeepCSVVeto20Tight/F");
  tree->Branch("bjetDeepCSVVeto30Loose", &bjetDeepCSVVeto30Loose, "bjetDeepCSVVeto30Loose/F");
  tree->Branch("bjetDeepCSVVeto30Medium", &bjetDeepCSVVeto30Medium, "bjetDeepCSVVeto30Medium/F");
  tree->Branch("bjetDeepCSVVeto30Tight", &bjetDeepCSVVeto30Tight, "bjetDeepCSVVeto30Tight/F");

  tree->Branch("topQuarkPt1", &topQuarkPt1, "topQuarkPt1/F");
  tree->Branch("topQuarkPt2", &topQuarkPt2, "topQuarkPt2/F");

  tree->Branch("tZTTGenPt", &tZTTGenPt, "tZTTGenPt/F");
  tree->Branch("tZTTGenPhi", &tZTTGenPhi, "tZTTGenPhi/F");
  tree->Branch("tZTTGenEta", &tZTTGenEta, "tZTTGenEta/F");
  tree->Branch("tZTTGenDR", &tZTTGenDR, "tZTTGenDR/F");
  tree->Branch("tGenDecayMode", &tGenDecayMode, "tGenDecayMode/F");
  tree->Branch("tGenEnergy", &tGenEnergy, "tGenEnergy/F");
  tree->Branch("tGenEta", &tGenEta, "tGenEta/F");
  tree->Branch("tGenJetEta", &tGenJetEta, "tGenJetEta/F");
  tree->Branch("tGenJetPt", &tGenJetPt, "tGenJetPt/F");
  tree->Branch("tGenMotherEnergy", &tGenMotherEnergy, "tGenMotherEnergy/F");
  tree->Branch("tGenMotherEta", &tGenMotherEta, "tGenMotherEta/F");
  tree->Branch("tGenMotherPdgId", &tGenMotherPdgId, "tGenMotherPdgId/F");
  tree->Branch("tGenMotherPhi", &tGenMotherPhi, "tGenMotherPhi/F");
  tree->Branch("tGenMotherPt", &tGenMotherPt, "tGenMotherPt/F");
  tree->Branch("tGenPdgId", &tGenPdgId, "tGenPdgId/F");
  tree->Branch("tGenPhi", &tGenPhi, "tGenPhi/F");
  tree->Branch("tGenPt", &tGenPt, "tGenPt/F");
  tree->Branch("tGenStatus", &tGenStatus, "tGenStatus/F");
  tree->Branch("eGenCharge", &eGenCharge, "eGenCharge/F");
  tree->Branch("eGenDirectPromptTauDecay", &eGenDirectPromptTauDecay, "eGenDirectPromptTauDecay/F");
  tree->Branch("eGenEnergy", &eGenEnergy, "eGenEnergy/F");
  tree->Branch("eGenEta", &eGenEta, "eGenEta/F");
  tree->Branch("eGenIsPrompt", &eGenIsPrompt, "eGenIsPrompt/F");
  tree->Branch("eGenMotherPdgId", &eGenMotherPdgId, "eGenMotherPdgId/F");
  tree->Branch("eGenParticle", &eGenParticle, "eGenParticle/F");
  tree->Branch("eGenPdgId", &eGenPdgId, "eGenPdgId/F");
  tree->Branch("eGenPhi", &eGenPhi, "eGenPhi/F");
  tree->Branch("eGenPrompt", &eGenPrompt, "eGenPrompt/F");
  tree->Branch("eGenPromptTauDecay", &eGenPromptTauDecay, "eGenPromptTauDecay/F");
  tree->Branch("eGenPt", &eGenPt, "eGenPt/F");
  tree->Branch("eGenTauDecay", &eGenTauDecay, "eGenTauDecay/F");
  tree->Branch("eGenVZ", &eGenVZ, "eGenVZ/F");
  tree->Branch("eGenVtxPVMatch", &eGenVtxPVMatch, "eGenVtxPVMatch/F");

  // Flags
  tree->Branch("Flag_BadChargedCandidateFilter"         , &Flag_BadChargedCandidateFilter         , "Flag_BadChargedCandidateFilter/F");
  tree->Branch("Flag_BadPFMuonFilter"                   , &Flag_BadPFMuonFilter                   , "Flag_BadPFMuonFilter/F");
  tree->Branch("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, "Flag_EcalDeadCellTriggerPrimitiveFilter/F");
  tree->Branch("Flag_HBHENoiseFilter"                   , &Flag_HBHENoiseFilter                   , "Flag_HBHENoiseFilter/F");
  tree->Branch("Flag_HBHENoiseIsoFilter"                , &Flag_HBHENoiseIsoFilter                , "Flag_HBHENoiseIsoFilter/F");
  tree->Branch("Flag_badMuons"                          , &Flag_badMuons                          , "Flag_badMuons/F");
  tree->Branch("Flag_duplicateMuons"                    , &Flag_duplicateMuons                    , "Flag_duplicateMuons/F");
  tree->Branch("Flag_ecalBadCalibFilter"                , &Flag_ecalBadCalibFilter                , "Flag_ecalBadCalibFilter/F");
  tree->Branch("Flag_eeBadScFilter"                     , &Flag_eeBadScFilter                     , "Flag_eeBadScFilter/F");
  tree->Branch("Flag_globalSuperTightHalo2016Filter"    , &Flag_globalSuperTightHalo2016Filter    , "Flag_globalSuperTightHalo2016Filter/F");
  tree->Branch("Flag_globalTightHalo2016Filter"         , &Flag_globalTightHalo2016Filter         , "Flag_globalTightHalo2016Filter/F");
  tree->Branch("Flag_goodVertices"                      , &Flag_goodVertices                      , "Flag_goodVertices/F");

  // 2016 placeholders
  tree->Branch("amcatNLO_weight", &amcatNLO_weight, "amcatNLO_weight/F");
  tree->Branch("eMatchesSingleE25Tight", &eMatchesSingleE25Tight, "eMatchesSingleE25Tight/F");
  tree->Branch("eMatchesEle25TightFilter", &eMatchesEle25TightFilter, "eMatchesEle25TightFilter/F");
  tree->Branch("singleE25eta2p1TightPass", &singleE25eta2p1TightPass, "singleE25eta2p1TightPass/F");
  tree->Branch("againstElectronTightMVA6_2", &tAgainstElectronTightMVA6, "againstElectronTightMVA6_2/F");
  tree->Branch("againstElectronVLooseMVA6_2", &tAgainstElectronVLooseMVA6, "againstElectronVLooseMVA6_2/F");
  tree->Branch("againstMuonTight3_2", &tAgainstMuonTight3, "againstMuonTight3_2/F");
  tree->Branch("againstMuonLoose3_2", &tAgainstMuonLoose3, "againstMuonLoose3_2/F");
  tree->Branch("byVLooseIsolationMVArun2v1DBoldDMwLT_2", &tByVLooseIsolationMVArun2v1DBoldDMwLT, "byVLooseIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("decayModeFindingNewDMs_2", &decayModeFindingNewDMs_2, "decayModeFindingNewDMs_2/F");

  // input branches
  // event info
  original->SetBranchAddress("evt", &evt);
  original->SetBranchAddress("run", &Run);
  original->SetBranchAddress("lumi", &Lumi);
  original->SetBranchAddress("rho", &rho);
  original->SetBranchAddress("genpX", &genpX);
  original->SetBranchAddress("genpY", &genpY);
  original->SetBranchAddress("vispX", &vispX);
  original->SetBranchAddress("vispY", &vispY);
  original->SetBranchAddress("nTruePU", &nTruePU);
  original->SetBranchAddress("nvtx", &nvtx);
  original->SetBranchAddress("NUP", &NUP);
  original->SetBranchAddress("GenWeight", &GenWeight);
  original->SetBranchAddress("numGenJets", &numGenJets);
  original->SetBranchAddress("genpT", &genpT);
  original->SetBranchAddress("genEta", &genEta);
  original->SetBranchAddress("genM", &genM);
  original->SetBranchAddress("Rivet_higgsPt", &Rivet_higgsPt);
  original->SetBranchAddress("Rivet_nJets30", &Rivet_nJets30);
  original->SetBranchAddress("vbfMassWoNoisyJets", &vbfMassWoNoisyJets);

  // electron branches
  original->SetBranchAddress("eMVANoisoWP80", &eMVANoisoWP80);
  original->SetBranchAddress("eZTTGenMatching", &eZTTGenMatching);
  original->SetBranchAddress("eCharge", &eCharge);

  // tau branches
  original->SetBranchAddress("tZTTGenMatching", &tZTTGenMatching);
  original->SetBranchAddress("tByLooseIsolationMVArun2v1DBoldDMwLT", &tByLooseIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tByMediumIsolationMVArun2v1DBoldDMwLT", &tByMediumIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tByTightIsolationMVArun2v1DBoldDMwLT", &tByTightIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tByVTightIsolationMVArun2v1DBoldDMwLT", &tByVTightIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tByVVTightIsolationMVArun2v1DBoldDMwLT", &tByVVTightIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTLoose"         , &tRerunMVArun2v2DBoldDMwLTLoose);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTMedium"        , &tRerunMVArun2v2DBoldDMwLTMedium);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTTight"         , &tRerunMVArun2v2DBoldDMwLTTight);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVTight"        , &tRerunMVArun2v2DBoldDMwLTVTight);
  original->SetBranchAddress("tRerunMVArun2v2DBoldDMwLTVVTight"       , &tRerunMVArun2v2DBoldDMwLTVVTight);

  // jet branches
  original->SetBranchAddress("jetVeto20WoNoisyJets", &jetVeto20);
  original->SetBranchAddress("jetVeto30WoNoisyJets", &jetVeto30);

  original->SetBranchAddress("j1ptWoNoisyJets", &jpt_1);
  original->SetBranchAddress("j1etaWoNoisyJets", &jeta_1);
  original->SetBranchAddress("j1phiWoNoisyJets", &jphi_1);
  original->SetBranchAddress("j1csvWoNoisyJets", &jcsv_1);
  original->SetBranchAddress("j2ptWoNoisyJets", &jpt_2);
  original->SetBranchAddress("j2etaWoNoisyJets", &jeta_2);
  original->SetBranchAddress("j2phiWoNoisyJets", &jphi_2);
  original->SetBranchAddress("j2csvWoNoisyJets", &jcsv_2);
  original->SetBranchAddress("jb1ptWoNoisyJets",  &bpt_1);
  original->SetBranchAddress("jb1etaWoNoisyJets", &beta_1);
  original->SetBranchAddress("jb1phiWoNoisyJets", &bphi_1);
  original->SetBranchAddress("jb1csvWoNoisyJets", &bcsv_1);
  original->SetBranchAddress("jb1hadronflavorWoNoisyJets", &bflavor_1);
  original->SetBranchAddress("jb2ptWoNoisyJets",  &bpt_2);
  original->SetBranchAddress("jb2etaWoNoisyJets", &beta_2);
  original->SetBranchAddress("jb2phiWoNoisyJets", &bphi_2);
  original->SetBranchAddress("jb2csvWoNoisyJets", &bcsv_2);
  original->SetBranchAddress("jb2hadronflavorWoNoisyJets", &bflavor_2);

  original->SetBranchAddress("bjetDeepCSVVeto20MediumWoNoisyJets", &bjetDeepCSVVeto20MediumWoNoisyJets);
  original->SetBranchAddress("bjetDeepCSVVeto20Tight", &bjetDeepCSVVeto20Tight);
  original->SetBranchAddress("bjetDeepCSVVeto30Loose", &bjetDeepCSVVeto30Loose);
  original->SetBranchAddress("bjetDeepCSVVeto30Medium", &bjetDeepCSVVeto30Medium);
  original->SetBranchAddress("bjetDeepCSVVeto30Tight", &bjetDeepCSVVeto30Tight);

  original->SetBranchAddress("topQuarkPt1", &topQuarkPt1);
  original->SetBranchAddress("topQuarkPt2", &topQuarkPt2);

  original->SetBranchAddress("tZTTGenPt", &tZTTGenPt);
  original->SetBranchAddress("tZTTGenPhi", &tZTTGenPhi);
  original->SetBranchAddress("tZTTGenEta", &tZTTGenEta);
  original->SetBranchAddress("tZTTGenDR", &tZTTGenDR);
  original->SetBranchAddress("tGenDecayMode", &tGenDecayMode);
  original->SetBranchAddress("tGenEnergy", &tGenEnergy);
  original->SetBranchAddress("tGenEta", &tGenEta);
  original->SetBranchAddress("tGenJetEta", &tGenJetEta);
  original->SetBranchAddress("tGenJetPt", &tGenJetPt);
  original->SetBranchAddress("tGenMotherEnergy", &tGenMotherEnergy);
  original->SetBranchAddress("tGenMotherEta", &tGenMotherEta);
  original->SetBranchAddress("tGenMotherPdgId", &tGenMotherPdgId);
  original->SetBranchAddress("tGenMotherPhi", &tGenMotherPhi);
  original->SetBranchAddress("tGenMotherPt", &tGenMotherPt);
  original->SetBranchAddress("tGenPdgId", &tGenPdgId);
  original->SetBranchAddress("tGenPhi", &tGenPhi);
  original->SetBranchAddress("tGenPt", &tGenPt);
  original->SetBranchAddress("tGenStatus", &tGenStatus);
  original->SetBranchAddress("eGenCharge", &eGenCharge);
  original->SetBranchAddress("eGenDirectPromptTauDecay", &eGenDirectPromptTauDecay);
  original->SetBranchAddress("eGenEnergy", &eGenEnergy);
  original->SetBranchAddress("eGenEta", &eGenEta);
  original->SetBranchAddress("eGenIsPrompt", &eGenIsPrompt);
  original->SetBranchAddress("eGenMotherPdgId", &eGenMotherPdgId);
  original->SetBranchAddress("eGenParticle", &eGenParticle);
  original->SetBranchAddress("eGenPdgId", &eGenPdgId);
  original->SetBranchAddress("eGenPhi", &eGenPhi);
  original->SetBranchAddress("eGenPrompt", &eGenPrompt);
  original->SetBranchAddress("eGenPromptTauDecay", &eGenPromptTauDecay);
  original->SetBranchAddress("eGenPt", &eGenPt);
  original->SetBranchAddress("eGenTauDecay", &eGenTauDecay);
  original->SetBranchAddress("eGenVZ", &eGenVZ);
  original->SetBranchAddress("eGenVtxPVMatch", &eGenVtxPVMatch);

  // MET branches
  original->SetBranchAddress("type1_pfMetEt", &type1_pfMetEt);
  original->SetBranchAddress("type1_pfMetPhi", &type1_pfMetPhi);
  original->SetBranchAddress("metSig", &metSig);
  original->SetBranchAddress("metcov00"   , &metcov00);
  original->SetBranchAddress("metcov01"   , &metcov01);
  original->SetBranchAddress("metcov10"   , &metcov10);
  original->SetBranchAddress("metcov11"   , &metcov11);
  original->SetBranchAddress("metcov00_DESYlike", &metcov00_v2);
  original->SetBranchAddress("metcov01_DESYlike", &metcov01_v2);
  original->SetBranchAddress("metcov10_DESYlike", &metcov10_v2);
  original->SetBranchAddress("metcov11_DESYlike", &metcov11_v2);

  // Flags
  original->SetBranchAddress("Flag_BadChargedCandidateFilter"         , &Flag_BadChargedCandidateFilter);
  original->SetBranchAddress("Flag_BadPFMuonFilter"                   , &Flag_BadPFMuonFilter);
  original->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter);
  original->SetBranchAddress("Flag_HBHENoiseFilter"                   , &Flag_HBHENoiseFilter);
  original->SetBranchAddress("Flag_HBHENoiseIsoFilter"                , &Flag_HBHENoiseIsoFilter);
  original->SetBranchAddress("Flag_badMuons"                          , &Flag_badMuons);
  original->SetBranchAddress("Flag_duplicateMuons"                    , &Flag_duplicateMuons);
  original->SetBranchAddress("Flag_ecalBadCalibFilter"                , &Flag_ecalBadCalibFilter);
  original->SetBranchAddress("Flag_eeBadScFilter"                     , &Flag_eeBadScFilter);
  original->SetBranchAddress("Flag_globalSuperTightHalo2016Filter"    , &Flag_globalSuperTightHalo2016Filter);
  original->SetBranchAddress("Flag_globalTightHalo2016Filter"         , &Flag_globalTightHalo2016Filter);
  original->SetBranchAddress("Flag_goodVertices"                      , &Flag_goodVertices);
}
