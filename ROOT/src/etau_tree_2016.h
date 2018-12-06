// Copyright 2018 Tyler Mitchell

#include <iostream>
#include <cmath>
#include <vector>
#include <utility>
#include "TTree.h"
#include "TLorentzVector.h"
#include "RecoilCorrector.h"
#include "base_tree.h"

class etau_tree2016: public virtual base_tree {
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
          Ele35WPTightPass, eMatchesEle24Tau30Filter, eMatchesEle24Tau30Path, Ele24Tau30Pass, tMatchesEle24Tau30Path, tMatchesEle24Tau30Filter, ePVDZ, ePVDXY, eMVANonTrigWP80, eMVANoisoWP90, ePassesConversionVeto, eMissingHits, tPVDZ,
          tByVLooseIsolationMVArun2v1DBoldDMwLT, tDecayModeFinding, tCharge, e_t_DR,
          tAgainstMuonLoose3, tAgainstElectronTightMVA6, muVetoZTTp001dxyzR0, eVetoZTTp001dxyzR0, dielectronVeto, eIsoDB03, tByIsolationMVArun2v1DBoldDMwLTraw, tRerunMVArun2v1DBoldDMwLTVLoose;

  // Constructed while running
  UInt_t run, lumi;
  Int_t gen_match_1, gen_match_2, njets, nbtag, njetspt20;
  Float_t jetVeto20, jetVeto30;
  Float_t met_px, met_py, extraelec_veto, extramuon_veto, dilepton_veto, pfmetcorr_ex, pfmetcorr_ey, genpX, genpY, vispX, vispY, met, metphi,
          pt_1, eta_1, phi_1, m_1, e_1, px_1, py_1, pz_1, pt_2, eta_2, phi_2, m_2, e_2, px_2, py_2, pz_2;
  Float_t pfmetcorr_ex_UESUp, pfmetcorr_ey_UESUp, pfmetcorr_ex_UESDown, pfmetcorr_ey_UESDown, pfmetcorr_ex_JESUp, pfmetcorr_ey_JESUp, pfmetcorr_ex_JESDown, pfmetcorr_ey_JESDown;

  // Tau Isolation
  Float_t tByLooseIsolationMVArun2v1DBoldDMwLT, tByMediumIsolationMVArun2v1DBoldDMwLT, tByTightIsolationMVArun2v1DBoldDMwLT, tByVTightIsolationMVArun2v1DBoldDMwLT, tByVVTightIsolationMVArun2v1DBoldDMwLT;

  // Jets
  Float_t vbfMass, jpt_1, jeta_1, jphi_1, jcsv_1, jpt_2, jeta_2, jphi_2, jcsv_2, bpt_1, beta_1, bphi_1, bcsv_1, bflavor_1, bpt_2, beta_2, bphi_2, bcsv_2, bflavor_2, bjetDeepCSVVeto20Tight, bjetDeepCSVVeto30Loose, bjetDeepCSVVeto30Medium, bjetDeepCSVVeto30Tight, topQuarkPt1, topQuarkPt2, bjetCISVVeto20Tight, bjetCISVVeto30Loose, bjetCISVVeto30Medium, bjetCISVVeto30Tight, bjetCISVVeto20Medium;

  // Gen Ino
  Float_t tZTTGenPt, tZTTGenPhi, tZTTGenEta, tZTTGenDR, tGenDecayMode, tGenEnergy, tGenEta, tGenJetEta, tGenJetPt, tGenMotherEnergy, tGenMotherEta, tGenMotherPdgId, tGenMotherPhi, tGenMotherPt, tGenPdgId, tGenPhi, tGenPt, tGenStatus,
          eGenCharge, eGenDirectPromptTauDecay, eGenEnergy, eGenEta, eGenIsPrompt, eGenMotherPdgId, eGenParticle, eGenPdgId, eGenPhi, eGenPrompt, eGenPromptTauDecay, eGenPt, eGenTauDecay, eGenVZ, eGenVtxPVMatch;

  // Others
  Float_t rho, metcov00, metcov01, metcov10, metcov11, NUP, genM, genpT, genEta, numGenJets, nTruePU, nvtx, eZTTGenMatching, eCharge, eMVANoisoWP80, GenWeight, metSig, metcov00_v2, metcov01_v2, metcov10_v2, metcov11_v2, Rivet_higgsPt, Rivet_nJets30,
  type1_pfMetEt, type1_pfMetPhi, bjetDeepCSVVeto20Medium;

  // Others Uncertainty
  Float_t type1_pfMet_shiftedPt_UnclusteredEnUp, type1_pfMet_shiftedPhi_UnclusteredEnUp, type1_pfMet_shiftedPt_UnclusteredEnDown, type1_pfMet_shiftedPhi_UnclusteredEnDown,
          type1_pfMet_shiftedPt_JetEnUp, type1_pfMet_shiftedPhi_JetEnUp, type1_pfMet_shiftedPt_JetEnDown, type1_pfMet_shiftedPhi_JetEnDown,
          type1_pfMet_shiftedPt_JetResUp, type1_pfMet_shiftedPhi_JetResUp, type1_pfMet_shiftedPt_JetResDown, type1_pfMet_shiftedPhi_JetResDown;

  // Flags
  Float_t Flag_BadChargedCandidateFilter, Flag_BadPFMuonFilter, Flag_EcalDeadCellTriggerPrimitiveFilter, Flag_HBHENoiseFilter, Flag_HBHENoiseIsoFilter, Flag_badMuons, Flag_duplicateMuons, Flag_ecalBadCalibFilter, Flag_eeBadScFilter,
          Flag_globalSuperTightHalo2016Filter, Flag_globalTightHalo2016Filter, Flag_goodVertices;

  // Systematics
  Float_t jetVeto30_JetEta0to3Down, jetVeto30_JetEta0to3Up, jetVeto30_JetEta0to5Down, jetVeto30_JetEta0to5Up, jetVeto30_JetEta3to5Down, jetVeto30_JetEta3to5Up,
          jetVeto30_JetAbsoluteStatDown, jetVeto30_JetAbsoluteStatUp, jetVeto30_JetClosureDown, jetVeto30_JetClosureUp, jetVeto30_JetEnDown, jetVeto30_JetFlavorQCDDown, jetVeto30_JetFlavorQCDUp, jetVeto30_JetFragmentationDown,
          jetVeto30_JetFragmentationUp, jetVeto30_JetPileUpDataMCDown, jetVeto30_JetPileUpDataMCUp, jetVeto30_JetPileUpPtBBDown, jetVeto30_JetPileUpPtBBUp, jetVeto30_JetPileUpPtEC1Down, jetVeto30_JetPileUpPtEC1Up, jetVeto30_JetPileUpPtEC2Down,
          jetVeto30_JetPileUpPtEC2Up, jetVeto30_JetPileUpPtHFDown, jetVeto30_JetPileUpPtHFUp, jetVeto30_JetPileUpPtRefDown, jetVeto30_JetPileUpPtRefUp, jetVeto30_JetRelativeBalDown, jetVeto30_JetRelativeBalUp, jetVeto30_JetRelativeFSRDown,
          jetVeto30_JetRelativeFSRUp, jetVeto30_JetRelativeJEREC1Down, jetVeto30_JetRelativeJEREC1Up, jetVeto30_JetRelativeJEREC2Down, jetVeto30_JetRelativeJEREC2Up, jetVeto30_JetRelativeJERHFDown, jetVeto30_JetRelativeJERHFUp,
          jetVeto30_JetRelativePtBBDown, jetVeto30_JetRelativePtBBUp, jetVeto30_JetRelativePtEC1Down, jetVeto30_JetRelativePtEC1Up, jetVeto30_JetRelativePtEC2Down, jetVeto30_JetRelativePtEC2Up, jetVeto30_JetRelativePtHFDown,
          jetVeto30_JetRelativePtHFUp, jetVeto30_JetRelativeSampleDown, jetVeto30_JetRelativeSampleUp, jetVeto30_JetRelativeStatECDown, jetVeto30_JetRelativeStatECUp, jetVeto30_JetRelativeStatFSRDown, jetVeto30_JetRelativeStatFSRUp,
          jetVeto30_JetRelativeStatHFDown, jetVeto30_JetRelativeStatHFUp, jetVeto30_JetSinglePionECALDown, jetVeto30_JetSinglePionECALUp, jetVeto30_JetSinglePionHCALDown, jetVeto30_JetSinglePionHCALUp, jetVeto30_JetTimePtEtaDown,
          jetVeto30_JetTimePtEtaUp, jetVeto30_JetTotalDown, jetVeto30_JetTotalUp;

  Float_t vbfMass_JetEta0to3Down, vbfMass_JetEta0to3Up, vbfMass_JetEta0to5Down, vbfMass_JetEta0to5Up, vbfMass_JetEta3to5Down, vbfMass_JetEta3to5Up,
          vbfMass_JetAbsoluteFlavMapDown, vbfMass_JetAbsoluteFlavMapUp, vbfMass_JetAbsoluteMPFBiasDown, vbfMass_JetAbsoluteMPFBiasUp, vbfMass_JetAbsoluteScaleDown, vbfMass_JetAbsoluteScaleUp,
          vbfMass_JetAbsoluteStatDown, vbfMass_JetAbsoluteStatUp, vbfMass_JetClosureDown, vbfMass_JetClosureUp, vbfMass_JetFlavorQCDDown, vbfMass_JetFlavorQCDUp, vbfMass_JetFragmentationDown,
          vbfMass_JetFragmentationUp, vbfMass_JetPileUpDataMCDown, vbfMass_JetPileUpDataMCUp, vbfMass_JetPileUpPtBBDown, vbfMass_JetPileUpPtBBUp, vbfMass_JetPileUpPtEC1Down, vbfMass_JetPileUpPtEC1Up, vbfMass_JetPileUpPtEC2Down,
          vbfMass_JetPileUpPtEC2Up, vbfMass_JetPileUpPtHFDown, vbfMass_JetPileUpPtHFUp, vbfMass_JetPileUpPtRefDown, vbfMass_JetPileUpPtRefUp, vbfMass_JetRelativeBalDown, vbfMass_JetRelativeBalUp, vbfMass_JetRelativeFSRDown,
          vbfMass_JetRelativeFSRUp, vbfMass_JetRelativeJEREC1Down, vbfMass_JetRelativeJEREC1Up, vbfMass_JetRelativeJEREC2Down, vbfMass_JetRelativeJEREC2Up, vbfMass_JetRelativeJERHFDown, vbfMass_JetRelativeJERHFUp,
          vbfMass_JetRelativePtBBDown, vbfMass_JetRelativePtBBUp, vbfMass_JetRelativePtEC1Down, vbfMass_JetRelativePtEC1Up, vbfMass_JetRelativePtEC2Down, vbfMass_JetRelativePtEC2Up, vbfMass_JetRelativePtHFDown,
          vbfMass_JetRelativePtHFUp, vbfMass_JetRelativeSampleDown, vbfMass_JetRelativeSampleUp, vbfMass_JetRelativeStatECDown, vbfMass_JetRelativeStatECUp, vbfMass_JetRelativeStatFSRDown, vbfMass_JetRelativeStatFSRUp,
          vbfMass_JetRelativeStatHFDown, vbfMass_JetRelativeStatHFUp, vbfMass_JetSinglePionECALDown, vbfMass_JetSinglePionECALUp, vbfMass_JetSinglePionHCALDown, vbfMass_JetSinglePionHCALUp, vbfMass_JetTimePtEtaDown,
          vbfMass_JetTimePtEtaUp, vbfMass_JetTotalDown, vbfMass_JetTotalUp;

  Float_t met_UESUp, met_UESDown, met_JESUp, met_JESDown, metphi_UESUp, metphi_UESDown, metphi_JESUp, metphi_JESDown;

  Float_t jetVeto30_JetAbsoluteFlavMapDown, jetVeto30_JetAbsoluteFlavMapUp, jetVeto30_JetAbsoluteMPFBiasDown, jetVeto30_JetAbsoluteMPFBiasUp, jetVeto30_JetAbsoluteScaleDown, jetVeto30_JetAbsoluteScaleUp;

  // 2016 Placeholder
  Float_t amcatNLO_weight, eMatchesSingleE25Tight, eMatchesEle25TightFilter, singleE25eta2p1TightPass, eMatchesEle25eta2p1TightPath, tAgainstElectronVLooseMVA6, tAgainstMuonTight3, decayModeFindingNewDMs_2;


  // Member functions
  etau_tree2016(TTree *orig, TTree *itree, bool isMC, bool isEmbed, Int_t rec);
  virtual ~etau_tree2016() {}
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
etau_tree2016::etau_tree2016(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec) :
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
  original->SetBranchAddress("eMVANonTrigWP80", &eMVANonTrigWP80);

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
  original->SetBranchAddress("tByIsolationMVArun2v1DBoldDMwLTraw", &tByIsolationMVArun2v1DBoldDMwLTraw);
  original->SetBranchAddress("tByVLooseIsolationMVArun2v1DBoldDMwLT", &tByVLooseIsolationMVArun2v1DBoldDMwLT);
  original->SetBranchAddress("tRerunMVArun2v1DBoldDMwLTVLoose", &tRerunMVArun2v1DBoldDMwLTVLoose);

  // 2016 triggers
  original->SetBranchAddress("singleE25eta2p1TightPass", &singleE25eta2p1TightPass);
  original->SetBranchAddress("eMatchesSingleE25Tight", &eMatchesSingleE25Tight);
  original->SetBranchAddress("eMatchesEle25TightFilter", &eMatchesEle25TightFilter);
  original->SetBranchAddress("eMatchesEle25eta2p1TightPath", &eMatchesEle25eta2p1TightPath);

  // other
  original->SetBranchAddress("dielectronVeto", &dielectronVeto);
  original->SetBranchAddress("eVetoZTTp001dxyzR0", &eVetoZTTp001dxyzR0);
  original->SetBranchAddress("muVetoZTTp001dxyzR0", &muVetoZTTp001dxyzR0);
  original->SetBranchAddress("e_t_DR", &e_t_DR);
}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void etau_tree2016::do_skimming(TH1F* cutflow) {
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
    if (isMC) {
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          tau *= 0.982;
        } else if (tDecayMode == 1) {
          tau *= 1.010;
        } else if (tDecayMode == 10) {
          tau *= 1.004;
        }
      } else if (tZTTGenMatching == 1 || tZTTGenMatching == 3) {
        if (tDecayMode == 1) {
          tau *= 1.095;
        }
      } else if (tZTTGenMatching == 2 || tZTTGenMatching == 4) {
        if (tDecayMode == 0) {
          tau *= 0.998;
        } else if (tDecayMode == 1) {
          tau *= 1.015;
        }
      }
    }

    float el_pt_min(26), tau_pt_min;
    if (!isMC || tZTTGenMatching > 4) {
      tau_pt_min = 29.5;
    } else {
      tau_pt_min = 27.0;
    }

    cutflow->Fill(1., 1.);

    // apply event selection
    auto Ele25 = singleE25eta2p1TightPass && eMatchesEle25eta2p1TightPath && eMatchesEle25TightFilter;

    if (isEmbed || Ele25) cutflow->Fill(2., 1.);
    else  continue;

    if (!isEmbed || (eMatchesSingleE25Tight && singleE25eta2p1TightPass)) cutflow->Fill(3., 1.);
    else  continue;

    if (ePt > el_pt_min && fabs(eEta) < 2.1 && fabs(ePVDZ) < 0.2 && fabs(ePVDXY) < 0.045) cutflow->Fill(4., 1.);  // electron kinematic selection
    else  continue;

    if (eMVANonTrigWP80 && ePassesConversionVeto && eMissingHits < 2) cutflow->Fill(5., 1.);  // electron quality selection
    else  continue;

    if (tau.Pt() > tau_pt_min && fabs(tau.Eta()) < 2.3 && fabs(tPVDZ) < 0.2) cutflow->Fill(6., 1.);  // tau kinematic selection
    else  continue;

    if ((tByVLooseIsolationMVArun2v1DBoldDMwLT || tRerunMVArun2v1DBoldDMwLTVLoose) && tDecayModeFinding > 0 && fabs(tCharge) < 2) cutflow->Fill(7., 1.);  // tau quality selection
    else  continue;

    if (tAgainstMuonLoose3 > 0.5 && tAgainstElectronTightMVA6 > 0.5) cutflow->Fill(8., 1.);  // tau against leptons
    else  continue;

    if (muVetoZTTp001dxyzR0 == 0 && eVetoZTTp001dxyzR0 < 2 && dielectronVeto == 0) cutflow->Fill(9., 1.);  // vetos
    else  continue;

    if (e_t_DR > 0.5) cutflow->Fill(10., 1.);
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
      tauCandidate  = std::make_pair(tPt,  tByIsolationMVArun2v1DBoldDMwLTraw);
    } else {  // not a new event
      std::pair<float, float> currEleCandidate(ePt, eIsoDB03);
      std::pair<float, float> currTauCandidate(tPt, tByIsolationMVArun2v1DBoldDMwLTraw);

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
TTree* etau_tree2016::fill_tree(RecoilCorrector recoilPFMetCorrector) {
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
    nbtag = bjetCISVVeto20Medium;
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
    MET_UESUp.SetPtEtaPhiM(type1_pfMet_shiftedPt_UnclusteredEnUp, 0, type1_pfMet_shiftedPhi_UnclusteredEnUp, 0);
    MET_UESDown.SetPtEtaPhiM(type1_pfMet_shiftedPt_UnclusteredEnDown, 0, type1_pfMet_shiftedPhi_UnclusteredEnDown, 0);
    MET_JESUp.SetPtEtaPhiM(type1_pfMet_shiftedPt_JetEnUp, 0, type1_pfMet_shiftedPhi_JetEnUp, 0);
    MET_JESDown.SetPtEtaPhiM(type1_pfMet_shiftedPt_JetEnDown, 0, type1_pfMet_shiftedPhi_JetEnDown, 0);

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
          genpX,          // generator Z/W/Higgs px (float)
          genpY,          // generator Z/W/Higgs py (float)
          vispX,          // generator visible Z/W/Higgs px (float)
          vispY,          // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,  // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex,   // corrected type I pf met px (float)
          pfmetcorr_ey);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          genpX,                // generator Z/W/Higgs px (float)
          genpY,                // generator Z/W/Higgs py (float)
          vispX,                // generator visible Z/W/Higgs px (float)
          vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          genpX,                // generator Z/W/Higgs px (float)
          genpY,                // generator Z/W/Higgs py (float)
          vispX,                // generator visible Z/W/Higgs px (float)
          vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,        // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          genpX,                  // generator Z/W/Higgs px (float)
          genpY,                  // generator Z/W/Higgs py (float)
          vispX,                  // generator visible Z/W/Higgs px (float)
          vispY,                  // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          genpX,                  // generator Z/W/Higgs px (float)
          genpY,                  // generator Z/W/Higgs py (float)
          vispX,                  // generator visible Z/W/Higgs px (float)
          vispY,                  // generator visible Z/W/Higgs py (float)
          jetVeto30 + 1,          // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESDown);  // corrected type I pf met py (float)

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

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESUp.Px(),       // uncorrected type I pf met px (float)
          MET_JESUp.Py(),       // uncorrected type I pf met py (float)
          genpX,                // generator Z/W/Higgs px (float)
          genpY,                // generator Z/W/Higgs py (float)
          vispX,                // generator visible Z/W/Higgs px (float)
          vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESUp.Px(),       // uncorrected type I pf met px (float)
          MET_UESUp.Py(),       // uncorrected type I pf met py (float)
          genpX,                // generator Z/W/Higgs px (float)
          genpY,                // generator Z/W/Higgs py (float)
          vispX,                // generator visible Z/W/Higgs px (float)
          vispY,                // generator visible Z/W/Higgs py (float)
          jetVeto30,            // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESUp,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESUp);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_JESDown.Px(),       // uncorrected type I pf met px (float)
          MET_JESDown.Py(),       // uncorrected type I pf met py (float)
          genpX,                  // generator Z/W/Higgs px (float)
          genpY,                  // generator Z/W/Higgs py (float)
          vispX,                  // generator visible Z/W/Higgs px (float)
          vispY,                  // generator visible Z/W/Higgs py (float)
          jetVeto30,              // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_JESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_JESDown);  // corrected type I pf met py (float)

      recoilPFMetCorrector.CorrectByMeanResolution(
          MET_UESDown.Px(),       // uncorrected type I pf met px (float)
          MET_UESDown.Py(),       // uncorrected type I pf met py (float)
          genpX,                  // generator Z/W/Higgs px (float)
          genpY,                  // generator Z/W/Higgs py (float)
          vispX,                  // generator visible Z/W/Higgs px (float)
          vispY,                  // generator visible Z/W/Higgs py (float)
          jetVeto30,              // number of jets (hadronic jet multiplicity) (int)
          pfmetcorr_ex_UESDown,   // corrected type I pf met px (float)
          pfmetcorr_ey_UESDown);  // corrected type I pf met py (float)
    }

    MET.SetPxPyPzE(pfmetcorr_ex, pfmetcorr_ey, 0, sqrt(pfmetcorr_ex * pfmetcorr_ex + pfmetcorr_ey * pfmetcorr_ey));
    MET_UESUp.SetPxPyPzE(pfmetcorr_ex_UESUp, pfmetcorr_ey_UESUp, 0, sqrt(pfmetcorr_ex_UESUp * pfmetcorr_ex_UESUp + pfmetcorr_ey_UESUp * pfmetcorr_ey_UESUp));
    MET_UESDown.SetPxPyPzE(pfmetcorr_ex_UESDown, pfmetcorr_ey_UESDown, 0, sqrt(pfmetcorr_ex_UESDown * pfmetcorr_ex_UESDown + pfmetcorr_ey_UESDown * pfmetcorr_ey_UESDown));
    MET_JESUp.SetPxPyPzE(pfmetcorr_ex_JESUp, pfmetcorr_ey_JESUp, 0, sqrt(pfmetcorr_ex_JESUp * pfmetcorr_ex_JESUp + pfmetcorr_ey_JESUp * pfmetcorr_ey_JESUp));
    MET_JESDown.SetPxPyPzE(pfmetcorr_ex_JESDown, pfmetcorr_ey_JESDown, 0, sqrt(pfmetcorr_ex_JESDown * pfmetcorr_ex_JESDown + pfmetcorr_ey_JESDown * pfmetcorr_ey_JESDown));

    if (isMC) {
      // met correction due to tau energy scale
      if (tZTTGenMatching == 5) {
        if (tDecayMode == 0) {
          MET = MET + tau - 0.982*tau;
          MET_JESUp = MET_JESUp+tau-0.982*tau;
          MET_JESDown = MET_JESDown+tau-0.982*tau;
          MET_UESUp = MET_UESUp+tau-0.982*tau;
          MET_UESDown = MET_UESDown+tau-0.982*tau;
          tau *= 0.982;
        } else if (tDecayMode == 1) {
          MET = MET + tau - 1.010*tau;
          MET_JESUp = MET_JESUp+tau-1.010*tau;
          MET_JESDown = MET_JESDown+tau-1.010*tau;
          MET_UESUp = MET_UESUp+tau-1.010*tau;
          MET_UESDown = MET_UESDown+tau-1.010*tau;
          tau *= 1.010;
        } else if (tDecayMode == 10) {
          MET = MET + tau - 1.004*tau;
          MET_JESUp = MET_JESUp+tau-1.004*tau;
          MET_JESDown = MET_JESDown+tau-1.004*tau;
          MET_UESUp = MET_UESUp+tau-1.004*tau;
          MET_UESDown = MET_UESDown+tau-1.004*tau;
          tau *= 1.004;
        }
      } else if (tZTTGenMatching == 1 || tZTTGenMatching == 3) {
        if (tDecayMode == 1) {
          MET = MET+tau-1.095*tau;
          MET_JESUp = MET_JESUp+tau-1.095*tau;
          MET_JESDown = MET_JESDown+tau-1.095*tau;
          MET_UESUp = MET_UESUp+tau-1.095*tau;
          MET_UESDown = MET_UESDown+tau-1.095*tau;
          tau *= 1.095;
        }
      } else if (tZTTGenMatching == 2 || tZTTGenMatching == 4) {
        if (tDecayMode == 0) {
          MET = MET+tau-0.998*tau;
          MET_JESUp = MET_JESUp+tau-0.998*tau;
          MET_JESDown = MET_JESDown+tau-0.998*tau;
          MET_UESUp = MET_UESUp+tau-0.998*tau;
          MET_UESDown = MET_UESDown+tau-0.998*tau;
          tau *= 0.998;
        } else if (tDecayMode == 1) {
          MET = MET+tau-1.015*tau;
          MET_JESUp = MET_JESUp+tau-1.015*tau;
          MET_JESDown = MET_JESDown+tau-1.015*tau;
          MET_UESUp = MET_UESUp+tau-1.015*tau;
          MET_UESDown = MET_UESDown+tau-1.015*tau;
          tau *= 1.015;
        }
      }
    }

    met = MET.Pt();
    metphi = MET.Phi();
    met_px = MET.Px();
    met_py = MET.Py();

    met_JESUp   = MET_JESUp.Pt();
    met_JESDown = MET_JESDown.Pt();
    met_UESUp   = MET_UESUp.Pt();
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
void etau_tree2016::set_branches() {
  // output file branches
  tree->Branch("evt", &evt);
  tree->Branch("run" , &run);
  tree->Branch("lumi", &lumi);
  tree->Branch("gen_match_1", &gen_match_1, "gen_match_1/I");
  tree->Branch("gen_match_2", &gen_match_2, "gen_match_2/I");
  tree->Branch("njets", &njets, "njets/I");
  tree->Branch("nbtag", &nbtag, "nbtag/I");
  tree->Branch("njetspt20", &njetspt20, "njetspt20/I");
  tree->Branch("vbfMass", &vbfMass, "vbfMass/F");

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
  tree->Branch("eMatchesSingleE25Tight"  , &eMatchesSingleE25Tight  , "eMatchesSingleE25Tight/F");
  tree->Branch("eMatchesEle25TightFilter", &eMatchesEle25TightFilter, "eMatchesEle25TightFilter/F");
  tree->Branch("singleE25eta2p1TightPass", &singleE25eta2p1TightPass, "singleE25eta2p1TightPass/F");

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
  tree->Branch("iso_2", &tByIsolationMVArun2v1DBoldDMwLTraw, "iso_2/F");
  tree->Branch("decayModeFinding_2", &tDecayModeFinding, "decayModeFinding_2/F");
  tree->Branch("l2_decayMode", &tDecayMode, "l2_decayMode/F");

  tree->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_2"  , &tByLooseIsolationMVArun2v1DBoldDMwLT  , "byLooseIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_2" , &tByMediumIsolationMVArun2v1DBoldDMwLT , "byMediumIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2"  , &tByTightIsolationMVArun2v1DBoldDMwLT  , "byTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byVTightIsolationMVArun2v1DBoldDMwLT_2" , &tByVTightIsolationMVArun2v1DBoldDMwLT , "byVTightIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("byVVTightIsolationMVArun2v1DBoldDMwLT_2", &tByVVTightIsolationMVArun2v1DBoldDMwLT, "byVVTightIsolationMVArun2v1DBoldDMwLT_2/F");

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

  // Systematics
  tree->Branch("met_JESUp", &met_JESUp, "met_JESUp/F");
  tree->Branch("met_JESDown", &met_JESDown, "met_JESDown/F");
  tree->Branch("met_UESUp", &met_UESUp, "met_UESUp/F");
  tree->Branch("met_UESDown", &met_UESDown, "met_UESDown/F");
  tree->Branch("metphi_JESUp", &metphi_JESUp, "metphi_JESUp/F");
  tree->Branch("metphi_JESDown", &metphi_JESDown, "metphi_JESDown/F");
  tree->Branch("metphi_UESUp", &metphi_UESUp, "metphi_UESUp/F");
  tree->Branch("metphi_UESDown", &metphi_UESDown, "metphi_UESDown/F");

  tree->Branch("pfmetcorr_ex_UESUp"                       , &pfmetcorr_ex_UESUp                       , "pfmetcorr_ex_UESUp/F");
  tree->Branch("pfmetcorr_ey_UESUp"                       , &pfmetcorr_ey_UESUp                       , "pfmetcorr_ey_UESUp/F");
  tree->Branch("pfmetcorr_ex_UESDown"                     , &pfmetcorr_ex_UESDown                     , "pfmetcorr_ex_UESDown/F");
  tree->Branch("pfmetcorr_ey_UESDown"                     , &pfmetcorr_ey_UESDown                     , "pfmetcorr_ey_UESDown/F");
  tree->Branch("pfmetcorr_ex_JESUp"                       , &pfmetcorr_ex_JESUp                       , "pfmetcorr_ex_JESUp/F");
  tree->Branch("pfmetcorr_ey_JESUp"                       , &pfmetcorr_ey_JESUp                       , "pfmetcorr_ey_JESUp/F");
  tree->Branch("pfmetcorr_ex_JESDown"                     , &pfmetcorr_ex_JESDown                     , "pfmetcorr_ex_JESDown/F");
  tree->Branch("pfmetcorr_ey_JESDown"                     , &pfmetcorr_ey_JESDown                     , "pfmetcorr_ey_JESDown/F");
  tree->Branch("type1_pfMet_shiftedPt_UnclusteredEnUp"    , &type1_pfMet_shiftedPt_UnclusteredEnUp    , "type1_pfMet_shiftedPt_UnclusteredEnUp/F");
  tree->Branch("type1_pfMet_shiftedPhi_UnclusteredEnUp"   , &type1_pfMet_shiftedPhi_UnclusteredEnUp   , "type1_pfMet_shiftedPhi_UnclusteredEnUp/F");
  tree->Branch("type1_pfMet_shiftedPt_UnclusteredEnDown"  , &type1_pfMet_shiftedPt_UnclusteredEnDown  , "type1_pfMet_shiftedPt_UnclusteredEnDown/F");
  tree->Branch("type1_pfMet_shiftedPhi_UnclusteredEnDown" , &type1_pfMet_shiftedPhi_UnclusteredEnDown , "type1_pfMet_shiftedPhi_UnclusteredEnDown/F");
  tree->Branch("type1_pfMet_shiftedPt_JetEnUp"            , &type1_pfMet_shiftedPt_JetEnUp            , "type1_pfMet_shiftedPt_JetEnUp/F");
  tree->Branch("type1_pfMet_shiftedPhi_JetEnUp"           , &type1_pfMet_shiftedPhi_JetEnUp           , "type1_pfMet_shiftedPhi_JetEnUp/F");
  tree->Branch("type1_pfMet_shiftedPt_JetEnDown"          , &type1_pfMet_shiftedPt_JetEnDown          , "type1_pfMet_shiftedPt_JetEnDown/F");
  tree->Branch("type1_pfMet_shiftedPhi_JetEnDown"         , &type1_pfMet_shiftedPhi_JetEnDown         , "type1_pfMet_shiftedPhi_JetEnDown/F");

  tree->Branch("jetVeto30_JetEta0to3Down", &jetVeto30_JetEta0to3Down, "jetVeto30_JetEta0to3Down/F");
  tree->Branch("jetVeto30_JetEta0to3Up", &jetVeto30_JetEta0to3Up);
  tree->Branch("jetVeto30_JetEta0to5Down", &jetVeto30_JetEta0to5Down);
  tree->Branch("jetVeto30_JetEta0to5Up", &jetVeto30_JetEta0to5Up);
  tree->Branch("jetVeto30_JetEta3to5Down", &jetVeto30_JetEta3to5Down);
  tree->Branch("jetVeto30_JetEta3to5Up", &jetVeto30_JetEta3to5Up);
  tree->Branch("jetVeto30_JetRelativeBalDown", &jetVeto30_JetRelativeBalDown);
  tree->Branch("jetVeto30_JetRelativeBalUp", &jetVeto30_JetRelativeBalUp);
  tree->Branch("jetVeto30_JetRelativeSampleDown", &jetVeto30_JetRelativeSampleDown);
  tree->Branch("jetVeto30_JetRelativeSampleUp", &jetVeto30_JetRelativeSampleUp);
  tree->Branch("jetVeto30_JetTotalDown", &jetVeto30_JetTotalDown);
  tree->Branch("jetVeto30_JetTotalUp", &jetVeto30_JetTotalUp);
  tree->Branch("jetVeto30_JetAbsoluteFlavMapDown", &jetVeto30_JetAbsoluteFlavMapDown);
  tree->Branch("jetVeto30_JetAbsoluteFlavMapUp", &jetVeto30_JetAbsoluteFlavMapUp);
  tree->Branch("jetVeto30_JetAbsoluteMPFBiasDown", &jetVeto30_JetAbsoluteMPFBiasDown);
  tree->Branch("jetVeto30_JetAbsoluteMPFBiasUp", &jetVeto30_JetAbsoluteMPFBiasUp);
  tree->Branch("jetVeto30_JetAbsoluteScaleDown", &jetVeto30_JetAbsoluteScaleDown);
  tree->Branch("jetVeto30_JetAbsoluteScaleUp", &jetVeto30_JetAbsoluteScaleUp);
  tree->Branch("jetVeto30_JetAbsoluteStatDown", &jetVeto30_JetAbsoluteStatDown);
  tree->Branch("jetVeto30_JetAbsoluteStatUp", &jetVeto30_JetAbsoluteStatUp);
  tree->Branch("jetVeto30_JetClosureDown", &jetVeto30_JetClosureDown);
  tree->Branch("jetVeto30_JetClosureUp", &jetVeto30_JetClosureUp);
  tree->Branch("jetVeto30_JetEnDown", &jetVeto30_JetEnDown);
  tree->Branch("jetVeto30_JetFlavorQCDDown", &jetVeto30_JetFlavorQCDDown);
  tree->Branch("jetVeto30_JetFlavorQCDUp", &jetVeto30_JetFlavorQCDUp);
  tree->Branch("jetVeto30_JetFragmentationDown", &jetVeto30_JetFragmentationDown);
  tree->Branch("jetVeto30_JetFragmentationUp", &jetVeto30_JetFragmentationUp);
  tree->Branch("jetVeto30_JetPileUpDataMCDown", &jetVeto30_JetPileUpDataMCDown);
  tree->Branch("jetVeto30_JetPileUpDataMCUp", &jetVeto30_JetPileUpDataMCUp);
  tree->Branch("jetVeto30_JetPileUpPtBBDown", &jetVeto30_JetPileUpPtBBDown);
  tree->Branch("jetVeto30_JetPileUpPtBBUp", &jetVeto30_JetPileUpPtBBUp);
  tree->Branch("jetVeto30_JetPileUpPtEC1Down", &jetVeto30_JetPileUpPtEC1Down);
  tree->Branch("jetVeto30_JetPileUpPtEC1Up", &jetVeto30_JetPileUpPtEC1Up);
  tree->Branch("jetVeto30_JetPileUpPtEC2Down", &jetVeto30_JetPileUpPtEC2Down);
  tree->Branch("jetVeto30_JetPileUpPtEC2Up", &jetVeto30_JetPileUpPtEC2Up);
  tree->Branch("jetVeto30_JetPileUpPtHFDown", &jetVeto30_JetPileUpPtHFDown);
  tree->Branch("jetVeto30_JetPileUpPtHFUp", &jetVeto30_JetPileUpPtHFUp);
  tree->Branch("jetVeto30_JetPileUpPtRefDown", &jetVeto30_JetPileUpPtRefDown);
  tree->Branch("jetVeto30_JetPileUpPtRefUp", &jetVeto30_JetPileUpPtRefUp);
  tree->Branch("jetVeto30_JetRelativeBalDown", &jetVeto30_JetRelativeBalDown);
  tree->Branch("jetVeto30_JetRelativeBalUp", &jetVeto30_JetRelativeBalUp);
  tree->Branch("jetVeto30_JetRelativeFSRDown", &jetVeto30_JetRelativeFSRDown);
  tree->Branch("jetVeto30_JetRelativeFSRUp", &jetVeto30_JetRelativeFSRUp);
  tree->Branch("jetVeto30_JetRelativeJEREC1Down", &jetVeto30_JetRelativeJEREC1Down);
  tree->Branch("jetVeto30_JetRelativeJEREC1Up", &jetVeto30_JetRelativeJEREC1Up);
  tree->Branch("jetVeto30_JetRelativeJEREC2Down", &jetVeto30_JetRelativeJEREC2Down);
  tree->Branch("jetVeto30_JetRelativeJEREC2Up", &jetVeto30_JetRelativeJEREC2Up);
  tree->Branch("jetVeto30_JetRelativeJERHFDown", &jetVeto30_JetRelativeJERHFDown);
  tree->Branch("jetVeto30_JetRelativeJERHFUp", &jetVeto30_JetRelativeJERHFUp);
  tree->Branch("jetVeto30_JetRelativePtBBDown", &jetVeto30_JetRelativePtBBDown);
  tree->Branch("jetVeto30_JetRelativePtBBUp", &jetVeto30_JetRelativePtBBUp);
  tree->Branch("jetVeto30_JetRelativePtEC1Down", &jetVeto30_JetRelativePtEC1Down);
  tree->Branch("jetVeto30_JetRelativePtEC1Up", &jetVeto30_JetRelativePtEC1Up);
  tree->Branch("jetVeto30_JetRelativePtEC2Down", &jetVeto30_JetRelativePtEC2Down);
  tree->Branch("jetVeto30_JetRelativePtEC2Up", &jetVeto30_JetRelativePtEC2Up);
  tree->Branch("jetVeto30_JetRelativePtHFDown", &jetVeto30_JetRelativePtHFDown);
  tree->Branch("jetVeto30_JetRelativePtHFUp", &jetVeto30_JetRelativePtHFUp);
  tree->Branch("jetVeto30_JetRelativeSampleDown", &jetVeto30_JetRelativeSampleDown);
  tree->Branch("jetVeto30_JetRelativeSampleUp", &jetVeto30_JetRelativeSampleUp);
  tree->Branch("jetVeto30_JetRelativeStatECDown", &jetVeto30_JetRelativeStatECDown);
  tree->Branch("jetVeto30_JetRelativeStatECUp", &jetVeto30_JetRelativeStatECUp);
  tree->Branch("jetVeto30_JetRelativeStatFSRDown", &jetVeto30_JetRelativeStatFSRDown);
  tree->Branch("jetVeto30_JetRelativeStatFSRUp", &jetVeto30_JetRelativeStatFSRUp);
  tree->Branch("jetVeto30_JetRelativeStatHFDown", &jetVeto30_JetRelativeStatHFDown);
  tree->Branch("jetVeto30_JetRelativeStatHFUp", &jetVeto30_JetRelativeStatHFUp);
  tree->Branch("jetVeto30_JetSinglePionECALDown", &jetVeto30_JetSinglePionECALDown);
  tree->Branch("jetVeto30_JetSinglePionECALUp", &jetVeto30_JetSinglePionECALUp);
  tree->Branch("jetVeto30_JetSinglePionHCALDown", &jetVeto30_JetSinglePionHCALDown);
  tree->Branch("jetVeto30_JetSinglePionHCALUp", &jetVeto30_JetSinglePionHCALUp);
  tree->Branch("jetVeto30_JetTimePtEtaDown", &jetVeto30_JetTimePtEtaDown);
  tree->Branch("jetVeto30_JetTimePtEtaUp", &jetVeto30_JetTimePtEtaUp);
  tree->Branch("jetVeto30_JetTotalDown", &jetVeto30_JetTotalDown);
  tree->Branch("jetVeto30_JetTotalUp", &jetVeto30_JetTotalUp);

  tree->Branch("vbfMass_JetEta0to3Down", &vbfMass_JetEta0to3Down);
  tree->Branch("vbfMass_JetEta0to3Up", &vbfMass_JetEta0to3Up);
  tree->Branch("vbfMass_JetEta0to5Down", &vbfMass_JetEta0to5Down);
  tree->Branch("vbfMass_JetEta0to5Up", &vbfMass_JetEta0to5Up);
  tree->Branch("vbfMass_JetEta3to5Down", &vbfMass_JetEta3to5Down);
  tree->Branch("vbfMass_JetEta3to5Up", &vbfMass_JetEta3to5Up);
  tree->Branch("vbfMass_JetRelativeSampleDown", &vbfMass_JetRelativeSampleDown);
  tree->Branch("vbfMass_JetRelativeSampleUp", &vbfMass_JetRelativeSampleUp);
  tree->Branch("vbfMass_JetTotalDown", &vbfMass_JetTotalDown);
  tree->Branch("vbfMass_JetTotalUp", &vbfMass_JetTotalUp);
  tree->Branch("vbfMass_JetAbsoluteFlavMapDown", &vbfMass_JetAbsoluteFlavMapDown);
  tree->Branch("vbfMass_JetAbsoluteFlavMapUp", &vbfMass_JetAbsoluteFlavMapUp);
  tree->Branch("vbfMass_JetAbsoluteMPFBiasDown", &vbfMass_JetAbsoluteMPFBiasDown);
  tree->Branch("vbfMass_JetAbsoluteMPFBiasUp", &vbfMass_JetAbsoluteMPFBiasUp);
  tree->Branch("vbfMass_JetAbsoluteScaleDown", &vbfMass_JetAbsoluteScaleDown);
  tree->Branch("vbfMass_JetAbsoluteScaleUp", &vbfMass_JetAbsoluteScaleUp);
  tree->Branch("vbfMass_JetAbsoluteStatDown", &vbfMass_JetAbsoluteStatDown);
  tree->Branch("vbfMass_JetAbsoluteStatUp", &vbfMass_JetAbsoluteStatUp);
  tree->Branch("vbfMass_JetClosureDown", &vbfMass_JetClosureDown);
  tree->Branch("vbfMass_JetClosureUp", &vbfMass_JetClosureUp);
  tree->Branch("vbfMass_JetFlavorQCDDown", &vbfMass_JetFlavorQCDDown);
  tree->Branch("vbfMass_JetFlavorQCDUp", &vbfMass_JetFlavorQCDUp);
  tree->Branch("vbfMass_JetFragmentationDown", &vbfMass_JetFragmentationDown);
  tree->Branch("vbfMass_JetFragmentationUp", &vbfMass_JetFragmentationUp);
  tree->Branch("vbfMass_JetPileUpDataMCDown", &vbfMass_JetPileUpDataMCDown);
  tree->Branch("vbfMass_JetPileUpDataMCUp", &vbfMass_JetPileUpDataMCUp);
  tree->Branch("vbfMass_JetPileUpPtBBDown", &vbfMass_JetPileUpPtBBDown);
  tree->Branch("vbfMass_JetPileUpPtBBUp", &vbfMass_JetPileUpPtBBUp);
  tree->Branch("vbfMass_JetPileUpPtEC1Down", &vbfMass_JetPileUpPtEC1Down);
  tree->Branch("vbfMass_JetPileUpPtEC1Up", &vbfMass_JetPileUpPtEC1Up);
  tree->Branch("vbfMass_JetPileUpPtEC2Down", &vbfMass_JetPileUpPtEC2Down);
  tree->Branch("vbfMass_JetPileUpPtEC2Up", &vbfMass_JetPileUpPtEC2Up);
  tree->Branch("vbfMass_JetPileUpPtHFDown", &vbfMass_JetPileUpPtHFDown);
  tree->Branch("vbfMass_JetPileUpPtHFUp", &vbfMass_JetPileUpPtHFUp);
  tree->Branch("vbfMass_JetPileUpPtRefDown", &vbfMass_JetPileUpPtRefDown);
  tree->Branch("vbfMass_JetPileUpPtRefUp", &vbfMass_JetPileUpPtRefUp);
  tree->Branch("vbfMass_JetRelativeBalDown", &vbfMass_JetRelativeBalDown);
  tree->Branch("vbfMass_JetRelativeBalUp", &vbfMass_JetRelativeBalUp);
  tree->Branch("vbfMass_JetRelativeFSRDown", &vbfMass_JetRelativeFSRDown);
  tree->Branch("vbfMass_JetRelativeFSRUp", &vbfMass_JetRelativeFSRUp);
  tree->Branch("vbfMass_JetRelativeJEREC1Down", &vbfMass_JetRelativeJEREC1Down);
  tree->Branch("vbfMass_JetRelativeJEREC1Up", &vbfMass_JetRelativeJEREC1Up);
  tree->Branch("vbfMass_JetRelativeJEREC2Down", &vbfMass_JetRelativeJEREC2Down);
  tree->Branch("vbfMass_JetRelativeJEREC2Up", &vbfMass_JetRelativeJEREC2Up);
  tree->Branch("vbfMass_JetRelativeJERHFDown", &vbfMass_JetRelativeJERHFDown);
  tree->Branch("vbfMass_JetRelativeJERHFUp", &vbfMass_JetRelativeJERHFUp);
  tree->Branch("vbfMass_JetRelativePtBBDown", &vbfMass_JetRelativePtBBDown);
  tree->Branch("vbfMass_JetRelativePtBBUp", &vbfMass_JetRelativePtBBUp);
  tree->Branch("vbfMass_JetRelativePtEC1Down", &vbfMass_JetRelativePtEC1Down);
  tree->Branch("vbfMass_JetRelativePtEC1Up", &vbfMass_JetRelativePtEC1Up);
  tree->Branch("vbfMass_JetRelativePtEC2Down", &vbfMass_JetRelativePtEC2Down);
  tree->Branch("vbfMass_JetRelativePtEC2Up", &vbfMass_JetRelativePtEC2Up);
  tree->Branch("vbfMass_JetRelativePtHFDown", &vbfMass_JetRelativePtHFDown);
  tree->Branch("vbfMass_JetRelativePtHFUp", &vbfMass_JetRelativePtHFUp);
  tree->Branch("vbfMass_JetRelativeSampleDown", &vbfMass_JetRelativeSampleDown);
  tree->Branch("vbfMass_JetRelativeSampleUp", &vbfMass_JetRelativeSampleUp);
  tree->Branch("vbfMass_JetRelativeStatECDown", &vbfMass_JetRelativeStatECDown);
  tree->Branch("vbfMass_JetRelativeStatECUp", &vbfMass_JetRelativeStatECUp);
  tree->Branch("vbfMass_JetRelativeStatFSRDown", &vbfMass_JetRelativeStatFSRDown);
  tree->Branch("vbfMass_JetRelativeStatFSRUp", &vbfMass_JetRelativeStatFSRUp);
  tree->Branch("vbfMass_JetRelativeStatHFDown", &vbfMass_JetRelativeStatHFDown);
  tree->Branch("vbfMass_JetRelativeStatHFUp", &vbfMass_JetRelativeStatHFUp);
  tree->Branch("vbfMass_JetSinglePionECALDown", &vbfMass_JetSinglePionECALDown);
  tree->Branch("vbfMass_JetSinglePionECALUp", &vbfMass_JetSinglePionECALUp);
  tree->Branch("vbfMass_JetSinglePionHCALDown", &vbfMass_JetSinglePionHCALDown);
  tree->Branch("vbfMass_JetSinglePionHCALUp", &vbfMass_JetSinglePionHCALUp);
  tree->Branch("vbfMass_JetTimePtEtaDown", &vbfMass_JetTimePtEtaDown);
  tree->Branch("vbfMass_JetTimePtEtaUp", &vbfMass_JetTimePtEtaUp);
  tree->Branch("vbfMass_JetTotalDown", &vbfMass_JetTotalDown);
  tree->Branch("vbfMass_JetTotalUp", &vbfMass_JetTotalUp);

  // 2016 placeholders

  tree->Branch("againstElectronTightMVA6_2", &tAgainstElectronTightMVA6, "againstElectronTightMVA6_2/F");
  tree->Branch("againstElectronVLooseMVA6_2", &tAgainstElectronVLooseMVA6, "againstElectronVLooseMVA6_2/F");
  tree->Branch("againstMuonTight3_2", &tAgainstMuonTight3, "againstMuonTight3_2/F");
  tree->Branch("againstMuonLoose3_2", &tAgainstMuonLoose3, "againstMuonLoose3_2/F");
  tree->Branch("byVLooseIsolationMVArun2v1DBoldDMwLT_2", &tByVLooseIsolationMVArun2v1DBoldDMwLT, "byVLooseIsolationMVArun2v1DBoldDMwLT_2/F");
  tree->Branch("decayModeFindingNewDMs_2", &decayModeFindingNewDMs_2, "decayModeFindingNewDMs_2/F");


  // input branches
  // event info
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
  original->SetBranchAddress("vbfMass", &vbfMass);

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

  // jet branches
  original->SetBranchAddress("jetVeto20", &jetVeto20);
  original->SetBranchAddress("jetVeto30", &jetVeto30);

  original->SetBranchAddress("j1pt", &jpt_1);
  original->SetBranchAddress("j1eta", &jeta_1);
  original->SetBranchAddress("j1phi", &jphi_1);
  original->SetBranchAddress("j1csv", &jcsv_1);
  original->SetBranchAddress("j2pt", &jpt_2);
  original->SetBranchAddress("j2eta", &jeta_2);
  original->SetBranchAddress("j2phi", &jphi_2);
  original->SetBranchAddress("j2csv", &jcsv_2);
  original->SetBranchAddress("jb1pt",  &bpt_1);
  original->SetBranchAddress("jb1eta", &beta_1);
  original->SetBranchAddress("jb1phi", &bphi_1);
  original->SetBranchAddress("jb1csv", &bcsv_1);
  original->SetBranchAddress("jb1hadronflavor", &bflavor_1);
  original->SetBranchAddress("jb2pt",  &bpt_2);
  original->SetBranchAddress("jb2eta", &beta_2);
  original->SetBranchAddress("jb2phi", &bphi_2);
  original->SetBranchAddress("jb2csv", &bcsv_2);
  original->SetBranchAddress("jb2hadronflavor", &bflavor_2);

  original->SetBranchAddress("bjetCISVVeto20Medium", &bjetCISVVeto20Medium);
  original->SetBranchAddress("bjetCISVVeto20Tight", &bjetCISVVeto20Tight);
  original->SetBranchAddress("bjetCISVVeto30Loose", &bjetCISVVeto30Loose);
  original->SetBranchAddress("bjetCISVVeto30Medium", &bjetCISVVeto30Medium);
  original->SetBranchAddress("bjetCISVVeto30Tight", &bjetCISVVeto30Tight);

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

  // Systematics
  original->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnUp"    , &type1_pfMet_shiftedPt_UnclusteredEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnUp"   , &type1_pfMet_shiftedPhi_UnclusteredEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_UnclusteredEnDown"  , &type1_pfMet_shiftedPt_UnclusteredEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_UnclusteredEnDown" , &type1_pfMet_shiftedPhi_UnclusteredEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPt_JetEnUp"            , &type1_pfMet_shiftedPt_JetEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnUp"           , &type1_pfMet_shiftedPhi_JetEnUp);
  original->SetBranchAddress("type1_pfMet_shiftedPt_JetEnDown"          , &type1_pfMet_shiftedPt_JetEnDown);
  original->SetBranchAddress("type1_pfMet_shiftedPhi_JetEnDown"         , &type1_pfMet_shiftedPhi_JetEnDown);

  original->SetBranchAddress("jetVeto30_JetEta0to3Down", &jetVeto30_JetEta0to3Down);
  original->SetBranchAddress("jetVeto30_JetEta0to3Up", &jetVeto30_JetEta0to3Up);
  original->SetBranchAddress("jetVeto30_JetEta0to5Down", &jetVeto30_JetEta0to5Down);
  original->SetBranchAddress("jetVeto30_JetEta0to5Up", &jetVeto30_JetEta0to5Up);
  original->SetBranchAddress("jetVeto30_JetEta3to5Down", &jetVeto30_JetEta3to5Down);
  original->SetBranchAddress("jetVeto30_JetEta3to5Up", &jetVeto30_JetEta3to5Up);
  original->SetBranchAddress("jetVeto30_JetRelativeSampleDown", &jetVeto30_JetRelativeSampleDown);
  original->SetBranchAddress("jetVeto30_JetRelativeSampleUp", &jetVeto30_JetRelativeSampleUp);
  original->SetBranchAddress("jetVeto30_JetTotalDown", &jetVeto30_JetTotalDown);
  original->SetBranchAddress("jetVeto30_JetTotalUp", &jetVeto30_JetTotalUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteFlavMapDown", &jetVeto30_JetAbsoluteFlavMapDown);
  original->SetBranchAddress("jetVeto30_JetAbsoluteFlavMapUp", &jetVeto30_JetAbsoluteFlavMapUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteMPFBiasDown", &jetVeto30_JetAbsoluteMPFBiasDown);
  original->SetBranchAddress("jetVeto30_JetAbsoluteMPFBiasUp", &jetVeto30_JetAbsoluteMPFBiasUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteScaleDown", &jetVeto30_JetAbsoluteScaleDown);
  original->SetBranchAddress("jetVeto30_JetAbsoluteScaleUp", &jetVeto30_JetAbsoluteScaleUp);
  original->SetBranchAddress("jetVeto30_JetAbsoluteStatDown", &jetVeto30_JetAbsoluteStatDown);
  original->SetBranchAddress("jetVeto30_JetAbsoluteStatUp", &jetVeto30_JetAbsoluteStatUp);
  original->SetBranchAddress("jetVeto30_JetClosureDown", &jetVeto30_JetClosureDown);
  original->SetBranchAddress("jetVeto30_JetClosureUp", &jetVeto30_JetClosureUp);
  original->SetBranchAddress("jetVeto30_JetFlavorQCDDown", &jetVeto30_JetFlavorQCDDown);
  original->SetBranchAddress("jetVeto30_JetFlavorQCDUp", &jetVeto30_JetFlavorQCDUp);
  original->SetBranchAddress("jetVeto30_JetFragmentationDown", &jetVeto30_JetFragmentationDown);
  original->SetBranchAddress("jetVeto30_JetFragmentationUp", &jetVeto30_JetFragmentationUp);
  original->SetBranchAddress("jetVeto30_JetPileUpDataMCDown", &jetVeto30_JetPileUpDataMCDown);
  original->SetBranchAddress("jetVeto30_JetPileUpDataMCUp", &jetVeto30_JetPileUpDataMCUp);
  original->SetBranchAddress("jetVeto30_JetPileUpPtBBDown", &jetVeto30_JetPileUpPtBBDown);
  original->SetBranchAddress("jetVeto30_JetPileUpPtBBUp", &jetVeto30_JetPileUpPtBBUp);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC1Down", &jetVeto30_JetPileUpPtEC1Down);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC1Up", &jetVeto30_JetPileUpPtEC1Up);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC2Down", &jetVeto30_JetPileUpPtEC2Down);
  original->SetBranchAddress("jetVeto30_JetPileUpPtEC2Up", &jetVeto30_JetPileUpPtEC2Up);
  original->SetBranchAddress("jetVeto30_JetPileUpPtHFDown", &jetVeto30_JetPileUpPtHFDown);
  original->SetBranchAddress("jetVeto30_JetPileUpPtHFUp", &jetVeto30_JetPileUpPtHFUp);
  original->SetBranchAddress("jetVeto30_JetPileUpPtRefDown", &jetVeto30_JetPileUpPtRefDown);
  original->SetBranchAddress("jetVeto30_JetPileUpPtRefUp", &jetVeto30_JetPileUpPtRefUp);
  original->SetBranchAddress("jetVeto30_JetRelativeBalDown", &jetVeto30_JetRelativeBalDown);
  original->SetBranchAddress("jetVeto30_JetRelativeBalUp", &jetVeto30_JetRelativeBalUp);
  original->SetBranchAddress("jetVeto30_JetRelativeFSRDown", &jetVeto30_JetRelativeFSRDown);
  original->SetBranchAddress("jetVeto30_JetRelativeFSRUp", &jetVeto30_JetRelativeFSRUp);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC1Down", &jetVeto30_JetRelativeJEREC1Down);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC1Up", &jetVeto30_JetRelativeJEREC1Up);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC2Down", &jetVeto30_JetRelativeJEREC2Down);
  original->SetBranchAddress("jetVeto30_JetRelativeJEREC2Up", &jetVeto30_JetRelativeJEREC2Up);
  original->SetBranchAddress("jetVeto30_JetRelativeJERHFDown", &jetVeto30_JetRelativeJERHFDown);
  original->SetBranchAddress("jetVeto30_JetRelativeJERHFUp", &jetVeto30_JetRelativeJERHFUp);
  original->SetBranchAddress("jetVeto30_JetRelativePtBBDown", &jetVeto30_JetRelativePtBBDown);
  original->SetBranchAddress("jetVeto30_JetRelativePtBBUp", &jetVeto30_JetRelativePtBBUp);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC1Down", &jetVeto30_JetRelativePtEC1Down);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC1Up", &jetVeto30_JetRelativePtEC1Up);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC2Down", &jetVeto30_JetRelativePtEC2Down);
  original->SetBranchAddress("jetVeto30_JetRelativePtEC2Up", &jetVeto30_JetRelativePtEC2Up);
  original->SetBranchAddress("jetVeto30_JetRelativePtHFDown", &jetVeto30_JetRelativePtHFDown);
  original->SetBranchAddress("jetVeto30_JetRelativePtHFUp", &jetVeto30_JetRelativePtHFUp);
  original->SetBranchAddress("jetVeto30_JetRelativeSampleDown", &jetVeto30_JetRelativeSampleDown);
  original->SetBranchAddress("jetVeto30_JetRelativeSampleUp", &jetVeto30_JetRelativeSampleUp);
  original->SetBranchAddress("jetVeto30_JetRelativeStatECDown", &jetVeto30_JetRelativeStatECDown);
  original->SetBranchAddress("jetVeto30_JetRelativeStatECUp", &jetVeto30_JetRelativeStatECUp);
  original->SetBranchAddress("jetVeto30_JetRelativeStatFSRDown", &jetVeto30_JetRelativeStatFSRDown);
  original->SetBranchAddress("jetVeto30_JetRelativeStatFSRUp", &jetVeto30_JetRelativeStatFSRUp);
  original->SetBranchAddress("jetVeto30_JetRelativeStatHFDown", &jetVeto30_JetRelativeStatHFDown);
  original->SetBranchAddress("jetVeto30_JetRelativeStatHFUp", &jetVeto30_JetRelativeStatHFUp);
  original->SetBranchAddress("jetVeto30_JetSinglePionECALDown", &jetVeto30_JetSinglePionECALDown);
  original->SetBranchAddress("jetVeto30_JetSinglePionECALUp", &jetVeto30_JetSinglePionECALUp);
  original->SetBranchAddress("jetVeto30_JetSinglePionHCALDown", &jetVeto30_JetSinglePionHCALDown);
  original->SetBranchAddress("jetVeto30_JetSinglePionHCALUp", &jetVeto30_JetSinglePionHCALUp);
  original->SetBranchAddress("jetVeto30_JetTimePtEtaDown", &jetVeto30_JetTimePtEtaDown);
  original->SetBranchAddress("jetVeto30_JetTimePtEtaUp", &jetVeto30_JetTimePtEtaUp);
  original->SetBranchAddress("jetVeto30_JetTotalDown", &jetVeto30_JetTotalDown);
  original->SetBranchAddress("jetVeto30_JetTotalUp", &jetVeto30_JetTotalUp);



  original->SetBranchAddress("vbfMass_JetEta0to3Down", &vbfMass_JetEta0to3Down);
  original->SetBranchAddress("vbfMass_JetEta0to3Up", &vbfMass_JetEta0to3Up);
  original->SetBranchAddress("vbfMass_JetEta0to5Down", &vbfMass_JetEta0to5Down);
  original->SetBranchAddress("vbfMass_JetEta0to5Up", &vbfMass_JetEta0to5Up);
  original->SetBranchAddress("vbfMass_JetEta3to5Down", &vbfMass_JetEta3to5Down);
  original->SetBranchAddress("vbfMass_JetEta3to5Up", &vbfMass_JetEta3to5Up);
  original->SetBranchAddress("vbfMass_JetRelativeSampleDown", &vbfMass_JetRelativeSampleDown);
  original->SetBranchAddress("vbfMass_JetRelativeSampleUp", &vbfMass_JetRelativeSampleUp);
  original->SetBranchAddress("vbfMass_JetTotalDown", &vbfMass_JetTotalDown);
  original->SetBranchAddress("vbfMass_JetTotalUp", &vbfMass_JetTotalUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteFlavMapDown", &vbfMass_JetAbsoluteFlavMapDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteFlavMapUp", &vbfMass_JetAbsoluteFlavMapUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteMPFBiasDown", &vbfMass_JetAbsoluteMPFBiasDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteMPFBiasUp", &vbfMass_JetAbsoluteMPFBiasUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteScaleDown", &vbfMass_JetAbsoluteScaleDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteScaleUp", &vbfMass_JetAbsoluteScaleUp);
  original->SetBranchAddress("vbfMass_JetAbsoluteStatDown", &vbfMass_JetAbsoluteStatDown);
  original->SetBranchAddress("vbfMass_JetAbsoluteStatUp", &vbfMass_JetAbsoluteStatUp);
  original->SetBranchAddress("vbfMass_JetClosureDown", &vbfMass_JetClosureDown);
  original->SetBranchAddress("vbfMass_JetClosureUp", &vbfMass_JetClosureUp);
  original->SetBranchAddress("vbfMass_JetFlavorQCDDown", &vbfMass_JetFlavorQCDDown);
  original->SetBranchAddress("vbfMass_JetFlavorQCDUp", &vbfMass_JetFlavorQCDUp);
  original->SetBranchAddress("vbfMass_JetFragmentationDown", &vbfMass_JetFragmentationDown);
  original->SetBranchAddress("vbfMass_JetFragmentationUp", &vbfMass_JetFragmentationUp);
  original->SetBranchAddress("vbfMass_JetPileUpDataMCDown", &vbfMass_JetPileUpDataMCDown);
  original->SetBranchAddress("vbfMass_JetPileUpDataMCUp", &vbfMass_JetPileUpDataMCUp);
  original->SetBranchAddress("vbfMass_JetPileUpPtBBDown", &vbfMass_JetPileUpPtBBDown);
  original->SetBranchAddress("vbfMass_JetPileUpPtBBUp", &vbfMass_JetPileUpPtBBUp);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC1Down", &vbfMass_JetPileUpPtEC1Down);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC1Up", &vbfMass_JetPileUpPtEC1Up);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC2Down", &vbfMass_JetPileUpPtEC2Down);
  original->SetBranchAddress("vbfMass_JetPileUpPtEC2Up", &vbfMass_JetPileUpPtEC2Up);
  original->SetBranchAddress("vbfMass_JetPileUpPtHFDown", &vbfMass_JetPileUpPtHFDown);
  original->SetBranchAddress("vbfMass_JetPileUpPtHFUp", &vbfMass_JetPileUpPtHFUp);
  original->SetBranchAddress("vbfMass_JetPileUpPtRefDown", &vbfMass_JetPileUpPtRefDown);
  original->SetBranchAddress("vbfMass_JetPileUpPtRefUp", &vbfMass_JetPileUpPtRefUp);
  original->SetBranchAddress("vbfMass_JetRelativeBalDown", &vbfMass_JetRelativeBalDown);
  original->SetBranchAddress("vbfMass_JetRelativeBalUp", &vbfMass_JetRelativeBalUp);
  original->SetBranchAddress("vbfMass_JetRelativeFSRDown", &vbfMass_JetRelativeFSRDown);
  original->SetBranchAddress("vbfMass_JetRelativeFSRUp", &vbfMass_JetRelativeFSRUp);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC1Down", &vbfMass_JetRelativeJEREC1Down);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC1Up", &vbfMass_JetRelativeJEREC1Up);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC2Down", &vbfMass_JetRelativeJEREC2Down);
  original->SetBranchAddress("vbfMass_JetRelativeJEREC2Up", &vbfMass_JetRelativeJEREC2Up);
  original->SetBranchAddress("vbfMass_JetRelativeJERHFDown", &vbfMass_JetRelativeJERHFDown);
  original->SetBranchAddress("vbfMass_JetRelativeJERHFUp", &vbfMass_JetRelativeJERHFUp);
  original->SetBranchAddress("vbfMass_JetRelativePtBBDown", &vbfMass_JetRelativePtBBDown);
  original->SetBranchAddress("vbfMass_JetRelativePtBBUp", &vbfMass_JetRelativePtBBUp);
  original->SetBranchAddress("vbfMass_JetRelativePtEC1Down", &vbfMass_JetRelativePtEC1Down);
  original->SetBranchAddress("vbfMass_JetRelativePtEC1Up", &vbfMass_JetRelativePtEC1Up);
  original->SetBranchAddress("vbfMass_JetRelativePtEC2Down", &vbfMass_JetRelativePtEC2Down);
  original->SetBranchAddress("vbfMass_JetRelativePtEC2Up", &vbfMass_JetRelativePtEC2Up);
  original->SetBranchAddress("vbfMass_JetRelativePtHFDown", &vbfMass_JetRelativePtHFDown);
  original->SetBranchAddress("vbfMass_JetRelativePtHFUp", &vbfMass_JetRelativePtHFUp);
  original->SetBranchAddress("vbfMass_JetRelativeSampleDown", &vbfMass_JetRelativeSampleDown);
  original->SetBranchAddress("vbfMass_JetRelativeSampleUp", &vbfMass_JetRelativeSampleUp);
  original->SetBranchAddress("vbfMass_JetRelativeStatECDown", &vbfMass_JetRelativeStatECDown);
  original->SetBranchAddress("vbfMass_JetRelativeStatECUp", &vbfMass_JetRelativeStatECUp);
  original->SetBranchAddress("vbfMass_JetRelativeStatFSRDown", &vbfMass_JetRelativeStatFSRDown);
  original->SetBranchAddress("vbfMass_JetRelativeStatFSRUp", &vbfMass_JetRelativeStatFSRUp);
  original->SetBranchAddress("vbfMass_JetRelativeStatHFDown", &vbfMass_JetRelativeStatHFDown);
  original->SetBranchAddress("vbfMass_JetRelativeStatHFUp", &vbfMass_JetRelativeStatHFUp);
  original->SetBranchAddress("vbfMass_JetSinglePionECALDown", &vbfMass_JetSinglePionECALDown);
  original->SetBranchAddress("vbfMass_JetSinglePionECALUp", &vbfMass_JetSinglePionECALUp);
  original->SetBranchAddress("vbfMass_JetSinglePionHCALDown", &vbfMass_JetSinglePionHCALDown);
  original->SetBranchAddress("vbfMass_JetSinglePionHCALUp", &vbfMass_JetSinglePionHCALUp);
  original->SetBranchAddress("vbfMass_JetTimePtEtaDown", &vbfMass_JetTimePtEtaDown);
  original->SetBranchAddress("vbfMass_JetTimePtEtaUp", &vbfMass_JetTimePtEtaUp);
  original->SetBranchAddress("vbfMass_JetTotalDown", &vbfMass_JetTotalDown);
  original->SetBranchAddress("vbfMass_JetTotalUp", &vbfMass_JetTotalUp);
}
