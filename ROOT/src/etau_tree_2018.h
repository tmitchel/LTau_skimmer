// Copyright 2019 Tyler Mitchell

#ifndef ROOT_SRC_ETAU_TREE_2018_H_
#define ROOT_SRC_ETAU_TREE_2018_H_

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "./base_tree.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "ltau_skimmer/ROOT/interface/etau_input_branches.h"
#include "ltau_skimmer/ROOT/interface/TauFESTool.h"

class etau_tree2018 : public virtual base_tree {
   private:
    TTree *tree, *original;
    etau_input_branches* in;
    bool isMC, isEmbed, isSignal;
    std::vector<Int_t> good_events;
    TLorentzVector ele, tau, MET, MET_reso_Up, MET_reso_Down, MET_resp_Up, MET_resp_Down;
    TLorentzVector MET_JERUp, MET_AbsoluteUp, MET_AbsoluteyearUp, MET_BBEC1Up, MET_BBEC1yearUp, MET_EC2Up, MET_EC2yearUp, MET_EnUp, MET_FlavorQCDUp,
        MET_HFUp, MET_HFyearUp, MET_RelBalUp, MET_RelSamUp, MET_ResUp, MET_TotalUp, MET_UESUp,
        MET_JERDown, MET_AbsoluteDown, MET_AbsoluteyearDown, MET_BBEC1Down, MET_BBEC1yearDown, MET_EC2Down, MET_EC2yearDown, MET_EnDown,
        MET_FlavorQCDDown, MET_HFDown, MET_HFyearDown, MET_RelBalDown, MET_RelSamDown, MET_ResDown, MET_TotalDown, MET_UESDown;
    TauFESTool tfes;

   public:
    // Member variables
    UInt_t Run, Lumi;
    Int_t recoil;
    Float_t placeholder;  // for all branches not present in 2018

    // // Constructed while running
    Int_t era;
    Int_t gen_match_1, gen_match_2, njets, nbtag, njetspt20;
    Float_t tes_dm0_sf, tes_dm1_sf, tes_dm10_sf, efake_dm0_sf, efake_dm1_sf, mfake_dm0_sf, mfake_dm1_sf;
    Float_t jetVeto20, jetVeto30, met, metphi, met_px, met_py, extraelec_veto, extramuon_veto, dilepton_veto, pfmetcorr_ex, pfmetcorr_ey;
    Float_t met_reso_Up, met_reso_Down, met_resp_Up, met_resp_Down, metphi_reso_Up, metphi_reso_Down, metphi_resp_Up, metphi_resp_Down;
    Float_t met_JERUp, met_AbsoluteUp, met_AbsoluteyearUp, met_BBEC1Up, met_BBEC1yearUp, met_EC2Up, met_EC2yearUp, met_EnUp, met_FlavorQCDUp,
        met_HFUp, met_HFyearUp, met_RelBalUp, met_RelSamUp, met_ResUp, met_TotalUp, met_UESUp,
        met_JERDown, met_AbsoluteDown, met_AbsoluteyearDown, met_BBEC1Down, met_BBEC1yearDown, met_EC2Down, met_EC2yearDown, met_EnDown,
        met_FlavorQCDDown, met_HFDown, met_HFyearDown, met_RelBalDown, met_RelSamDown, met_ResDown, met_TotalDown, met_UESDown;
    Float_t metphi_JERUp, metphi_AbsoluteUp, metphi_AbsoluteyearUp, metphi_BBEC1Up, metphi_BBEC1yearUp, metphi_EC2Up, metphi_EC2yearUp, metphi_EnUp,
        metphi_FlavorQCDUp, metphi_HFUp, metphi_HFyearUp, metphi_RelBalUp, metphi_RelSamUp, metphi_ResUp, metphi_TotalUp, metphi_UESUp,
        metphi_JERDown, metphi_AbsoluteDown, metphi_AbsoluteyearDown, metphi_BBEC1Down, metphi_BBEC1yearDown, metphi_EC2Down, metphi_EC2yearDown,
        metphi_EnDown, metphi_FlavorQCDDown, metphi_HFDown, metphi_HFyearDown, metphi_RelBalDown, metphi_RelSamDown, metphi_ResDown, metphi_TotalDown,
        metphi_UESDown;
    Float_t tes_syst_up, tes_syst_down, ftes_syst_up, ftes_syst_down;
    Float_t pt_1, eta_1, phi_1, m_1, e_1, px_1, py_1, pz_1, pt_2, eta_2, phi_2, m_2, e_2, px_2, py_2, pz_2;

    std::vector<TLorentzVector*> mets;

    // Member functions
    etau_tree2018(TTree* orig, TTree* itree, bool isMC, bool isEmbed, bool IsSignal, Int_t rec);
    virtual ~etau_tree2018() {}
    void do_skimming(TH1F*);
    void set_branches();
    void do_met_corr_nom(Float_t, TLorentzVector, TLorentzVector*);
    void do_recoil_corr(RecoilCorrector*, TLorentzVector*, int);
    TTree* fill_tree(RecoilCorrector, MEtSys);
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
etau_tree2018::etau_tree2018(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, bool IsSignal, Int_t rec)
    : tree(itree),
      original(Original),
      in(new etau_input_branches(Original)),
      isMC(IsMC),
      isEmbed(IsEmbed),
      isSignal(IsSignal),
      tfes("2018ReReco", "DeepTau2017v2p1", "ltau_skimmer/ROOT/data/", "et", isEmbed),
      recoil(rec),
      era(2018) {}

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

    bool isData = !isEmbed && !isMC;
    Int_t nevt = (Int_t)original->GetEntries();
    for (auto ievt = 0; ievt < nevt; ievt++) {
        original->GetEntry(ievt);
        evt_now = in->evt;

        // TLorentzVector ele, tau;
        ele.SetPtEtaPhiM(in->ePt, in->eEta, in->ePhi, in->eMass);
        tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);

        // electron energy scale
        if (isMC) {
            ele *= in->eCorrectedEt / ele.Energy();
        }

        // apply TES
        if (isMC || isEmbed) {
            tau *= tfes.getFES(in->tDecayMode, in->tEta, in->tZTTGenMatching);
            tau *= tfes.getTES(in->tPt, in->tDecayMode, in->tZTTGenMatching);
        }

        cutflow->Fill(1., 1.);
        // apply event selection

        auto Ele32 = in->eMatchesEle32Filter && in->eMatchesEle32Path && in->Ele32WPTightPass;
        auto Ele35 = in->eMatchesEle35Filter && in->eMatchesEle35Path && in->Ele35WPTightPass;
        auto Cross_v1 = in->run < 317509 && isData && in->Ele24LooseTau30Pass &&
                        in->eMatchesEle24Tau30Filter && in->eMatchesEle24Tau30Path && in->tMatchesEle24Tau30Filter && in->tMatchesEle24Tau30Path;
        auto Cross_v2 = (isMC || (in->run >= 317509 && isData))&& in->Ele24LooseHPSTau30TightIDPass && 
                        in->eMatchesEle24HPSTau30Filter && in->eMatchesEle24HPSTau30Path && in->tMatchesEle24HPSTau30Filter && in->tMatchesEle24HPSTau30Path;
        auto Ele32_emb = in->eMatchEmbeddedFilterEle32 || fabs(ele.Eta()) > 1.479;
        auto Ele35_emb = in->eMatchEmbeddedFilterEle35 || fabs(ele.Eta()) > 1.479;
        // auto Cross_emb = (in->eMatchEmbeddedFilterEle24Tau30 && in->tMatchEmbeddedFilterEle24Tau30) || fabs(ele.Eta()) > 1.479;

        // embedded has it's own trigger paths
        if (isEmbed) {
            if (Ele35 && ele.Pt() > 36) {
                cutflow->Fill(2., 1.);
            } else if (Ele32 && ele.Pt() > 33) {
                cutflow->Fill(2., 1.);
            } else if (ele.Pt() > 25 && ele.Pt() < 33 && tau.Pt() > 35 && fabs(tau.Eta()) < 2.1) {
                cutflow->Fill(2., 1.);
            } else {
                continue;
            }
        } else {
            if (Ele35 && ele.Pt() > 36) {
                cutflow->Fill(2., 1.);
            } else if (Ele32 && ele.Pt() > 33) {
                cutflow->Fill(2., 1.);
            } else if ((Cross_v1 || Cross_v2) && ele.Pt() > 25 && ele.Pt() < 33 && tau.Pt() > 35 && fabs(tau.Eta()) < 2.1) {
                cutflow->Fill(2., 1.);
            } else {
                continue;
            }
        }

        if (ele.Pt() > 25. && fabs(ele.Eta()) < 2.1 && fabs(in->ePVDZ) < 0.2 && fabs(in->ePVDXY) < 0.045)
            cutflow->Fill(4., 1.);  // electron kinematic selection
        else
            continue;

        if (in->eMVANoisoWP90 && in->ePassesConversionVeto && in->eMissingHits < 2)
            cutflow->Fill(5., 1.);  // electron quality selection
        else
            continue;

        if (tau.Pt() > 30. && fabs(tau.Eta()) < 2.3 && fabs(in->tPVDZ) < 0.2)
            cutflow->Fill(6., 1.);  // tau kinematic selection
        else
            continue;

        if (in->tVVVLooseDeepTau2017v2p1VSjet > 0.5 && in->tDecayModeFindingNewDMs &&
            in->tDecayMode != 5 && in->tDecayMode != 6  && fabs(in->tCharge) < 2)
            cutflow->Fill(7., 1.);  // tau quality selection
        else
            continue;

        if (in->tVVVLooseDeepTau2017v2p1VSmu > 0.5 && in->tVVVLooseDeepTau2017v2p1VSe > 0.5)
            cutflow->Fill(8., 1.);  // tau against leptons
        else
            continue;

        if (in->muVetoZTTp001dxyzR0 == 0 && in->eVetoZTTp001dxyzR0 < 2 && in->dielectronVeto == 0)
            cutflow->Fill(9., 1.);  // vetos
        else
            continue;

        if (ele.DeltaR(tau) > 0.5) {
            cutflow->Fill(10., 1.);
        } else {
            continue;
        }

        // if (in->eRelPFIsoRho < 0.15) {
        //     cutflow->Fill(11., 1.);
        // } else {
        //     continue;
        // }

        // if ((isMC && !isEmbed) || in->bjetDeepCSVVeto20Medium_2018_DR0p5 < 1) {
        //     cutflow->Fill(12., 1.);
        // } else {
        //     continue;
        // }

        // implement new sorting per
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
        if (evt_now != evt_before) {  // new event, save the tau candidates
            //   since it is new event, do we have the best entry to save? If yes, save it!
            if (best_evt > -1) good_events.push_back(best_evt);

            //  this is a new event, so the first tau pair is the best! :)
            best_evt = ievt;
            eleCandidate = std::make_pair(ele.Pt(), in->eRelPFIsoRho);
            tauCandidate = std::make_pair(tau.Pt(), in->tDeepTau2017v2p1VSjetraw);
        } else {  // not a new event
            std::pair<float, float> currEleCandidate(ele.Pt(), in->eRelPFIsoRho);
            std::pair<float, float> currTauCandidate(tau.Pt(), in->tDeepTau2017v2p1VSjetraw);

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
                }      // tau1 has the same pT
            }          // tau1 has the same isolation
        }              // not a new event
        evt_before = evt_now;
    }
    if (best_evt > -1) good_events.push_back(best_evt);
}

void etau_tree2018::do_met_corr_nom(Float_t sf, TLorentzVector tau, TLorentzVector* met) {
    *met = (*met) + tau - sf * tau;  // update input met}
}

void etau_tree2018::do_recoil_corr(RecoilCorrector* recoilPFMetCorrector, TLorentzVector* met, int njets) {
    float pfmetcorr_ex, pfmetcorr_ey;
    recoilPFMetCorrector->CorrectByMeanResolution(met->Px(),      // uncorrected type I pf met px (float)
                                                  met->Py(),      // uncorrected type I pf met py (float)
                                                  in->genpX,      // generator Z/W/Higgs px (float)
                                                  in->genpY,      // generator Z/W/Higgs py (float)
                                                  in->vispX,      // generator visible Z/W/Higgs px (float)
                                                  in->vispY,      // generator visible Z/W/Higgs py (float)
                                                  njets,          // number of jets (hadronic jet multiplicity) (int)
                                                  pfmetcorr_ex,   // corrected type I pf met px (float, to be updated)
                                                  pfmetcorr_ey);  // corrected type I pf met py (float, to be updated)
    met->SetPxPyPzE(pfmetcorr_ex, pfmetcorr_ey, 0, sqrt(pfmetcorr_ex * pfmetcorr_ex + pfmetcorr_ey * pfmetcorr_ey));  // update MET
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
TTree* etau_tree2018::fill_tree(RecoilCorrector recoilPFMetCorrector, MEtSys metSys) {
    set_branches();  // get all the branches set up

    // loop through all events pasing skimming/sorting
    for (auto& ievt : good_events) {
        original->GetEntry(ievt);

        // remove anti-iso region from signal
        // if (isSignal && in->tMediumDeepTau2017v2p1VSjet < 0.5) {
        //     continue;
        // }

        Run = in->run;
        Lumi = in->lumi;

        // convert from Float_t in FSA to Int_t for analyzer
        gen_match_1 = in->eZTTGenMatching;
        gen_match_2 = in->tZTTGenMatching;
        njets = in->jetVeto30;
        nbtag = in->bjetDeepCSVVeto20Medium_2018_DR0p5;
        njetspt20 = in->jetVeto20;

        met_px = in->type1_pfMetEt * cos(in->type1_pfMetPhi);
        met_py = in->type1_pfMetEt * sin(in->type1_pfMetPhi);

        extraelec_veto = in->eVetoZTTp001dxyzR0 > 1;
        extramuon_veto = in->muVetoZTTp001dxyzR0 > 0;
        dilepton_veto = in->dielectronVeto > 0;

        // TLorentzVector ele, tau;
        ele.SetPtEtaPhiM(in->ePt, in->eEta, in->ePhi, in->eMass);
        tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);

        MET.SetPtEtaPhiM(in->type1_pfMetEt, 0, in->type1_pfMetPhi, 0);
        MET_reso_Up.SetPtEtaPhiM(in->type1_pfMetEt, 0, in->type1_pfMetPhi, 0);
        MET_reso_Down.SetPtEtaPhiM(in->type1_pfMetEt, 0, in->type1_pfMetPhi, 0);
        MET_resp_Up.SetPtEtaPhiM(in->type1_pfMetEt, 0, in->type1_pfMetPhi, 0);
        MET_resp_Down.SetPtEtaPhiM(in->type1_pfMetEt, 0, in->type1_pfMetPhi, 0);

        MET_JERUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JERUp, 0, in->type1_pfMet_shiftedPhi_JERUp, 0);
        MET_AbsoluteUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetAbsoluteUp, 0, in->type1_pfMet_shiftedPhi_JetAbsoluteUp, 0);
        MET_AbsoluteyearUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetAbsoluteyearUp, 0, in->type1_pfMet_shiftedPhi_JetAbsoluteyearUp, 0);
        MET_BBEC1Up.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetBBEC1Up, 0, in->type1_pfMet_shiftedPhi_JetBBEC1Up, 0);
        MET_BBEC1yearUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetBBEC1yearUp, 0, in->type1_pfMet_shiftedPhi_JetBBEC1yearUp, 0);
        MET_EC2Up.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEC2Up, 0, in->type1_pfMet_shiftedPhi_JetEC2Up, 0);
        MET_EC2yearUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEC2yearUp, 0, in->type1_pfMet_shiftedPhi_JetEC2yearUp, 0);
        MET_EnUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEnUp, 0, in->type1_pfMet_shiftedPhi_JetEnUp, 0);
        MET_FlavorQCDUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEnUp, 0, in->type1_pfMet_shiftedPhi_JetEnUp, 0);
        MET_HFUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetHFUp, 0, in->type1_pfMet_shiftedPhi_JetHFUp, 0);
        MET_HFyearUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetHFyearUp, 0, in->type1_pfMet_shiftedPhi_JetHFyearUp, 0);
        MET_RelBalUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetRelativeBalUp, 0, in->type1_pfMet_shiftedPhi_JetRelativeBalUp, 0);
        MET_RelSamUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetRelativeSampleUp, 0, in->type1_pfMet_shiftedPhi_JetRelativeSampleUp, 0);
        MET_ResUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetResUp, 0, in->type1_pfMet_shiftedPhi_JetResUp, 0);
        MET_TotalUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetTotalUp, 0, in->type1_pfMet_shiftedPhi_JetTotalUp, 0);
        MET_UESUp.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_UnclusteredEnUp, 0, in->type1_pfMet_shiftedPhi_UnclusteredEnUp, 0);

        MET_JERDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JERDown, 0, in->type1_pfMet_shiftedPhi_JERDown, 0);
        MET_AbsoluteDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetAbsoluteDown, 0, in->type1_pfMet_shiftedPhi_JetAbsoluteDown, 0);
        MET_AbsoluteyearDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetAbsoluteyearDown, 0, in->type1_pfMet_shiftedPhi_JetAbsoluteyearDown, 0);
        MET_BBEC1Down.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetBBEC1Down, 0, in->type1_pfMet_shiftedPhi_JetBBEC1Down, 0);
        MET_BBEC1yearDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetBBEC1yearDown, 0, in->type1_pfMet_shiftedPhi_JetBBEC1yearDown, 0);
        MET_EC2Down.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEC2Down, 0, in->type1_pfMet_shiftedPhi_JetEC2Down, 0);
        MET_EC2yearDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEC2yearDown, 0, in->type1_pfMet_shiftedPhi_JetEC2yearDown, 0);
        MET_EnDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEnDown, 0, in->type1_pfMet_shiftedPhi_JetEnDown, 0);
        MET_FlavorQCDDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetEnDown, 0, in->type1_pfMet_shiftedPhi_JetEnDown, 0);
        MET_HFDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetHFDown, 0, in->type1_pfMet_shiftedPhi_JetHFDown, 0);
        MET_HFyearDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetHFyearDown, 0, in->type1_pfMet_shiftedPhi_JetHFyearDown, 0);
        MET_RelBalDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetRelativeBalDown, 0, in->type1_pfMet_shiftedPhi_JetRelativeBalDown, 0);
        MET_RelSamDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetRelativeSampleDown, 0, in->type1_pfMet_shiftedPhi_JetRelativeSampleDown, 0);
        MET_ResDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetResDown, 0, in->type1_pfMet_shiftedPhi_JetResDown, 0);
        MET_TotalDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_JetTotalDown, 0, in->type1_pfMet_shiftedPhi_JetTotalDown, 0);
        MET_UESDown.SetPtEtaPhiM(in->type1_pfMet_shiftedPt_UnclusteredEnDown, 0, in->type1_pfMet_shiftedPhi_UnclusteredEnDown, 0);

        mets = {&MET,
            &MET_JERUp, &MET_JERDown,
            &MET_AbsoluteUp, &MET_AbsoluteDown,
            &MET_AbsoluteyearUp, &MET_AbsoluteyearDown,
            &MET_BBEC1Up, &MET_BBEC1Down,
            &MET_BBEC1yearUp, &MET_BBEC1yearDown,
            &MET_EC2Up, &MET_EC2Down,
            &MET_EC2yearUp, &MET_EC2yearDown,
            &MET_EnUp, &MET_EnDown,
            &MET_FlavorQCDUp, &MET_FlavorQCDDown,
            &MET_HFUp, &MET_HFDown,
            &MET_HFyearUp, &MET_HFyearDown,
            &MET_RelBalUp, &MET_RelBalDown,
            &MET_RelSamUp, &MET_RelSamDown,
            &MET_ResUp, &MET_ResDown,
            &MET_TotalUp, &MET_TotalDown,
            &MET_UESUp, &MET_UESDown};

        auto jet_for_correction = in->jetVeto30;
        if (recoil == 1) {
            jet_for_correction += 1;
        }

        // do recoil corrections on all met
        if (recoil > 0) {
            for (unsigned i = 0; i < mets.size(); i++) {
                do_recoil_corr(&recoilPFMetCorrector, mets.at(i), jet_for_correction);
            }

            float pfmetcorr_recoil_ex, pfmetcorr_recoil_ey;
            metSys.ApplyMEtSys(MET_resp_Up.Px(), MET_resp_Up.Py(), in->genpX, in->genpY, in->vispX, in->vispY, jet_for_correction,
                               MEtSys::ProcessType::BOSON, MEtSys::SysType::Response, MEtSys::SysShift::Up, pfmetcorr_recoil_ex, pfmetcorr_recoil_ey);
            MET_resp_Up.SetPxPyPzE(pfmetcorr_recoil_ex, pfmetcorr_recoil_ey, 0,
                                   sqrt(pfmetcorr_recoil_ex * pfmetcorr_recoil_ex + pfmetcorr_recoil_ey * pfmetcorr_recoil_ey));
            metSys.ApplyMEtSys(MET_resp_Down.Px(), MET_resp_Down.Py(), in->genpX, in->genpY, in->vispX, in->vispY, jet_for_correction,
                               MEtSys::ProcessType::BOSON, MEtSys::SysType::Response, MEtSys::SysShift::Down, pfmetcorr_recoil_ex, pfmetcorr_recoil_ey);
            MET_resp_Down.SetPxPyPzE(pfmetcorr_recoil_ex, pfmetcorr_recoil_ey, 0,
                                     sqrt(pfmetcorr_recoil_ex * pfmetcorr_recoil_ex + pfmetcorr_recoil_ey * pfmetcorr_recoil_ey));
            metSys.ApplyMEtSys(MET_reso_Up.Px(), MET_reso_Up.Py(), in->genpX, in->genpY, in->vispX, in->vispY, jet_for_correction,
                               MEtSys::ProcessType::BOSON, MEtSys::SysType::Resolution, MEtSys::SysShift::Up, pfmetcorr_recoil_ex, pfmetcorr_recoil_ey);
            MET_reso_Up.SetPxPyPzE(pfmetcorr_recoil_ex, pfmetcorr_recoil_ey, 0,
                                   sqrt(pfmetcorr_recoil_ex * pfmetcorr_recoil_ex + pfmetcorr_recoil_ey * pfmetcorr_recoil_ey));
            metSys.ApplyMEtSys(MET_reso_Down.Px(), MET_reso_Down.Py(), in->genpX, in->genpY, in->vispX, in->vispY, jet_for_correction,
                               MEtSys::ProcessType::BOSON, MEtSys::SysType::Resolution, MEtSys::SysShift::Down, pfmetcorr_recoil_ex, pfmetcorr_recoil_ey);
            MET_reso_Down.SetPxPyPzE(pfmetcorr_recoil_ex, pfmetcorr_recoil_ey, 0,
                                     sqrt(pfmetcorr_recoil_ex * pfmetcorr_recoil_ex + pfmetcorr_recoil_ey * pfmetcorr_recoil_ey));
        }

        tes_syst_up = 0;
        tes_syst_down = 0;
        ftes_syst_up = 0;
        ftes_syst_down = 0;
        if (isMC || isEmbed) {
            auto fes_sf = tfes.getFES(in->tDecayMode, in->tEta, in->tZTTGenMatching);
            auto tes_sf = tfes.getTES(in->tPt, in->tDecayMode, in->tZTTGenMatching);
            for (unsigned i = 0; i < mets.size(); i++) {
                do_met_corr_nom(fes_sf * tes_sf, tau, mets.at(i));  // correct for change in tau energy scale
                do_met_corr_nom(in->eCorrectedEt / ele.Energy(), ele, mets.at(i));  // correct for change in electron energy scale
            }
            tau = tau * fes_sf * tes_sf;
            ftes_syst_up = tfes.getFES(in->tDecayMode, in->tEta, in->tZTTGenMatching, "up");
            ftes_syst_down = tfes.getFES(in->tDecayMode, in->tEta, in->tZTTGenMatching, "down");
            tes_syst_up = tfes.getTES(in->tPt, in->tDecayMode, in->tZTTGenMatching, "up");
            tes_syst_down = tfes.getTES(in->tPt, in->tDecayMode, in->tZTTGenMatching, "down");
        }

        // electron energy scale
        if (isMC) {
            ele *= in->eCorrectedEt / ele.Energy();
        }

        met = MET.Pt();
        metphi = MET.Phi();
        met_px = MET.Px();
        met_py = MET.Py();

        // systematics
        met_resp_Up = MET_resp_Up.Pt();
        met_resp_Down = MET_resp_Down.Pt();
        met_reso_Up = MET_reso_Up.Pt();
        met_reso_Down = MET_reso_Down.Pt();
        metphi_resp_Up = MET_resp_Up.Phi();
        metphi_resp_Down = MET_resp_Down.Phi();
        metphi_reso_Up = MET_reso_Up.Phi();
        metphi_reso_Down = MET_reso_Down.Phi();

        met_JERUp = MET_JERUp.Pt();
        met_AbsoluteUp = MET_AbsoluteUp.Pt();
        met_AbsoluteyearUp = MET_AbsoluteyearUp.Pt();
        met_BBEC1Up = MET_BBEC1Up.Pt();
        met_BBEC1yearUp = MET_BBEC1yearUp.Pt();
        met_EC2Up = MET_EC2Up.Pt();
        met_EC2yearUp = MET_EC2yearUp.Pt();
        met_EnUp = MET_EnUp.Pt();
        met_FlavorQCDUp = MET_FlavorQCDUp.Pt();
        met_HFUp = MET_HFUp.Pt();
        met_HFyearUp = MET_HFyearUp.Pt();
        met_RelBalUp = MET_RelBalUp.Pt();
        met_RelSamUp = MET_RelSamUp.Pt();
        met_ResUp = MET_ResUp.Pt();
        met_TotalUp = MET_TotalUp.Pt();
        met_UESUp = MET_UESUp.Pt();
        met_JERDown = MET_JERDown.Pt();
        met_AbsoluteDown = MET_AbsoluteDown.Pt();
        met_AbsoluteyearDown = MET_AbsoluteyearDown.Pt();
        met_BBEC1Down = MET_BBEC1Down.Pt();
        met_BBEC1yearDown = MET_BBEC1yearDown.Pt();
        met_EC2Down = MET_EC2Down.Pt();
        met_EC2yearDown = MET_EC2yearDown.Pt();
        met_EnDown = MET_EnDown.Pt();
        met_FlavorQCDDown = MET_FlavorQCDDown.Pt();
        met_HFDown = MET_HFDown.Pt();
        met_HFyearDown = MET_HFyearDown.Pt();
        met_RelBalDown = MET_RelBalDown.Pt();
        met_RelSamDown = MET_RelSamDown.Pt();
        met_ResDown = MET_ResDown.Pt();
        met_TotalDown = MET_TotalDown.Pt();
        met_UESDown = MET_UESDown.Pt();

        metphi_JERUp = MET_JERUp.Phi();
        metphi_AbsoluteUp = MET_AbsoluteUp.Phi();
        metphi_AbsoluteyearUp = MET_AbsoluteyearUp.Phi();
        metphi_BBEC1Up = MET_BBEC1Up.Phi();
        metphi_BBEC1yearUp = MET_BBEC1yearUp.Phi();
        metphi_EC2Up = MET_EC2Up.Phi();
        metphi_EC2yearUp = MET_EC2yearUp.Phi();
        metphi_EnUp = MET_EnUp.Phi();
        metphi_FlavorQCDUp = MET_FlavorQCDUp.Phi();
        metphi_HFUp = MET_HFUp.Phi();
        metphi_HFyearUp = MET_HFyearUp.Phi();
        metphi_RelBalUp = MET_RelBalUp.Phi();
        metphi_RelSamUp = MET_RelSamUp.Phi();
        metphi_ResUp = MET_ResUp.Phi();
        metphi_TotalUp = MET_TotalUp.Phi();
        metphi_UESUp = MET_UESUp.Phi();
        metphi_JERDown = MET_JERDown.Phi();
        metphi_AbsoluteDown = MET_AbsoluteDown.Phi();
        metphi_AbsoluteyearDown = MET_AbsoluteyearDown.Phi();
        metphi_BBEC1Down = MET_BBEC1Down.Phi();
        metphi_BBEC1yearDown = MET_BBEC1yearDown.Phi();
        metphi_EC2Down = MET_EC2Down.Phi();
        metphi_EC2yearDown = MET_EC2yearDown.Phi();
        metphi_EnDown = MET_EnDown.Phi();
        metphi_FlavorQCDDown = MET_FlavorQCDDown.Phi();
        metphi_HFDown = MET_HFDown.Phi();
        metphi_HFyearDown = MET_HFyearDown.Phi();
        metphi_RelBalDown = MET_RelBalDown.Phi();
        metphi_RelSamDown = MET_RelSamDown.Phi();
        metphi_ResDown = MET_ResDown.Phi();
        metphi_TotalDown = MET_TotalDown.Phi();
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
    // new branches
    tree->Branch("era", &era);
    tree->Branch("run", &Run);
    tree->Branch("lumi", &Lumi);
    tree->Branch("gen_match_1", &gen_match_1);
    tree->Branch("gen_match_2", &gen_match_2);
    tree->Branch("extraelec_veto", &extraelec_veto);
    tree->Branch("extramuon_veto", &extramuon_veto);
    tree->Branch("dilepton_veto", &dilepton_veto);
    tree->Branch("met", &met);
    tree->Branch("metphi", &metphi);
    tree->Branch("met_px", &met_px);
    tree->Branch("met_py", &met_py);
    tree->Branch("met_reso_Up", &met_reso_Up);
    tree->Branch("met_reso_Down", &met_reso_Down);
    tree->Branch("met_resp_Up", &met_resp_Up);
    tree->Branch("met_resp_Down", &met_resp_Down);
    tree->Branch("metphi_reso_Up", &metphi_reso_Up);
    tree->Branch("metphi_reso_Down", &metphi_reso_Down);
    tree->Branch("metphi_resp_Up", &metphi_resp_Up);
    tree->Branch("metphi_resp_Down", &metphi_resp_Down);
    tree->Branch("met_JERUp", &met_JERUp);
    tree->Branch("met_AbsoluteUp", &met_AbsoluteUp);
    tree->Branch("met_AbsoluteyearUp", &met_AbsoluteyearUp);
    tree->Branch("met_BBEC1Up", &met_BBEC1Up);
    tree->Branch("met_BBEC1yearUp", &met_BBEC1yearUp);
    tree->Branch("met_EC2Up", &met_EC2Up);
    tree->Branch("met_EC2yearUp", &met_EC2yearUp);
    tree->Branch("met_EnUp", &met_EnUp);
    tree->Branch("met_FlavorQCDUp", &met_FlavorQCDUp);
    tree->Branch("met_HFUp", &met_HFUp);
    tree->Branch("met_HFyearUp", &met_HFyearUp);
    tree->Branch("met_RelBalUp", &met_RelBalUp);
    tree->Branch("met_RelSamUp", &met_RelSamUp);
    tree->Branch("met_ResUp", &met_ResUp);
    tree->Branch("met_TotalUp", &met_TotalUp);
    tree->Branch("met_UESUp", &met_UESUp);
    tree->Branch("met_JERDown", &met_JERDown);
    tree->Branch("met_AbsoluteDown", &met_AbsoluteDown);
    tree->Branch("met_AbsoluteyearDown", &met_AbsoluteyearDown);
    tree->Branch("met_BBEC1Down", &met_BBEC1Down);
    tree->Branch("met_BBEC1yearDown", &met_BBEC1yearDown);
    tree->Branch("met_EC2Down", &met_EC2Down);
    tree->Branch("met_EC2yearDown", &met_EC2yearDown);
    tree->Branch("met_EnDown", &met_EnDown);
    tree->Branch("met_FlavorQCDDown", &met_FlavorQCDDown);
    tree->Branch("met_HFDown", &met_HFDown);
    tree->Branch("met_HFyearDown", &met_HFyearDown);
    tree->Branch("met_RelBalDown", &met_RelBalDown);
    tree->Branch("met_RelSamDown", &met_RelSamDown);
    tree->Branch("met_ResDown", &met_ResDown);
    tree->Branch("met_TotalDown", &met_TotalDown);
    tree->Branch("met_UESDown", &met_UESDown);
    tree->Branch("metphi_JERUp", &metphi_JERUp);
    tree->Branch("metphi_AbsoluteUp", &metphi_AbsoluteUp);
    tree->Branch("metphi_AbsoluteyearUp", &metphi_AbsoluteyearUp);
    tree->Branch("metphi_BBEC1Up", &metphi_BBEC1Up);
    tree->Branch("metphi_BBEC1yearUp", &metphi_BBEC1yearUp);
    tree->Branch("metphi_EC2Up", &metphi_EC2Up);
    tree->Branch("metphi_EC2yearUp", &metphi_EC2yearUp);
    tree->Branch("metphi_EnUp", &metphi_EnUp);
    tree->Branch("metphi_FlavorQCDUp", &metphi_FlavorQCDUp);
    tree->Branch("metphi_HFUp", &metphi_HFUp);
    tree->Branch("metphi_HFyearUp", &metphi_HFyearUp);
    tree->Branch("metphi_RelBalUp", &metphi_RelBalUp);
    tree->Branch("metphi_RelSamUp", &metphi_RelSamUp);
    tree->Branch("metphi_ResUp", &metphi_ResUp);
    tree->Branch("metphi_TotalUp", &metphi_TotalUp);
    tree->Branch("metphi_UESUp", &metphi_UESUp);
    tree->Branch("metphi_JERDown", &metphi_JERDown);
    tree->Branch("metphi_AbsoluteDown", &metphi_AbsoluteDown);
    tree->Branch("metphi_AbsoluteyearDown", &metphi_AbsoluteyearDown);
    tree->Branch("metphi_BBEC1Down", &metphi_BBEC1Down);
    tree->Branch("metphi_BBEC1yearDown", &metphi_BBEC1yearDown);
    tree->Branch("metphi_EC2Down", &metphi_EC2Down);
    tree->Branch("metphi_EC2yearDown", &metphi_EC2yearDown);
    tree->Branch("metphi_EnDown", &metphi_EnDown);
    tree->Branch("metphi_FlavorQCDDown", &metphi_FlavorQCDDown);
    tree->Branch("metphi_HFDown", &metphi_HFDown);
    tree->Branch("metphi_HFyearDown", &metphi_HFyearDown);
    tree->Branch("metphi_RelBalDown", &metphi_RelBalDown);
    tree->Branch("metphi_RelSamDown", &metphi_RelSamDown);
    tree->Branch("metphi_ResDown", &metphi_ResDown);
    tree->Branch("metphi_TotalDown", &metphi_TotalDown);
    tree->Branch("metphi_UESDown", &metphi_UESDown);
    tree->Branch("ftes_syst_up", &ftes_syst_up);
    tree->Branch("ftes_syst_down", &ftes_syst_down);
    tree->Branch("tes_syst_up", &tes_syst_up);
    tree->Branch("tes_syst_down", &tes_syst_down);

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
    tree->Branch("q_1", &in->eCharge);
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
    tree->Branch("Ele24LooseHPSTau30Pass", &in->Ele24LooseHPSTau30Pass);
    tree->Branch("Ele24LooseHPSTau30TightIDPass", &in->Ele24LooseHPSTau30TightIDPass);
    tree->Branch("Ele24LooseTau30Pass", &in->Ele24LooseTau30Pass);
    tree->Branch("Ele24LooseTau30TightIDPass", &in->Ele24LooseTau30TightIDPass);
    // tree->Branch("Ele27WPTightPass", &in->Ele27WPTightPass);
    // tree->Branch("Ele32WPTightPass", &in->Ele32WPTightPass);
    // tree->Branch("Ele35WPTightPass", &in->Ele35WPTightPass);
    // tree->Branch("Ele38WPTightPass", &in->Ele38WPTightPass);
    // tree->Branch("Ele40WPTightPass", &in->Ele40WPTightPass);
    tree->Branch("EmbPtWeight", &in->EmbPtWeight);
    // tree->Branch("Eta", &in->Eta);
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
    // tree->Branch("Ht", &in->Ht);
    // tree->Branch("IsoMu24Pass", &in->IsoMu24Pass);
    // tree->Branch("IsoMu27Pass", &in->IsoMu27Pass);
    // tree->Branch("LT", &in->LT);
    // tree->Branch("Mass", &in->Mass);
    // tree->Branch("MassError", &in->MassError);
    // tree->Branch("MassErrord1", &in->MassErrord1);
    // tree->Branch("MassErrord2", &in->MassErrord2);
    // tree->Branch("MassErrord3", &in->MassErrord3);
    // tree->Branch("MassErrord4", &in->MassErrord4);
    // tree->Branch("Mt", &in->Mt);
    // tree->Branch("Mu20LooseHPSTau27Pass", &in->Mu20LooseHPSTau27Pass);
    // tree->Branch("Mu20LooseHPSTau27TightIDPass", &in->Mu20LooseHPSTau27TightIDPass);
    // tree->Branch("Mu20LooseTau27Pass", &in->Mu20LooseTau27Pass);
    // tree->Branch("Mu20LooseTau27TightIDPass", &in->Mu20LooseTau27TightIDPass);
    // tree->Branch("Mu50Pass", &in->Mu50Pass);
    tree->Branch("NUP", &in->NUP);
//    tree->Branch("Phi", &in->Phi);
    // tree->Branch("Pt", &in->Pt);
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
    // tree->Branch("bjetDeepCSVVeto20Loose_2016_DR0p5", &in->bjetDeepCSVVeto20Loose_2016_DR0p5);
    // tree->Branch("bjetDeepCSVVeto20Loose_2017_DR0p5", &in->bjetDeepCSVVeto20Loose_2017_DR0p5);
    tree->Branch("bjetDeepCSVVeto20Loose_2018_DR0p5", &in->bjetDeepCSVVeto20Loose_2018_DR0p5);
    // tree->Branch("bjetDeepCSVVeto20Medium_2016_DR0", &in->bjetDeepCSVVeto20Medium_2016_DR0);
    // tree->Branch("bjetDeepCSVVeto20Medium_2016_DR0p5", &in->bjetDeepCSVVeto20Medium_2016_DR0p5);
    // tree->Branch("bjetDeepCSVVeto20Medium_2017_DR0", &in->bjetDeepCSVVeto20Medium_2017_DR0);
    // tree->Branch("bjetDeepCSVVeto20Medium_2017_DR0p5", &in->bjetDeepCSVVeto20Medium_2017_DR0p5);
    tree->Branch("bjetDeepCSVVeto20Medium_2018_DR0", &in->bjetDeepCSVVeto20Medium_2018_DR0);
    tree->Branch("bjetDeepCSVVeto20Medium_2018_DR0p5", &in->bjetDeepCSVVeto20Medium_2018_DR0p5);
    // tree->Branch("bjetDeepCSVVeto20Tight_2016_DR0p5", &in->bjetDeepCSVVeto20Tight_2016_DR0p5);
    // tree->Branch("bjetDeepCSVVeto20Tight_2017_DR0p5", &in->bjetDeepCSVVeto20Tight_2017_DR0p5);
    tree->Branch("bjetDeepCSVVeto20Tight_2018_DR0p5", &in->bjetDeepCSVVeto20Tight_2018_DR0p5);
    // tree->Branch("bweight_2016", &in->bweight_2016);
    // tree->Branch("bweight_2017", &in->bweight_2017);
    // tree->Branch("bweight_2018", &in->bweight_2018);
    tree->Branch("charge", &in->charge);
    tree->Branch("dielectronVeto", &in->dielectronVeto);
    // tree->Branch("dimu9ele9Pass", &in->dimu9ele9Pass);
    tree->Branch("dimuonVeto", &in->dimuonVeto);
    // tree->Branch("eCBIDLoose", &in->eCBIDLoose);
    // tree->Branch("eCBIDMedium", &in->eCBIDMedium);
    // tree->Branch("eCBIDTight", &in->eCBIDTight);
    // tree->Branch("eCBIDVeto", &in->eCBIDVeto);
    tree->Branch("eCharge", &in->eCharge);
    // tree->Branch("eChargeIdLoose", &in->eChargeIdLoose);
    // tree->Branch("eChargeIdMed", &in->eChargeIdMed);
    // tree->Branch("eChargeIdTight", &in->eChargeIdTight);
    // tree->Branch("eComesFromHiggs", &in->eComesFromHiggs);
    tree->Branch("eCorrectedEt", &in->eCorrectedEt);
    // tree->Branch("eE1x5", &in->eE1x5);
    // tree->Branch("eE2x5Max", &in->eE2x5Max);
    // tree->Branch("eE5x5", &in->eE5x5);
    tree->Branch("eEcalIsoDR03", &in->eEcalIsoDR03);
    tree->Branch("eEnergyError", &in->eEnergyError);
    tree->Branch("eEnergyScaleDown", &in->eEnergyScaleDown);
    tree->Branch("eEnergyScaleGainDown", &in->eEnergyScaleGainDown);
    tree->Branch("eEnergyScaleGainUp", &in->eEnergyScaleGainUp);
    tree->Branch("eEnergyScaleStatDown", &in->eEnergyScaleStatDown);
    tree->Branch("eEnergyScaleStatUp", &in->eEnergyScaleStatUp);
    tree->Branch("eEnergyScaleSystDown", &in->eEnergyScaleSystDown);
    tree->Branch("eEnergyScaleSystUp", &in->eEnergyScaleSystUp);
    tree->Branch("eEnergyScaleUp", &in->eEnergyScaleUp);
    tree->Branch("eEnergySigmaDown", &in->eEnergySigmaDown);
    tree->Branch("eEnergySigmaPhiDown", &in->eEnergySigmaPhiDown);
    tree->Branch("eEnergySigmaPhiUp", &in->eEnergySigmaPhiUp);
    tree->Branch("eEnergySigmaRhoDown", &in->eEnergySigmaRhoDown);
    tree->Branch("eEnergySigmaRhoUp", &in->eEnergySigmaRhoUp);
    tree->Branch("eEnergySigmaUp", &in->eEnergySigmaUp);
    tree->Branch("eEta", &in->eEta);
    tree->Branch("eGenCharge", &in->eGenCharge);
    // tree->Branch("eGenDirectPromptTauDecay", &in->eGenDirectPromptTauDecay);
    tree->Branch("eGenEnergy", &in->eGenEnergy);
    tree->Branch("eGenEta", &in->eGenEta);
    tree->Branch("eGenIsPrompt", &in->eGenIsPrompt);
    tree->Branch("eGenMotherPdgId", &in->eGenMotherPdgId);
    tree->Branch("eGenParticle", &in->eGenParticle);
    tree->Branch("eGenPdgId", &in->eGenPdgId);
    tree->Branch("eGenPhi", &in->eGenPhi);
    tree->Branch("eGenPrompt", &in->eGenPrompt);
    tree->Branch("eGenPromptTauDecay", &in->eGenPromptTauDecay);
    tree->Branch("eGenPt", &in->eGenPt);
    tree->Branch("eGenTauDecay", &in->eGenTauDecay);
    tree->Branch("eGenVZ", &in->eGenVZ);
    tree->Branch("eGenVtxPVMatch", &in->eGenVtxPVMatch);
    // tree->Branch("eHadronicDepth1OverEm", &in->eHadronicDepth1OverEm);
    // tree->Branch("eHadronicDepth2OverEm", &in->eHadronicDepth2OverEm);
    // tree->Branch("eHadronicOverEM", &in->eHadronicOverEM);
    // tree->Branch("eHcalIsoDR03", &in->eHcalIsoDR03);
    // tree->Branch("eIP3D", &in->eIP3D);
    // tree->Branch("eIP3DErr", &in->eIP3DErr);
    // tree->Branch("eIsoDB03", &in->eIsoDB03);
    // tree->Branch("eJetArea", &in->eJetArea);
    // tree->Branch("eJetBtag", &in->eJetBtag);
    // tree->Branch("eJetDR", &in->eJetDR);
    // tree->Branch("eJetEtaEtaMoment", &in->eJetEtaEtaMoment);
    // tree->Branch("eJetEtaPhiMoment", &in->eJetEtaPhiMoment);
    // tree->Branch("eJetEtaPhiSpread", &in->eJetEtaPhiSpread);
    // tree->Branch("eJetHadronFlavour", &in->eJetHadronFlavour);
    // tree->Branch("eJetPFCISVBtag", &in->eJetPFCISVBtag);
    // tree->Branch("eJetPartonFlavour", &in->eJetPartonFlavour);
    // tree->Branch("eJetPhiPhiMoment", &in->eJetPhiPhiMoment);
    // tree->Branch("eJetPt", &in->eJetPt);
    // tree->Branch("eLowestMll", &in->eLowestMll);
    tree->Branch("eMVAIsoWP80", &in->eMVAIsoWP80);
    tree->Branch("eMVAIsoWP90", &in->eMVAIsoWP90);
    // tree->Branch("eMVAIsoWPHZZ", &in->eMVAIsoWPHZZ);
    // tree->Branch("eMVAIsoWPLoose", &in->eMVAIsoWPLoose);
    tree->Branch("eMVANoisoWP80", &in->eMVANoisoWP80);
    tree->Branch("eMVANoisoWP90", &in->eMVANoisoWP90);
    tree->Branch("eMVANoisoWPLoose", &in->eMVANoisoWPLoose);
    tree->Branch("eMass", &in->eMass);
    tree->Branch("eMatchesEle24HPSTau30Filter", &in->eMatchesEle24HPSTau30Filter);
    tree->Branch("eMatchesEle24HPSTau30Path", &in->eMatchesEle24HPSTau30Path);
    tree->Branch("eMatchesEle24Tau30Filter", &in->eMatchesEle24Tau30Filter);
    tree->Branch("eMatchesEle24Tau30Path", &in->eMatchesEle24Tau30Path);
    tree->Branch("eMatchesEle25Filter", &in->eMatchesEle25Filter);
    tree->Branch("eMatchesEle25Path", &in->eMatchesEle25Path);
    tree->Branch("eMatchesEle27Filter", &in->eMatchesEle27Filter);
    tree->Branch("eMatchesEle27Path", &in->eMatchesEle27Path);
    tree->Branch("eMatchesEle32Filter", &in->eMatchesEle32Filter);
    tree->Branch("eMatchesEle32Path", &in->eMatchesEle32Path);
    tree->Branch("eMatchesEle35Filter", &in->eMatchesEle35Filter);
    tree->Branch("eMatchesEle35Path", &in->eMatchesEle35Path);
    tree->Branch("eMatchEmbeddedFilterEle32", &in->eMatchEmbeddedFilterEle32);
    tree->Branch("eMatchEmbeddedFilterEle35", &in->eMatchEmbeddedFilterEle35);
    tree->Branch("eMatchEmbeddedFilterEle24Tau30", &in->eMatchEmbeddedFilterEle24Tau30);
    // tree->Branch("eMissingHits", &in->eMissingHits);
    // tree->Branch("eNearMuonVeto", &in->eNearMuonVeto);
    // tree->Branch("eNearestMuonDR", &in->eNearestMuonDR);
    // tree->Branch("eNearestZMass", &in->eNearestZMass);
    // tree->Branch("ePFChargedIso", &in->ePFChargedIso);
    // tree->Branch("ePFNeutralIso", &in->ePFNeutralIso);
    // tree->Branch("ePFPUChargedIso", &in->ePFPUChargedIso);
    // tree->Branch("ePFPhotonIso", &in->ePFPhotonIso);
    tree->Branch("ePVDXY", &in->ePVDXY);
    tree->Branch("ePVDZ", &in->ePVDZ);
    // tree->Branch("ePassesConversionVeto", &in->ePassesConversionVeto);
    tree->Branch("ePhi", &in->ePhi);
    tree->Branch("ePt", &in->ePt);
    tree->Branch("eRelIso", &in->eRelIso);
    tree->Branch("eRelPFIsoDB", &in->eRelPFIsoDB);
    tree->Branch("eRelPFIsoRho", &in->eRelPFIsoRho);
    // tree->Branch("eRho", &in->eRho);
    // tree->Branch("eSCEnergy", &in->eSCEnergy);
    // tree->Branch("eSCEta", &in->eSCEta);
    // tree->Branch("eSCEtaWidth", &in->eSCEtaWidth);
    // tree->Branch("eSCPhi", &in->eSCPhi);
    // tree->Branch("eSCPhiWidth", &in->eSCPhiWidth);
    // tree->Branch("eSCPreshowerEnergy", &in->eSCPreshowerEnergy);
    // tree->Branch("eSCRawEnergy", &in->eSCRawEnergy);
    // tree->Branch("eSIP2D", &in->eSIP2D);
    // tree->Branch("eSIP3D", &in->eSIP3D);
    // tree->Branch("eSigmaIEtaIEta", &in->eSigmaIEtaIEta);
    // tree->Branch("eTrkIsoDR03", &in->eTrkIsoDR03);
    // tree->Branch("eVZ", &in->eVZ);
    // tree->Branch("eVetoHZZPt5", &in->eVetoHZZPt5);
    // tree->Branch("eVetoZTTp001dxyz", &in->eVetoZTTp001dxyz);
    // tree->Branch("eVetoZTTp001dxyzR0", &in->eVetoZTTp001dxyzR0);
    tree->Branch("eZTTGenMatching", &in->eZTTGenMatching);
    // tree->Branch("e_t_DR", &in->e_t_DR);
    // tree->Branch("e_t_Mass", &in->e_t_Mass);
    // tree->Branch("e_t_doubleL1IsoTauMatch", &in->e_t_doubleL1IsoTauMatch);
    // tree->Branch("edeltaEtaSuperClusterTrackAtVtx", &in->edeltaEtaSuperClusterTrackAtVtx);
    // tree->Branch("edeltaPhiSuperClusterTrackAtVtx", &in->edeltaPhiSuperClusterTrackAtVtx);
    // tree->Branch("eeSuperClusterOverP", &in->eeSuperClusterOverP);
    // tree->Branch("eecalEnergy", &in->eecalEnergy);
    // tree->Branch("efBrem", &in->efBrem);
    // tree->Branch("etrackMomentumAtVtxP", &in->etrackMomentumAtVtxP);
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
    tree->Branch("j1eta", &in->j1eta);
    tree->Branch("j1hadronflavor", &in->j1hadronflavor);
    tree->Branch("j1phi", &in->j1phi);
    tree->Branch("j1pt", &in->j1pt);
    tree->Branch("j2csv", &in->j2csv);
    tree->Branch("j2eta", &in->j2eta);
    tree->Branch("j2hadronflavor", &in->j2hadronflavor);
    tree->Branch("j2phi", &in->j2phi);
    tree->Branch("j2pt", &in->j2pt);
    tree->Branch("jb1eta", &in->jb1eta);
    tree->Branch("deepcsvb1_btagscore", &in->deepcsvb1_btagscore);
    tree->Branch("jb1hadronflavor", &in->jb1hadronflavor);
    tree->Branch("jb1phi", &in->jb1phi);
    tree->Branch("jb1pt", &in->jb1pt);
    tree->Branch("jb2eta", &in->jb2eta);
    tree->Branch("deepcsvb2_btagscore", &in->deepcsvb2_btagscore);
    tree->Branch("jb2hadronflavor", &in->jb2hadronflavor);
    tree->Branch("jb2phi", &in->jb2phi);
    tree->Branch("jb2pt", &in->jb2pt);
    tree->Branch("jetVeto20", &in->jetVeto20);
    tree->Branch("jetVeto20_JetEnDown", &in->jetVeto20_JetEnDown);
    tree->Branch("jetVeto20_JetEnUp", &in->jetVeto20_JetEnUp);
    tree->Branch("jetVeto30", &in->jetVeto30);
    tree->Branch("jetVeto30_JERDown", &in->jetVeto30_JERDown);
    tree->Branch("jetVeto30_JERUp", &in->jetVeto30_JERUp);
    tree->Branch("jetVeto30_JetAbsoluteDown", &in->jetVeto30_JetAbsoluteDown);
    tree->Branch("jetVeto30_JetAbsoluteUp", &in->jetVeto30_JetAbsoluteUp);
    tree->Branch("jetVeto30_JetAbsoluteyearDown", &in->jetVeto30_JetAbsoluteyearDown);
    tree->Branch("jetVeto30_JetAbsoluteyearUp", &in->jetVeto30_JetAbsoluteyearUp);
    tree->Branch("jetVeto30_JetBBEC1Down", &in->jetVeto30_JetBBEC1Down);
    tree->Branch("jetVeto30_JetBBEC1Up", &in->jetVeto30_JetBBEC1Up);
    tree->Branch("jetVeto30_JetBBEC1yearDown", &in->jetVeto30_JetBBEC1yearDown);
    tree->Branch("jetVeto30_JetBBEC1yearUp", &in->jetVeto30_JetBBEC1yearUp);
    tree->Branch("jetVeto30_JetEC2Down", &in->jetVeto30_JetEC2Down);
    tree->Branch("jetVeto30_JetEC2Up", &in->jetVeto30_JetEC2Up);
    tree->Branch("jetVeto30_JetEC2yearDown", &in->jetVeto30_JetEC2yearDown);
    tree->Branch("jetVeto30_JetEC2yearUp", &in->jetVeto30_JetEC2yearUp);
    tree->Branch("jetVeto30_JetEnDown", &in->jetVeto30_JetEnDown);
    tree->Branch("jetVeto30_JetEnUp", &in->jetVeto30_JetEnUp);
    tree->Branch("jetVeto30_JetFlavorQCDDown", &in->jetVeto30_JetFlavorQCDDown);
    tree->Branch("jetVeto30_JetFlavorQCDUp", &in->jetVeto30_JetFlavorQCDUp);
    tree->Branch("jetVeto30_JetHFDown", &in->jetVeto30_JetHFDown);
    tree->Branch("jetVeto30_JetHFUp", &in->jetVeto30_JetHFUp);
    tree->Branch("jetVeto30_JetHFyearDown", &in->jetVeto30_JetHFyearDown);
    tree->Branch("jetVeto30_JetHFyearUp", &in->jetVeto30_JetHFyearUp);
    tree->Branch("jetVeto30_JetRelativeBalDown", &in->jetVeto30_JetRelativeBalDown);
    tree->Branch("jetVeto30_JetRelativeBalUp", &in->jetVeto30_JetRelativeBalUp);
    tree->Branch("jetVeto30_JetRelativeSampleDown", &in->jetVeto30_JetRelativeSampleDown);
    tree->Branch("jetVeto30_JetRelativeSampleUp", &in->jetVeto30_JetRelativeSampleUp);
    tree->Branch("jetVeto30_JetTotalDown", &in->jetVeto30_JetTotalDown);
    tree->Branch("jetVeto30_JetTotalUp", &in->jetVeto30_JetTotalUp);
    tree->Branch("lumi", &in->lumi);
    tree->Branch("metSig", &in->metSig);
    tree->Branch("metcov00", &in->metcov00);
    tree->Branch("metcov01", &in->metcov01);
    tree->Branch("metcov10", &in->metcov10);
    tree->Branch("metcov11", &in->metcov11);
    // tree->Branch("mu12e23DZPass", &in->mu12e23DZPass);
    // tree->Branch("mu12e23Pass", &in->mu12e23Pass);
    // tree->Branch("mu23e12DZPass", &in->mu23e12DZPass);
    // tree->Branch("mu23e12Pass", &in->mu23e12Pass);
    // tree->Branch("mu8diele12DZPass", &in->mu8diele12DZPass);
    // tree->Branch("mu8diele12Pass", &in->mu8diele12Pass);
    // tree->Branch("mu8e23DZPass", &in->mu8e23DZPass);
    // tree->Branch("mu8e23Pass", &in->mu8e23Pass);
    // tree->Branch("muGlbIsoVetoPt10", &in->muGlbIsoVetoPt10);
    // tree->Branch("muVeto5", &in->muVeto5);
    // tree->Branch("muVetoZTTp001dxyz", &in->muVetoZTTp001dxyz);
    // tree->Branch("muVetoZTTp001dxyzR0", &in->muVetoZTTp001dxyzR0);
    tree->Branch("nTruePU", &in->nTruePU);
    tree->Branch("npNLO", &in->npNLO);
    tree->Branch("numGenJets", &in->numGenJets);
    tree->Branch("nvtx", &in->nvtx);
    tree->Branch("prefiring_weight", &in->prefiring_weight);
    tree->Branch("prefiring_weight_up", &in->prefiring_weight_up);
    tree->Branch("prefiring_weight_down", &in->prefiring_weight_down);
    tree->Branch("processID", &in->processID);
    // tree->Branch("puppiMetEt", &in->puppiMetEt);
    // tree->Branch("puppiMetPhi", &in->puppiMetPhi);
    // tree->Branch("pvChi2", &in->pvChi2);
    // tree->Branch("pvDX", &in->pvDX);
    // tree->Branch("pvDY", &in->pvDY);
    // tree->Branch("pvDZ", &in->pvDZ);
    // tree->Branch("pvIsFake", &in->pvIsFake);
    // tree->Branch("pvIsValid", &in->pvIsValid);
    // tree->Branch("pvNormChi2", &in->pvNormChi2);
    // tree->Branch("pvRho", &in->pvRho);
    // tree->Branch("pvX", &in->pvX);
    // tree->Branch("pvY", &in->pvY);
    // tree->Branch("pvZ", &in->pvZ);
    // tree->Branch("pvndof", &in->pvndof);
    tree->Branch("raw_pfMetEt", &in->raw_pfMetEt);
    tree->Branch("raw_pfMetPhi", &in->raw_pfMetPhi);
    // tree->Branch("recoilDaught", &in->recoilDaught);
    // tree->Branch("recoilWithMet", &in->recoilWithMet);
    tree->Branch("rho", &in->rho);
    tree->Branch("run", &in->run);
    tree->Branch("singleE25eta2p1TightPass", &in->singleE25eta2p1TightPass);
    tree->Branch("singleIsoMu22Pass", &in->singleIsoMu22Pass);
    tree->Branch("singleIsoMu22eta2p1Pass", &in->singleIsoMu22eta2p1Pass);
    tree->Branch("singleIsoTkMu22Pass", &in->singleIsoTkMu22Pass);
    tree->Branch("singleIsoTkMu22eta2p1Pass", &in->singleIsoTkMu22eta2p1Pass);
    tree->Branch("singleMu19eta2p1LooseTau20Pass", &in->singleMu19eta2p1LooseTau20Pass);
    tree->Branch("singleMu19eta2p1LooseTau20singleL1Pass", &in->singleMu19eta2p1LooseTau20singleL1Pass);
    // tree->Branch("tAgainstElectronLooseMVA6", &in->tAgainstElectronLooseMVA6);
    // tree->Branch("tAgainstElectronLooseMVA62018", &in->tAgainstElectronLooseMVA62018);
    // tree->Branch("tAgainstElectronMVA6Raw", &in->tAgainstElectronMVA6Raw);
    // tree->Branch("tAgainstElectronMVA6Raw2018", &in->tAgainstElectronMVA6Raw2018);
    // tree->Branch("tAgainstElectronMVA6category", &in->tAgainstElectronMVA6category);
    // tree->Branch("tAgainstElectronMVA6category2018", &in->tAgainstElectronMVA6category2018);
    // tree->Branch("tAgainstElectronMediumMVA6", &in->tAgainstElectronMediumMVA6);
    // tree->Branch("tAgainstElectronMediumMVA62018", &in->tAgainstElectronMediumMVA62018);
    // tree->Branch("tAgainstElectronTightMVA6", &in->tAgainstElectronTightMVA6);
    // tree->Branch("tAgainstElectronTightMVA62018", &in->tAgainstElectronTightMVA62018);
    // tree->Branch("tAgainstElectronVLooseMVA6", &in->tAgainstElectronVLooseMVA6);
    // tree->Branch("tAgainstElectronVLooseMVA62018", &in->tAgainstElectronVLooseMVA62018);
    // tree->Branch("tAgainstElectronVTightMVA6", &in->tAgainstElectronVTightMVA6);
    // tree->Branch("tAgainstElectronVTightMVA62018", &in->tAgainstElectronVTightMVA62018);
    // tree->Branch("tAgainstMuonLoose3", &in->tAgainstMuonLoose3);
    // tree->Branch("tAgainstMuonTight3", &in->tAgainstMuonTight3);
    tree->Branch("tCharge", &in->tCharge);
    // tree->Branch("tChargedIsoPtSum", &in->tChargedIsoPtSum);
    // tree->Branch("tChargedIsoPtSumdR03", &in->tChargedIsoPtSumdR03);
    // tree->Branch("tComesFromHiggs", &in->tComesFromHiggs);
    tree->Branch("tDecayMode", &in->tDecayMode);
    tree->Branch("tDecayModeFinding", &in->tDecayModeFinding);
    tree->Branch("tDecayModeFindingNewDMs", &in->tDecayModeFindingNewDMs);
    tree->Branch("tDeepTau2017v2p1VSeraw", &in->tDeepTau2017v2p1VSeraw);
    tree->Branch("tDeepTau2017v2p1VSjetraw", &in->tDeepTau2017v2p1VSjetraw);
    tree->Branch("tDeepTau2017v2p1VSmuraw", &in->tDeepTau2017v2p1VSmuraw);
    tree->Branch("tEta", &in->tEta);
    // tree->Branch("tFootprintCorrection", &in->tFootprintCorrection);
    // tree->Branch("tFootprintCorrectiondR03", &in->tFootprintCorrectiondR03);
    // tree->Branch("tGenCharge", &in->tGenCharge);
    // tree->Branch("tGenDecayMode", &in->tGenDecayMode);
    // tree->Branch("tGenEnergy", &in->tGenEnergy);
    tree->Branch("tGenEta", &in->tGenEta);
    // tree->Branch("tGenJetEta", &in->tGenJetEta);
    // tree->Branch("tGenJetPt", &in->tGenJetPt);
    // tree->Branch("tGenMotherEnergy", &in->tGenMotherEnergy);
    // tree->Branch("tGenMotherEta", &in->tGenMotherEta);
    // tree->Branch("tGenMotherPdgId", &in->tGenMotherPdgId);
    // tree->Branch("tGenMotherPhi", &in->tGenMotherPhi);
    // tree->Branch("tGenMotherPt", &in->tGenMotherPt);
    // tree->Branch("tGenPdgId", &in->tGenPdgId);
    tree->Branch("tGenPhi", &in->tGenPhi);
    tree->Branch("tGenPt", &in->tGenPt);
    // tree->Branch("tGenStatus", &in->tGenStatus);
    // tree->Branch("tJetArea", &in->tJetArea);
    // tree->Branch("tJetBtag", &in->tJetBtag);
    // tree->Branch("tJetDR", &in->tJetDR);
    // tree->Branch("tJetEtaEtaMoment", &in->tJetEtaEtaMoment);
    // tree->Branch("tJetEtaPhiMoment", &in->tJetEtaPhiMoment);
    // tree->Branch("tJetEtaPhiSpread", &in->tJetEtaPhiSpread);
    // tree->Branch("tJetHadronFlavour", &in->tJetHadronFlavour);
    // tree->Branch("tJetPFCISVBtag", &in->tJetPFCISVBtag);
    // tree->Branch("tJetPartonFlavour", &in->tJetPartonFlavour);
    // tree->Branch("tJetPhiPhiMoment", &in->tJetPhiPhiMoment);
    // tree->Branch("tJetPt", &in->tJetPt);
    // tree->Branch("tL1IsoTauMatch", &in->tL1IsoTauMatch);
    // tree->Branch("tL1IsoTauPt", &in->tL1IsoTauPt);
    // tree->Branch("tLeadTrackPt", &in->tLeadTrackPt);
    tree->Branch("tLooseDeepTau2017v2p1VSe", &in->tLooseDeepTau2017v2p1VSe);
    tree->Branch("tLooseDeepTau2017v2p1VSjet", &in->tLooseDeepTau2017v2p1VSjet);
    tree->Branch("tLooseDeepTau2017v2p1VSmu", &in->tLooseDeepTau2017v2p1VSmu);
    tree->Branch("tLowestMll", &in->tLowestMll);
    tree->Branch("tMass", &in->tMass);
    tree->Branch("tMatchesEle24HPSTau30Filter", &in->tMatchesEle24HPSTau30Filter);
    tree->Branch("tMatchesEle24HPSTau30Path", &in->tMatchesEle24HPSTau30Path);
    tree->Branch("tMatchesEle24Tau30Filter", &in->tMatchesEle24Tau30Filter);
    tree->Branch("tMatchesEle24Tau30Path", &in->tMatchesEle24Tau30Path);
    // tree->Branch("tMatchesIsoMu19Tau20Filter", &in->tMatchesIsoMu19Tau20Filter);
    // tree->Branch("tMatchesIsoMu19Tau20Path", &in->tMatchesIsoMu19Tau20Path);
    // tree->Branch("tMatchesIsoMu19Tau20SingleL1Filter", &in->tMatchesIsoMu19Tau20SingleL1Filter);
    // tree->Branch("tMatchesIsoMu19Tau20SingleL1Path", &in->tMatchesIsoMu19Tau20SingleL1Path);
    // tree->Branch("tMatchesIsoMu20HPSTau27Filter", &in->tMatchesIsoMu20HPSTau27Filter);
    // tree->Branch("tMatchesIsoMu20HPSTau27Path", &in->tMatchesIsoMu20HPSTau27Path);
    // tree->Branch("tMatchesIsoMu20Tau27Filter", &in->tMatchesIsoMu20Tau27Filter);
    // tree->Branch("tMatchesIsoMu20Tau27Path", &in->tMatchesIsoMu20Tau27Path);
    tree->Branch("tMatchEmbeddedFilterEle24Tau30", &in->tMatchEmbeddedFilterEle24Tau30);
    tree->Branch("tMediumDeepTau2017v2p1VSe", &in->tMediumDeepTau2017v2p1VSe);
    tree->Branch("tMediumDeepTau2017v2p1VSjet", &in->tMediumDeepTau2017v2p1VSjet);
    tree->Branch("tMediumDeepTau2017v2p1VSmu", &in->tMediumDeepTau2017v2p1VSmu);
    // tree->Branch("tNChrgHadrIsolationCands", &in->tNChrgHadrIsolationCands);
    // tree->Branch("tNChrgHadrSignalCands", &in->tNChrgHadrSignalCands);
    // tree->Branch("tNGammaSignalCands", &in->tNGammaSignalCands);
    // tree->Branch("tNNeutralHadrSignalCands", &in->tNNeutralHadrSignalCands);
    // tree->Branch("tNSignalCands", &in->tNSignalCands);
    // tree->Branch("tNearestZMass", &in->tNearestZMass);
    // tree->Branch("tNeutralIsoPtSum", &in->tNeutralIsoPtSum);
    // tree->Branch("tNeutralIsoPtSumWeight", &in->tNeutralIsoPtSumWeight);
    // tree->Branch("tNeutralIsoPtSumWeightdR03", &in->tNeutralIsoPtSumWeightdR03);
    // tree->Branch("tNeutralIsoPtSumdR03", &in->tNeutralIsoPtSumdR03);
    tree->Branch("tPVDXY", &in->tPVDXY);
    tree->Branch("tPVDZ", &in->tPVDZ);
    tree->Branch("tPhi", &in->tPhi);
    // tree->Branch("tPhotonPtSumOutsideSignalCone", &in->tPhotonPtSumOutsideSignalCone);
    // tree->Branch("tPhotonPtSumOutsideSignalConedR03", &in->tPhotonPtSumOutsideSignalConedR03);
    tree->Branch("tPt", &in->tPt);
    // tree->Branch("tPuCorrPtSum", &in->tPuCorrPtSum);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTLoose", &in->tRerunMVArun2v2DBoldDMwLTLoose);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTMedium", &in->tRerunMVArun2v2DBoldDMwLTMedium);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTTight", &in->tRerunMVArun2v2DBoldDMwLTTight);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTVLoose", &in->tRerunMVArun2v2DBoldDMwLTVLoose);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTVTight", &in->tRerunMVArun2v2DBoldDMwLTVTight);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTVVLoose", &in->tRerunMVArun2v2DBoldDMwLTVVLoose);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTVVTight", &in->tRerunMVArun2v2DBoldDMwLTVVTight);
    // tree->Branch("tRerunMVArun2v2DBoldDMwLTraw", &in->tRerunMVArun2v2DBoldDMwLTraw);
    tree->Branch("tTightDeepTau2017v2p1VSe", &in->tTightDeepTau2017v2p1VSe);
    tree->Branch("tTightDeepTau2017v2p1VSjet", &in->tTightDeepTau2017v2p1VSjet);
    tree->Branch("tTightDeepTau2017v2p1VSmu", &in->tTightDeepTau2017v2p1VSmu);
    tree->Branch("tVLooseDeepTau2017v2p1VSe", &in->tVLooseDeepTau2017v2p1VSe);
    tree->Branch("tVLooseDeepTau2017v2p1VSjet", &in->tVLooseDeepTau2017v2p1VSjet);
    tree->Branch("tVLooseDeepTau2017v2p1VSmu", &in->tVLooseDeepTau2017v2p1VSmu);
    tree->Branch("tVTightDeepTau2017v2p1VSe", &in->tVTightDeepTau2017v2p1VSe);
    tree->Branch("tVTightDeepTau2017v2p1VSjet", &in->tVTightDeepTau2017v2p1VSjet);
    tree->Branch("tVTightDeepTau2017v2p1VSmu", &in->tVTightDeepTau2017v2p1VSmu);
    tree->Branch("tVVLooseDeepTau2017v2p1VSe", &in->tVVLooseDeepTau2017v2p1VSe);
    tree->Branch("tVVLooseDeepTau2017v2p1VSjet", &in->tVVLooseDeepTau2017v2p1VSjet);
    tree->Branch("tVVLooseDeepTau2017v2p1VSmu", &in->tVVLooseDeepTau2017v2p1VSmu);
    tree->Branch("tVVTightDeepTau2017v2p1VSe", &in->tVVTightDeepTau2017v2p1VSe);
    tree->Branch("tVVTightDeepTau2017v2p1VSjet", &in->tVVTightDeepTau2017v2p1VSjet);
    tree->Branch("tVVTightDeepTau2017v2p1VSmu", &in->tVVTightDeepTau2017v2p1VSmu);
    tree->Branch("tVVVLooseDeepTau2017v2p1VSe", &in->tVVVLooseDeepTau2017v2p1VSe);
    tree->Branch("tVVVLooseDeepTau2017v2p1VSjet", &in->tVVVLooseDeepTau2017v2p1VSjet);
    tree->Branch("tVVVLooseDeepTau2017v2p1VSmu", &in->tVVVLooseDeepTau2017v2p1VSmu);
    tree->Branch("tVZ", &in->tVZ);
    tree->Branch("tZTTGenDR", &in->tZTTGenDR);
    tree->Branch("tZTTGenEta", &in->tZTTGenEta);
    tree->Branch("tZTTGenMatching", &in->tZTTGenMatching);
    tree->Branch("tZTTGenPhi", &in->tZTTGenPhi);
    tree->Branch("tZTTGenPt", &in->tZTTGenPt);
    // tree->Branch("tauVetoPt20Loose3HitsVtx", &in->tauVetoPt20Loose3HitsVtx);
    // tree->Branch("tauVetoPt20TightMVALTVtx", &in->tauVetoPt20TightMVALTVtx);
    tree->Branch("topQuarkPt1", &in->topQuarkPt1);
    tree->Branch("topQuarkPt2", &in->topQuarkPt2);
    // tree->Branch("tripleEPass", &in->tripleEPass);
    // tree->Branch("tripleMu10_5_5Pass", &in->tripleMu10_5_5Pass);
    // tree->Branch("tripleMu12_10_5Pass", &in->tripleMu12_10_5Pass);
    tree->Branch("vbfDeta", &in->vbfDeta);
    tree->Branch("vbfJetVeto20", &in->vbfJetVeto20);
    tree->Branch("vbfJetVeto30", &in->vbfJetVeto30);
    tree->Branch("vbfMass", &in->vbfMass);
    tree->Branch("vbfMass_JERDown", &in->vbfMass_JERDown);
    tree->Branch("vbfMass_JERUp", &in->vbfMass_JERUp);
    tree->Branch("vbfMass_JetAbsoluteDown", &in->vbfMass_JetAbsoluteDown);
    tree->Branch("vbfMass_JetAbsoluteUp", &in->vbfMass_JetAbsoluteUp);
    tree->Branch("vbfMass_JetAbsoluteyearDown", &in->vbfMass_JetAbsoluteyearDown);
    tree->Branch("vbfMass_JetAbsoluteyearUp", &in->vbfMass_JetAbsoluteyearUp);
    tree->Branch("vbfMass_JetBBEC1Down", &in->vbfMass_JetBBEC1Down);
    tree->Branch("vbfMass_JetBBEC1Up", &in->vbfMass_JetBBEC1Up);
    tree->Branch("vbfMass_JetBBEC1yearDown", &in->vbfMass_JetBBEC1yearDown);
    tree->Branch("vbfMass_JetBBEC1yearUp", &in->vbfMass_JetBBEC1yearUp);
    tree->Branch("vbfMass_JetEC2Down", &in->vbfMass_JetEC2Down);
    tree->Branch("vbfMass_JetEC2Up", &in->vbfMass_JetEC2Up);
    tree->Branch("vbfMass_JetEC2yearDown", &in->vbfMass_JetEC2yearDown);
    tree->Branch("vbfMass_JetEC2yearUp", &in->vbfMass_JetEC2yearUp);
    tree->Branch("vbfMass_JetEnDown", &in->vbfMass_JetEnDown);
    tree->Branch("vbfMass_JetEnUp", &in->vbfMass_JetEnUp);
    tree->Branch("vbfMass_JetFlavorQCDDown", &in->vbfMass_JetFlavorQCDDown);
    tree->Branch("vbfMass_JetFlavorQCDUp", &in->vbfMass_JetFlavorQCDUp);
    tree->Branch("vbfMass_JetHFDown", &in->vbfMass_JetHFDown);
    tree->Branch("vbfMass_JetHFUp", &in->vbfMass_JetHFUp);
    tree->Branch("vbfMass_JetHFyearDown", &in->vbfMass_JetHFyearDown);
    tree->Branch("vbfMass_JetHFyearUp", &in->vbfMass_JetHFyearUp);
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

    tree->Branch("sm_weight_nlo", &in->sm_weight_nlo);
    tree->Branch("ps_weight_nlo", &in->ps_weight_nlo);
    tree->Branch("mm_weight_nlo", &in->mm_weight_nlo);

    tree->Branch("idx", &in->idx);
}

#endif  // ROOT_SRC_ETAU_TREE_2018_H_
