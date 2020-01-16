// Copyright 2019 Tyler Mitchell

#ifndef ROOT_SRC_SYNC_MUTAU_TREE_2018_H_
#define ROOT_SRC_SYNC_MUTAU_TREE_2018_H_

#include <cmath>
#include <iostream>
#include <utility>
#include <vector>
#include "./base_tree.h"
#include "RooFunctor.h"
#include "RooMsgService.h"
#include "RooRealVar.h"
#include "RooWorkspace.h"
#include "data/LumiReweightingStandAlone.h"
#include "HTT-utilities/RecoilCorrections/interface/MEtSys.h"
#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "ltau_skimmer/ROOT/interface/mutau_input_branches.h"

class sync_mutau_tree2018 : public virtual base_tree {
   private:
    TTree *tree, *original;
    mutau_input_branches* in;
    bool isMC, isEmbed, isData;
    std::vector<Int_t> good_events;
    TLorentzVector mu, tau, MET;


   public:
    // Member variables
    UInt_t Run, Lumi;
    Int_t recoil;
    Float_t placeholder;  // for all branches not present in 2018

    // // Constructed while running
    Bool_t trg_singlemuon, trg_mutaucross, flagFilter;
    Float_t puweight, trigweight_1, idisoweight_1, idisoweight_2, trackingweight_1;

    Int_t era;
    Int_t gen_match_1, gen_match_2, njets, nbtag, njetspt20;
    Float_t tes_dm0_sf, tes_dm1_sf, tes_dm10_sf, efake_dm0_sf, efake_dm1_sf, mfake_dm0_sf, mfake_dm1_sf;
    Float_t jetVeto20, jetVeto30, met, metphi, met_px, met_py, extraelec_veto, extramuon_veto, dilepton_veto, pfmetcorr_ex, pfmetcorr_ey;
    Float_t pt_1, eta_1, phi_1, m_1, e_1, px_1, py_1, pz_1, pt_2, eta_2, phi_2, m_2, e_2, px_2, py_2, pz_2;

    std::vector<TLorentzVector*> mets;

    // Member functions
    sync_mutau_tree2018(TTree* orig, TTree* itree, bool isMC, bool isEmbed, Int_t rec);
    virtual ~sync_mutau_tree2018() {}
    void do_skimming(TH1F*);
    void set_branches();
    Float_t get_tes_sf(Float_t);
    Float_t get_efake_sf(Float_t);
    Float_t get_mfake_sf(Float_t);
    void do_met_corr_nom(Float_t, energy_scale, TLorentzVector, TLorentzVector*);
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
sync_mutau_tree2018::sync_mutau_tree2018(TTree* Original, TTree* itree, bool IsMC, bool IsEmbed, Int_t rec)
    : tree(itree),
      original(Original),
      in(new mutau_input_branches(Original)),
      isMC(IsMC),
      isEmbed(IsEmbed),
      recoil(rec),
      era(2018),
      tes_dm0_sf(0.987),
      tes_dm1_sf(0.995),
      tes_dm10_sf(0.988),
      efake_dm0_sf(0.968),
      efake_dm1_sf(1.026),
      mfake_dm0_sf(0.998),
      mfake_dm1_sf(0.990) {
    // set embedded TES
    if (isEmbed) {
        tes_dm0_sf = 0.975;
        tes_dm1_sf = 0.975 * 1.051;
        tes_dm10_sf = 0.975 * 0.975 * 0.975;
        efake_dm0_sf = 1.;
        efake_dm1_sf = 1.;
        mfake_dm0_sf = 1.;
        mfake_dm1_sf = 1.;
    }
}

//////////////////////////////////////////////////////////////////
// Purpose: Skim original then apply Isolation-based sorting.   //
//          Good events will be placed in the good_events       //
//          vector for later                                    //
//////////////////////////////////////////////////////////////////
void sync_mutau_tree2018::do_skimming(TH1F* cutflow) {
    // declare variables for sorting
    ULong64_t evt_now(0);
    ULong64_t evt_before(1);
    int best_evt(-1);
    std::pair<float, float> muCandidate, tauCandidate;

    std::cout << "Starting the skim..." << std::endl;

    isData = !isEmbed && !isMC;
    Int_t nevt = (Int_t)original->GetEntries();
    for (auto ievt = 0; ievt < nevt; ievt++) {
        original->GetEntry(ievt);
        evt_now = in->evt;

        // TLorentzVector ele, tau;
        mu.SetPtEtaPhiM(in->mPt, in->mEta, in->mPhi, in->mMass);
        tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);

        // apply TES
        if (isMC || isEmbed) {
            if (in->tZTTGenMatching == 5) {
                tau *= get_tes_sf(in->tDecayMode);
            } else if (in->tZTTGenMatching == 1 || in->tZTTGenMatching == 3) {
                tau *= get_efake_sf(in->tDecayMode);
            } else if (in->tZTTGenMatching == 2 || in->tZTTGenMatching == 4) {
                tau *= get_mfake_sf(in->tDecayMode);
            }
        }

        cutflow->Fill(1., 1.);
        // apply event selection

        if (mu.Pt() > 21. && fabs(mu.Eta()) < 2.1 && fabs(in->mPVDZ) < 0.2 && fabs(in->mPVDXY) < 0.045)
            cutflow->Fill(3., 1.);  // electron kinematic selection
        else
            continue;

        if (in->mPFIDMedium)
            cutflow->Fill(4., 1.);  // muon quality selection
        else
            continue;

        if (tau.Pt() > 20. && fabs(tau.Eta()) < 2.3 && fabs(in->tPVDZ) < 0.2)
            cutflow->Fill(6., 1.);  // tau kinematic selection
        else
            continue;

        if (in->tRerunMVArun2v2DBoldDMwLTVVLoose && in->tDecayModeFinding > 0.5 && fabs(in->tCharge) < 2)
            cutflow->Fill(7., 1.);  // tau quality selection
        else
            continue;

        if (mu.DeltaR(tau) > 0.5) {
            cutflow->Fill(10., 1.);
        } else {
            continue;
        }

        // implement new sorting per
        // https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2017#Baseline_Selection
        if (evt_now != evt_before) {  // new event, save the tau candidates
            // since it is new event, do we have the best entry to save? If yes, save it!
            if (best_evt > -1) good_events.push_back(best_evt);

            //  this is a new event, so the first tau pair is the best! :)
            best_evt = ievt;
            muCandidate = std::make_pair(mu.Pt(), in->mRelPFIsoDBDefault);
            tauCandidate = std::make_pair(tau.Pt(), in->tRerunMVArun2v2DBoldDMwLTraw);
        } else {  // not a new event
            std::pair<float, float> currEleCandidate(mu.Pt(), in->mRelPFIsoDBDefault);
            std::pair<float, float> currTauCandidate(tau.Pt(), in->tRerunMVArun2v2DBoldDMwLTraw);

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
                }      // tau1 has the same pT
            }          // tau1 has the same isolation
        }              // not a new event
        evt_before = evt_now;
    }
    if (best_evt > -1) good_events.push_back(best_evt);
}

Float_t sync_mutau_tree2018::get_tes_sf(Float_t decayMode) {
    if (decayMode == 0) {
        return tes_dm0_sf;
    } else if (decayMode == 1) {
        return tes_dm1_sf;
    } else if (decayMode == 10) {
        return tes_dm10_sf;
    }
    return 1.;
}

Float_t sync_mutau_tree2018::get_efake_sf(Float_t decayMode) {
    if (decayMode == 0) {
        return efake_dm0_sf;
    } else if (decayMode == 1) {
        return efake_dm1_sf;
    }
    return 1.;
}

Float_t sync_mutau_tree2018::get_mfake_sf(Float_t decayMode) {
    if (decayMode == 0) {
        return mfake_dm0_sf;
    } else if (decayMode == 1) {
        return mfake_dm1_sf;
    }
    return 1.;
}

void sync_mutau_tree2018::do_met_corr_nom(Float_t decayMode, energy_scale escale, TLorentzVector tau, TLorentzVector* met) {
    double sf(1.);
    if (escale == tes) {
        sf = get_tes_sf(decayMode);
    } else if (escale == efake) {
        sf = get_efake_sf(decayMode);
    } else if (escale == mfake) {
        sf = get_mfake_sf(decayMode);
    } else {
        std::cerr << "Not a valid energy correction" << std::endl;
    }
    *met = (*met) + tau - sf * tau;  // update input met}
}

void sync_mutau_tree2018::do_recoil_corr(RecoilCorrector* recoilPFMetCorrector, TLorentzVector* met, int njets) {
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
TTree* sync_mutau_tree2018::fill_tree(RecoilCorrector recoilPFMetCorrector, MEtSys metSys) {
    std::cout << "setting branches..." << std::endl;
    set_branches();  // get all the branches set up
    std::cout << "branches set." << std::endl;

    // legacy sf's
    TFile htt_sf_file("data/htt_scalefactors_legacy_2018.root");
    RooWorkspace *htt_sf = reinterpret_cast<RooWorkspace*>(htt_sf_file.Get("w"));
    htt_sf_file.Close();

    auto lumi_weights =
        new reweight::LumiReWeighting("data/pu_distributions_mc_2018.root", "data/pu_distributions_data_2018.root", "pileup", "pileup");

    // loop through all events pasing skimming/sorting
    for (auto& ievt : good_events) {
        original->GetEntry(ievt);

        // TLorentzVector mu, tau;
        mu.SetPtEtaPhiM(in->mPt, in->mEta, in->mPhi, in->mMass);
        tau.SetPtEtaPhiM(in->tPt, in->tEta, in->tPhi, in->tMass);
        met_px = in->type1_pfMetEt * cos(in->type1_pfMetPhi);
        met_py = in->type1_pfMetEt * sin(in->type1_pfMetPhi);
        MET.SetPtEtaPhiM(in->type1_pfMetEt, 0, in->type1_pfMetPhi, 0);

        mets = {&MET};
        auto jet_for_correction = in->jetVeto30;
        if (recoil == 1) {
            jet_for_correction += 1;
        }

        // do recoil corrections on all met
        for (unsigned i = 0; i < mets.size(); i++) {
            do_recoil_corr(&recoilPFMetCorrector, mets.at(i), jet_for_correction);
        }

        if (isMC || isEmbed) {
            // met correction due to tau energy scale
            if (in->tZTTGenMatching == 5) {
                for (unsigned i = 0; i < mets.size(); i++) {
                    do_met_corr_nom(in->tDecayMode, tes, tau, mets.at(i));
                }
                tau *= get_tes_sf(in->tDecayMode);
            } else if (in->tZTTGenMatching == 1 || in->tZTTGenMatching == 3) {
                // electron -> tau fake energy scale
                for (unsigned i = 0; i < mets.size(); i++) {
                    do_met_corr_nom(in->tDecayMode, efake, tau, mets.at(i));
                }
                tau *= get_efake_sf(in->tDecayMode);
            } else if (in->tZTTGenMatching == 2 || in->tZTTGenMatching == 4) {
                // muon -> tau fake energy scale
                for (unsigned i = 0; i < mets.size(); i++) {
                    do_met_corr_nom(in->tDecayMode, mfake, tau, mets.at(i));
                }
                tau *= get_mfake_sf(in->tDecayMode);
            }
        }

        auto Mu24 = in->IsoMu24Pass && in->mMatchesIsoMu24Path && in->mMatchesIsoMu24Filter;
        auto Mu27 = in->IsoMu27Pass && in->mMatchesIsoMu27Path && in->mMatchesIsoMu27Filter;
        auto Cross_base = in->mMatchesIsoMu20HPSTau27Filter && in->mMatchesIsoMu20HPSTau27Path && in->tMatchesIsoMu20HPSTau27Filter &&
                          in->tMatchesIsoMu20HPSTau27Path;
        auto Cross_v1 = Cross_base && in->Mu20LooseHPSTau27Pass && in->run < 317509 && isData;
        auto Cross_v2 = Cross_base && in->Mu20LooseHPSTau27TightIDPass && (isMC || (in->run > 317509 && isData));
        auto Mu24_emb = in->mMatchEmbeddedFilterMu24;
        auto Mu27_emb = in->mMatchEmbeddedFilterMu27;
        // auto Cross_emb = in->mMatchEmbeddedFilterMu20Tau27_2018 && in->tMatchEmbeddedFilterMu20HPSTau27;

        trg_singlemuon = false;
        trg_mutaucross = false;
        if (isEmbed) {
            if (Mu27_emb && mu.Pt() > 28) {
                trg_singlemuon = true;
            } else if (Mu24_emb && mu.Pt() > 25) {
                trg_singlemuon = true;
            } else if (mu.Pt() > 21 && mu.Pt() < 25 && tau.Pt() > 32 && fabs(tau.Eta()) < 2.1) {
                trg_mutaucross = true;
            } else {
                // nothing
            }
        } else {
            if (Mu27 && mu.Pt() > 28) {
                trg_singlemuon = true;
            } else if (Mu24 && mu.Pt() > 25) {
                trg_singlemuon = true;
            } else if ((Cross_v1 || Cross_v2) && mu.Pt() > 21 && mu.Pt() < 25 && tau.Pt() > 32 && fabs(tau.Eta()) < 2.1) {
                trg_mutaucross = true;
            } else {
                // nothing
            }
        }

        Run = in->run;
        Lumi = in->lumi;
        dilepton_veto = in->dimuonVeto > 0;
        extraelec_veto = in->eVetoZTTp001dxyzR0 > 0;
        extramuon_veto = in->muVetoZTTp001dxyzR0 > 1;

        puweight = lumi_weights->weight(in->nvtx);
        if (trg_singlemuon) {
            trigweight_1 = htt_sf->function("m_trg_ic_ratio")->getVal();
        } else if (trg_mutaucross) {
            trigweight_1 = htt_sf->function("m_trg_20_ic_ratio")->getVal() * htt_sf->function("t_trg_pog_deeptau_medium_mutau_ratio")->getVal();
        }
        idisoweight_1 = htt_sf->function("m_idiso_ic_ratio")->getVal();
        idisoweight_2 = htt_sf->function("t_deeptauid_pt_medium")->getVal();
        trackingweight_1 = htt_sf->function("m_trk_ratio")->getVal();

        flagFilter = !(in->Flag_goodVertices || in->Flag_globalSuperTightHalo2016Filter || in->Flag_HBHENoiseFilter
                || in->Flag_HBHENoiseIsoFilter || in->Flag_EcalDeadCellTriggerPrimitiveFilter || in->Flag_BadPFMuonFilter
                || in->Flag_eeBadScFilter || in->Flag_ecalBadCalibFilter);

        pt_1 = mu.Pt();
        phi_1 = mu.Phi();
        eta_1 = mu.Eta();
        m_1 = mu.M();

        px_1 = mu.Px();
        py_1 = mu.Py();
        pz_1 = mu.Pz();
        e_1 = mu.E();
        m_2 = tau.M();
        px_2 = tau.Px();
        py_2 = tau.Py();
        pz_2 = tau.Pz();
        e_2 = tau.E();
        pt_2 = tau.Pt();
        phi_2 = tau.Phi();
        eta_2 = tau.Eta();

        // convert from Float_t in FSA to Int_t for analyzer
        gen_match_1 = in->mZTTGenMatching;
        gen_match_2 = in->tZTTGenMatching;
        njets = in->jetVeto30;
        nbtag = in->bjetDeepCSVVeto20Medium_2018_DR0p5;
        njetspt20 = in->jetVeto20;

        met = MET.Pt();
        metphi = MET.Phi();
        met_px = MET.Px();
        met_py = MET.Py();

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
void sync_mutau_tree2018::set_branches() {
    // new branches
    tree->Branch("run", &Run);
    tree->Branch("era", &era);
    tree->Branch("lumi", &Lumi);
    tree->Branch("evt", &in->evt);
    tree->Branch("dilepton_veto", &dilepton_veto);
    tree->Branch("extraelec_veto", &extraelec_veto);
    tree->Branch("extramuon_veto", &extramuon_veto);
    tree->Branch("trg_singlemuon", trg_singlemuon);
    tree->Branch("trg_mutaucross", trg_mutaucross);
    tree->Branch("pt_1", &pt_1);
    tree->Branch("phi_1", &phi_1);
    tree->Branch("eta_1", &eta_1);
    tree->Branch("m_1", &m_1);
    tree->Branch("q_1", &in->mCharge);
    tree->Branch("d0_1", &in->mPVDXY);
    tree->Branch("dZ_1", &in->mPVDZ);
    tree->Branch("iso_1", &in->mRelPFIsoDBDefault);
    tree->Branch("gen_match_1", &gen_match_1);
    tree->Branch("pt_2", &pt_2);
    tree->Branch("phi_2", &phi_2);
    tree->Branch("eta_2", &eta_2);
    tree->Branch("m_2", &m_2);
    tree->Branch("q_2", &in->tCharge);
    tree->Branch("d0_2", &in->tPVDXY);
    tree->Branch("dZ_2", &in->tPVDZ);
    tree->Branch("iso_2", &in->tDeepTau2017v2p1VSjetraw);
    tree->Branch("gen_match_2", &gen_match_2);
    tree->Branch("tDeepTau2017v2p1VSmuraw", &in->tDeepTau2017v2p1VSmuraw);
    tree->Branch("tLooseDeepTau2017v2p1VSmu", &in->tLooseDeepTau2017v2p1VSmu);
    tree->Branch("tMediumDeepTau2017v2p1VSmu", &in->tMediumDeepTau2017v2p1VSmu);
    tree->Branch("tTightDeepTau2017v2p1VSmu", &in->tTightDeepTau2017v2p1VSmu);
    tree->Branch("tVLooseDeepTau2017v2p1VSmu", &in->tVLooseDeepTau2017v2p1VSmu);
    tree->Branch("tVTightDeepTau2017v2p1VSmu", &in->tVTightDeepTau2017v2p1VSmu);
    tree->Branch("tVVLooseDeepTau2017v2p1VSmu", &in->tVVLooseDeepTau2017v2p1VSmu);
    tree->Branch("tVVTightDeepTau2017v2p1VSmu", &in->tVVTightDeepTau2017v2p1VSmu);
    tree->Branch("tVVVLooseDeepTau2017v2p1VSmu", &in->tVVVLooseDeepTau2017v2p1VSmu);
    tree->Branch("met", &met);
    tree->Branch("metphi", &metphi);
    tree->Branch("metcov00", &in->metcov00);
    tree->Branch("metcov01", &in->metcov01);
    tree->Branch("metcov10", &in->metcov10);
    tree->Branch("metcov11", &in->metcov11);
    tree->Branch("mjj", &in->vbfMass);
    tree->Branch("jdeta", &in->vbfDeta);
    tree->Branch("nbtag", &nbtag);
    tree->Branch("njets", &njets);
    tree->Branch("njetspt20", &njetspt20);
    tree->Branch("jpt_1", &in->j1pt);
    tree->Branch("jeta_1", &in->j1eta);
    tree->Branch("jphi_1", &in->j1phi);
    tree->Branch("jpt_2", &in->j2pt);
    tree->Branch("jeta_2", &in->j2eta);
    tree->Branch("jphi_2", &in->j2phi);
    tree->Branch("NUP", &in->NUP);
    tree->Branch("rho", &in->rho);

    tree->Branch("idx", &in->idx);
}

#endif  // ROOT_SRC_SYNC_MUTAU_TREE_2018_H_
