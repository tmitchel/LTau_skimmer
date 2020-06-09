// Copyright Tyler Mitchell [2020]

#ifndef ROOT_INTERFACE_TAUFESTOOL_H_
#define ROOT_INTERFACE_TAUFESTOOL_H_

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "TFile.h"
#include "TGraph.h"
#include "TH1F.h"

typedef std::vector<Float_t> fvector;
std::vector<Float_t> percentage_to_weight(std::vector<Float_t>);
std::vector<Float_t> percentage_to_decimal(std::vector<Float_t>);

class TauFESTool {
   public:
    TauFESTool(std::string, std::string, std::string, std::string, bool);
    Float_t getFES(Float_t, Float_t, Float_t, std::string);
    Float_t getTES(Float_t, Float_t, Float_t, std::string);

   private:
    std::string year, id, path, channel;
    Float_t low_pt, high_pt;
    bool isEmbed;
    std::unordered_map<std::string, fvector> mc_fes_sfs;
    std::unordered_map<std::string, fvector> mc_tes_sfs, mc_tes_highpt_sfs, mc_tes_int_sfs;
    std::unordered_map<std::string, fvector> embed_fes_sfs;
    std::unordered_map<std::string, fvector> embed_tes_sfs;
    std::unordered_map<int, int> dm_map;
};

TauFESTool::TauFESTool(std::string _year, std::string _id, std::string _path, std::string _channel, bool _isEmbed)
    : year(_year), id(_id), path(_path), channel(_channel), low_pt(34.), high_pt(170.), isEmbed(_isEmbed),
      mc_fes_sfs{{"nominal", fvector{}}, {"up", fvector{}}, {"down", fvector{}}},
      mc_tes_sfs{{"nominal", fvector{}}, {"up", fvector{}}, {"down", fvector{}}},
      mc_tes_highpt_sfs{{"nominal", fvector{}}, {"up", fvector{}}, {"down", fvector{}}},
      mc_tes_int_sfs{{"nominal", fvector{}}, {"up", fvector{}}, {"down", fvector{}}},
      embed_fes_sfs{{"nominal", fvector{}}, {"up", fvector{}}, {"down", fvector{}}},
      embed_tes_sfs{{"nominal", fvector{}}, {"up", fvector{}}, {"down", fvector{}}},
      dm_map{{0, 0}, {1, 1}, {10, 2}, {11, 3}} {
    std::string fes_filepath = "$CMSSW_BASE/src/" + path + "TauFES_eta-dm_" + id + "VSe_" + year + ".root";
    auto fes_file = new TFile(fes_filepath.c_str(), "READ");
    TGraph* fes_graph = reinterpret_cast<TGraph*>(fes_file->Get("fes"));

    std::string tes_filepath = "$CMSSW_BASE/src/" + path + "TauES_dm_" + id + "VSjet_" + year + ".root";
    auto tes_file = new TFile(tes_filepath.c_str(), "READ");
    TH1F* tes_hist = reinterpret_cast<TH1F*>(tes_file->Get("tes"));

    std::string tes_highpt_filepath = "$CMSSW_BASE/src/" + path + "TauES_dm_" + id + "VSjet_" + year + "_ptgt100.root";
    auto tes_highpt_file = new TFile(tes_highpt_filepath.c_str(), "READ");
    TH1F* tes_highpt_hist = reinterpret_cast<TH1F*>(tes_highpt_file->Get("tes"));

    // lepton -> tau fake energy scale corrections (not embedded)
    if (channel == "et") {
        // i = 0 : barrel dm 0
        // i = 1 : barrel dm 1
        // i = 2 : endcap dm 0
        // i = 3 : endcap dm 1
        for (auto i = 0; i < 4; i++) {
            mc_fes_sfs["nominal"].push_back(fes_graph->GetY()[i]);
            mc_fes_sfs["up"].push_back(fes_graph->GetErrorYhigh(i) - fes_graph->GetY()[i]);
            mc_fes_sfs["down"].push_back(fes_graph->GetErrorYlow(i) - fes_graph->GetY()[i]);
        }
    } else if (channel == "mt") {
        // mutau has no correction, but 1% uncorrelated in DM uncertainty
        mc_fes_sfs["nominal"]  = {1.0,  1.0,  1.0,  1.0};
        mc_fes_sfs["up"]       = percentage_to_decimal({ 1.0,  1.0,  1.0,  1.0});
        mc_fes_sfs["down"]     = percentage_to_decimal({-1.0, -1.0, -1.0, -1.0});
    } else {
        std::cerr << "Don't know the provided channel " << channel << std::endl;
        return;
    }

    // lepton -> tau fake energy scale corrections (embedded)
    if (channel == "et") {  // {eta < 1.479, eta > 1.479}
        if (year == "2016Legacy") {
            embed_fes_sfs["nominal"]    = percentage_to_weight({-0.24, -0.70});
        } else if (year == "2017ReReco") {
            embed_fes_sfs["nominal"]    = percentage_to_weight({-0.07, -1.13});
        } else if (year == "2018ReReco") {
            embed_fes_sfs["nominal"]    = percentage_to_weight({-0.33, -0.56});
        }
        // systematics are the same for all years
        embed_fes_sfs["up"]         = percentage_to_decimal({ 0.50,  1.25});
        embed_fes_sfs["down"]       = percentage_to_decimal({-0.50, -1.25});
    } else if (channel == "mt") {
        // mutau has no correction and no uncertainty
        embed_fes_sfs["nominal"]  = {1.0, 1.0};
        embed_fes_sfs["up"]       = {0.0, 0.0};
        embed_fes_sfs["down"]     = {0.0, 0.0};
    } else {
        std::cerr << "Don't know the provided channel " << channel << std::endl;
        return;
    }

    // genuine tau energy scale corrections (not embedded)
    std::vector<int> temp_dms = {0, 1, 10, 11};
    for (auto dm : temp_dms) {
        auto bin = dm + 1;  // dm + 1 bin offset in histogram
        // nominal weight always comes from low pT measurement
        auto nom = tes_hist->GetBinContent(bin);
        mc_tes_sfs["nominal"].push_back(nom);
        mc_tes_highpt_sfs["nominal"].push_back(nom);

        // need to store shift / nominal to apply later
        mc_tes_sfs["up"].push_back(tes_hist->GetBinError(bin) / nom);
        mc_tes_highpt_sfs["up"].push_back(tes_highpt_hist->GetBinError(bin) / nom);
        // just slope for interpolation later
        mc_tes_int_sfs["up"].push_back((tes_highpt_hist->GetBinError(bin) - tes_hist->GetBinError(bin)) / (high_pt - low_pt));

        mc_tes_sfs["down"].push_back(-1. * tes_hist->GetBinError(bin) / nom);
        mc_tes_highpt_sfs["down"].push_back(-1. * tes_highpt_hist->GetBinError(bin) / nom);
        // just slope for interpolation later
        mc_tes_int_sfs["down"].push_back((tes_highpt_hist->GetBinError(bin) - tes_hist->GetBinError(bin)) / (high_pt - low_pt));
    }

    // genuine tau energy scale corrections (embedded)
    if (year == "2016Legacy") {
        embed_tes_sfs["nominal"]    = percentage_to_weight({ -0.20, -0.22, -1.26, -1.26});
        embed_tes_sfs["up"]         = percentage_to_decimal({ 0.46,  0.22,  0.33,  0.33});
        embed_tes_sfs["down"]       = percentage_to_decimal({-0.46, -0.25, -0.51, -0.51});
    } else if (year == "2017ReReco") {
        embed_tes_sfs["nominal"]    = percentage_to_weight({ -0.04, -1.20, -0.75, -0.75});
        embed_tes_sfs["up"]         = percentage_to_decimal({ 0.41,  0.52,  0.44,  0.44});
        embed_tes_sfs["down"]       = percentage_to_decimal({-0.42, -0.21, -0.46, -0.46});
    } else if (year == "2018ReReco") {
        embed_tes_sfs["nominal"]    = percentage_to_weight({ -0.33, -0.57, -0.74, -0.74});
        embed_tes_sfs["up"]         = percentage_to_decimal({ 0.39,  0.37,  0.32,  0.32});
        embed_tes_sfs["down"]       = percentage_to_decimal({-0.39, -0.31, -0.32, -0.32});
    }
}

Float_t TauFESTool::getTES(Float_t pt, Float_t dm, Float_t gen_match, std::string syst = "") {
    if (gen_match != 5) {
        return syst == "" ? 1. : 0.;  // 0 for syst, 1 for nominal
    } else if (dm == 5 || dm == 6) {
      return syst == "" ? 1. : 0.;  // these will be removed anyways
    }

    // handle small fraction of 1prong + 2pizero
    if (dm == 2) {
      dm = 1;
    }

    auto index = dm_map.at(dm);
    // handle default case
    if (syst == "") {
        syst = "nominal";
    }

    if (isEmbed) {
        return embed_tes_sfs.at(syst).at(index);
    }

    if (syst != "nominal") {
        // handle pT dependence of systematic
        if (pt <= 34) {
            return mc_tes_sfs.at(syst).at(index);
        } else if (pt >= 170) {
            return mc_tes_highpt_sfs.at(syst).at(index);
        } else {
            auto shift = mc_tes_sfs.at(syst).at(index) + (mc_tes_int_sfs.at(syst).at(index) / (pt - low_pt));
            shift /= mc_tes_sfs.at("nominal").at(index);
            if (syst == "down") {
                shift *= -1;
            }
            return shift;
        }
    }

    return mc_tes_sfs.at("nominal").at(index);
}

Float_t TauFESTool::getFES(Float_t dm, Float_t eta, Float_t gen_match, std::string syst = "") {
    if (gen_match > 4) {  // don't apply to jets or genuine taus
        return syst == "" || syst == "nomnial" ? 1. : 0.;
    } else if (dm != 0 && dm != 1) {
      return syst == "" || syst == "nomnial" ? 1. : 0.;  // these will be removed anyways
    }

    // handle default case
    if (syst == "") {
        syst = "nominal";
    }

    if (isEmbed) {
        if (fabs(eta) < 1.479) {
            return embed_fes_sfs[syst].at(0);
        } else {
            return embed_fes_sfs[syst].at(1);
        }
    }

    // handle barrel vs endcap
    int index = 0;
    if (fabs(eta) > 1.479) {
        index = 2;
    }

    // decay mode 0 or 1
    index += dm;

    return mc_fes_sfs[syst].at(index);
}

std::vector<Float_t> percentage_to_weight(std::vector<Float_t> input) {
    std::vector<Float_t> output;
    for (auto in : input) {
        output.push_back(1 + (0.01 * in));
    }
    return output;
}

std::vector<Float_t> percentage_to_decimal(std::vector<Float_t> input) {
    std::vector<Float_t> output;
    for (auto in : input) {
        output.push_back(0.01 * in);
    }
    return output;
}

#endif  // ROOT_INTERFACE_TAUFESTOOL_H_
