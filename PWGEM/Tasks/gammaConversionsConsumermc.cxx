// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \brief perform photon conversion analysis on V0 candidates from aod::StoredV0Datas
/// dependencies: o2-analysis-lf-lambdakzerobuilder
/// \author stephan.friedrich.stiefelmaier@cern.ch

#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/runDataProcessing.h"

#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Centrality.h"
// todo: probably move somewhere else
#include "gammaTables.h"
#include "range.hpp"

#include <TH1.h>
#include <TH1F.h>
#include <TVector3.h>

#include <cmath>
#include <iostream>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

// using collisionEvSelIt = soa::Join<aod::Collisions, aod::EvSels>::iterator;
struct GammaConversionsConsumermc {

  Configurable<bool> fDoEventSel{"fDoEventSel", 0, "demand sel7 for events"};

  Configurable<float> fCentMin{"fCentMin", 0.0, "lower bound of centrality selection"};
  Configurable<float> fCentMax{"fCentMax", 100.0, "upper bound of centrality selection"};

  Configurable<float> fSinglePtCut{"fSinglePtCut", 0.04, "minimum daughter track pt"};

  Configurable<float> fEtaCut{"fEtaCut", 0.8, "accepted eta range"};

  Configurable<float> fMinR{"fMinR", 5., "minimum conversion radius of the V0s"};
  Configurable<float> fMaxR{"fMaxR", 180., "maximum conversion radius of the V0s"};

  Configurable<float> fPIDnSigmaBelowElectronLine{"fPIDnSigmaBelowElectronLine", -3., "minimum sigma electron PID for V0 daughter tracks"};

  Configurable<float> fPIDnSigmaAboveElectronLine{"fPIDnSigmaAboveElectronLine", 3., "maximum sigma electron PID for V0 daughter tracks"};

  Configurable<float> fPIDnSigmaAbovePionLine{"fPIDnSigmaAbovePionLine", 3., "minimum sigma to be over the pion line for low momentum tracks"}; //case 4: 3.0sigma, 1.0 sigma at high momentum

  Configurable<float> fPIDnSigmaAbovePionLineHighP{"fPIDnSigmaAbovePionLineHighP", 1., "minimum sigma to be over the pion line for high momentum tracks"};

  Configurable<float> fPIDMinPnSigmaAbovePionLine{"fPIDMinPnSigmaAbovePionLine", 0.4, "minimum track momentum to apply any pion rejection"}; //case 7:  // 0.4 GeV

  Configurable<float> fPIDMaxPnSigmaAbovePionLine{"fPIDMaxPnSigmaAbovePionLine", 8., "border between low and high momentum pion rejection"}; //case 7:  // 8. GeV

  Configurable<float> fMinTPCFoundOverFindableCls{"fMinTPCNClsFoundOverFindable", 0.6, "minimum ratio found tpc clusters over findable"}; //case 9:  // 0.6

  Configurable<float> fMinTPCCrossedRowsOverFindableCls{"fMinTPCCrossedRowsOverFindableCls", 0.0, "minimum ratio TPC crossed rows over findable clusters"};

  Configurable<float> fQtPtMax{"fQtPtMax", 0.11, "up to fQtMax, multiply the pt of the V0s by this value to get the maximum qt "};

  Configurable<float> fQtMax{"fQtMax", 0.040, "maximum qt"};
  Configurable<float> fMaxPhotonAsymmetry{"fMaxPhotonAsymmetry", 0.95, "maximum photon asymetry"};
  Configurable<float> fPsiPairCut{"fPsiPairCut", 0.1, "maximum psi angle of the track pair"};

  Configurable<float> fCosPAngleCut{"fCosPAngleCut", 0.85, "mimimum cosinus of the pointing angle"}; // case 4

  // SFS todo: struct MyHistogramSpec : public HistogramSpec {
  struct MyHistogramSpec {
    MyHistogramSpec(char const* const name_, char const* const title_, HistogramConfigSpec config_, bool callSumw2_ = false)
      : name(name_),
        hash(compile_time_hash(name_)),
        title(title_),
        config(config_),
        callSumw2(callSumw2_)
    {
    }

    MyHistogramSpec(char const* const name_, HistogramConfigSpec config_, bool callSumw2_ = false)
      : name(name_),
        hash(compile_time_hash(name_)),
        title(name_),
        config(config_),
        callSumw2(callSumw2_)
    {
    }

    MyHistogramSpec()
      : name(""),
        hash(0),
        config()
    {
    }
    MyHistogramSpec(MyHistogramSpec const& other) = default;
    MyHistogramSpec(MyHistogramSpec&& other) = default;

    std::string name{};
    uint32_t hash{};
    std::string title{};
    HistogramConfigSpec config{};
    bool callSumw2{}; // wether or not hist needs heavy error structure produced by Sumw2()
  };

  HistogramRegistry fHistogramRegistry{"fHistogramRegistry"};

  // track histograms
  std::vector<MyHistogramSpec> fTrackHistoDefinitions{
    {"hTPCFoundOverFindableCls", {HistType::kTH1F, {{800, 0.9f, 1.01f}}}},
    {"hTPCCrossedRowsOverFindableCls", {HistType::kTH1F, {{800, 0.8f, 1.5f}}}},

    {"hTPCdEdx", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, 0.f, 200.f}}}},
    {"hTPCdEdxSigEl", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, -10.f, 10.f}}}},
    {"hTPCdEdxSigPi", {HistType::kTH2F, {{800, 0.03f, 20.f}, {800, -10.f, 10.f}}}}};

  // v0 histograms
  std::vector<MyHistogramSpec> fV0HistoDefinitions{
    {"hPt", {HistType::kTH1F, {{800, 0.0f, 25.0f}}}},
    {"hEta", {HistType::kTH1F, {{800, -2.f, 2.f}}}},
    {"hPhi", {HistType::kTH1F, {{800, 0.f, 2.f * M_PI}}}},
    {"hConvPointR", {HistType::kTH1F, {{800, 0.f, 200.f}}}},
    {"hArmenteros", {HistType::kTH2F, {{800, -1.f, 1.f}, {800, 0.f, 0.25f}}}},
    {"hPsiPt", {HistType::kTH2F, {{800, -2.f, 2.f}, {800, 0.f, 10.f}}}},
    {"hCosPAngle", {HistType::kTH1F, {{800, 0.99f, 1.005f}}}},

    {"IsPhotonSelected", {HistType::kTH1F, {{13, -0.0f, 12.5f}}}} // only in afterCuts
  };

  // only in mc
  // resolution histos
  std::vector<MyHistogramSpec> fV0ResolutionHistoDefinitions{
    {"hPtRes", "hPtRes_Rec-MC", {HistType::kTH1F, {{800, -0.5f, 0.5f}}}},
    {"hEtaRes", "hEtaRes_Rec-MC", {HistType::kTH1F, {{800, -0.05f, 0.05f}}}},
    {"hPhiRes", "hPhiRes_Rec-MC", {HistType::kTH1F, {{800, -0.05f, 0.05f}}}},
    {"hConvPointRRes", "hConvPointRRes_Rec-MC", {HistType::kTH1F, {{800, -200.f, 200.f}}}},
    {"hConvPointAbsoluteDistanceRes", "hConvPointAbsoluteDistanceRes", {HistType::kTH1F, {{800, -0.0f, 200.f}}}},
  };

  // Reconstructed info of MC validated V0s histos fV0ReconstructedInfoMCValidatedHistos
  std::vector<MyHistogramSpec> fV0ReconstructedInfoMCValidatedHistoDefinitions{
    {"hValidatedPtRec", {HistType::kTH1F, {{800, 0.0f, 25.f}}}},
    {"hValidatedRRec", {HistType::kTH1F, {{800, 0.0f, 250.f}}}},
  };

  // MC info only histos
  std::vector<MyHistogramSpec> fV0MCInfoOnlyHistoDefinitions{
    {"hValidatedPtTrue", {HistType::kTH1F, {{800, -0.0f, 25.f}}}},
  };
  
  typedef std::map<std::string, HistPtr> mapStringHistPtr;
  mapStringHistPtr fTrackHistos;
  mapStringHistPtr fV0ResolutionHistos;
  mapStringHistPtr fV0ReconstructedInfoMCValidatedHistos;
  
  std::vector<mapStringHistPtr> fRecTrueV0Histos{2};
  enum eRecTrueEnum {kRec, kTrue};
  std::vector<std::string> fRecTrueStrings{"Rec", "True"};
        

  std::string fPrefixReconstructedInfoHistos{"reconstructedInformationOnly/"};
  std::string fPrefixMCInfoNeededHistos{"MCinformationNeeded/"};

  std::string fFullNameIsPhotonSelectedHisto{fPrefixReconstructedInfoHistos + "v0/afterCuts/IsPhotonSelectedRec"};

  std::map<std::string, size_t> fPhotonCutIndeces{
    {"kPhotonIn", 0},
    {"kTrackEta", 1},
    {"kTrackPt", 2},
    {"kElectronPID", 3},
    {"kPionRejLowMom", 4},
    {"kPionRejHighMom", 5},
    {"kTPCFoundOverFindableCls", 6},
    {"kTPCCrossedRowsOverFindableCls", 7},
    {"kV0Radius", 8},
    {"kArmenteros", 9},
    {"kPsiPair", 10},
    {"kCosinePA", 11},
    {"kPhotonOut", 12}};

  //~ lXaxis = registry.get<TH1>(HIST("hEvents"))->GetXaxis();
  //~ lXaxis->SetBinLabel(1, "in");
  //~ lXaxis->SetBinLabel(2, "int7||sel7");
  //~ lXaxis->SetBinLabel(3, "cent");
  //~ lXaxis->SetBinLabel(4, "out");

  void init(InitContext const&)
  {
    for (auto bac : std::vector<std::string>{"beforeCuts/", "afterCuts/"}) {
      
      // reconstructed data histograms
      {
        // track histograms
        {
          std::string lPath(fPrefixReconstructedInfoHistos + "track/" + bac);
          for (auto& tHisto : fTrackHistoDefinitions) {
            std::string lFullName(lPath + tHisto.name);
            LOGF(debug, "adding %s", lFullName);
            fTrackHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          }
        }

        // v0 histograms
        //~ {
          //~ std::string lPath(fPrefixReconstructedInfoHistos + "v0/" + bac);
          //~ for (auto tHisto : fV0HistoDefinitions) {
            //~ if (tHisto.name == "IsPhotonSelected" && bac == "beforeCuts/") {
              //~ continue;
            //~ }
            //~ std::string lFullName(lPath + tHisto.name);
                        //~ LOGF(debug, "adding %s", lFullName);

            //~ fV0Histos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          //~ }
        //~ }
      }
      
      { 
        std::vector<eRecTrueEnum> lRecTrue{kRec, kTrue};
        std::vector<std::string> lPaths{fPrefixReconstructedInfoHistos + "v0/" + bac , fPrefixMCInfoNeededHistos + "v0/MCinformationOnly/" + bac};

        for (auto iRecTrue : lRecTrue){
          //~ auto const &lRecTrueStr = lRecTrueStrings[iRecTrue];
          std::string lPath = lPaths[iRecTrue];
          auto &lV0Histos = fRecTrueV0Histos[iRecTrue];
          
          for (auto tHisto : fV0HistoDefinitions) {
            
            if (tHisto.name == "IsPhotonSelected" && ((iRecTrue == kTrue) || bac == "beforeCuts/")) {
              continue;
            }
            
            std::string lFullName(lPath + tHisto.name + fRecTrueStrings[iRecTrue]);
            LOGF(debug, "adding %s", lFullName);
            LOGF(warning, "adding %s", lFullName); // todo remove this


            lV0Histos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          }
        }
      }
      

      // MC information histos
      {
        // v0 Info only histos
        //~ {
          //~ std::string lPath(fPrefixMCInfoNeededHistos + "v0/MCinformationOnly/" + bac);
          //~ for (auto tHisto : fV0MCInfoOnlyHistoDefinitions) {
            //~ std::string lFullName(lPath + tHisto.name);
            //~ LOGF(debug, "adding %s", lFullName);
            //~ fV0MCInfoOnlyHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          //~ }
        //~ }
        
        // reconstructed info of MC validated V0s histos
        {
          {
            std::string lPath(fPrefixMCInfoNeededHistos + "v0/reconstructedInfoOfMCvalidated/" + bac);
            for (auto tHisto : fV0ReconstructedInfoMCValidatedHistoDefinitions) {
              std::string lFullName(lPath + tHisto.name);
              LOGF(debug, "adding %s", lFullName);
              fV0ReconstructedInfoMCValidatedHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
            }
          }
        }
        
        // v0 Resolution histos
        {
          std::string lPath(fPrefixMCInfoNeededHistos + "v0/resolutions/" + bac);
          for (auto tHisto : fV0ResolutionHistoDefinitions) {
            std::string lFullName(lPath + tHisto.name);
            LOGF(debug, "adding %s", lFullName);
            fV0ResolutionHistos.insert(std::pair{lFullName, fHistogramRegistry.add(lFullName.data(), tHisto.title.data(), tHisto.config)});
          }
        }
      }
    }

    {
      // do some labeling
      auto lIsPhotonSelectedHisto = fRecTrueV0Histos[kRec].find(fFullNameIsPhotonSelectedHisto);
      if (lIsPhotonSelectedHisto != fRecTrueV0Histos[kRec].end()) {
        TAxis* lXaxis = std::get<std::shared_ptr<TH1>>(lIsPhotonSelectedHisto->second)->GetXaxis();
        for (auto& lPairIt : fPhotonCutIndeces) {
          lXaxis->SetBinLabel(lPairIt.second + 1, lPairIt.first.data());
        }
      }
    }
  }

  // SFS todo: think about if this is actually too expensive. Going the other way round with the indices as keys wouldnt require lookups at inserting but pbly produce a but of code duplication at the definition of the cut names
  size_t gMax_size = (size_t)-1;
  size_t getPhotonCutIndex(std::string const& theKey)
  {
    auto lPairIt = fPhotonCutIndeces.find(theKey);
    if (lPairIt != fPhotonCutIndeces.end()) {
      return lPairIt->second;
    }
    return gMax_size;
  }

  template <typename TCOLL, typename TV0, typename TTRACKS>
  bool processPhoton(TCOLL const &theCollision, TV0 const &theV0, TTRACKS const &theV0Tracks)
  {
    auto lV0Tracks = theV0Tracks.sliceBy(aod::v0data::v0Id, theV0.v0Id());
    float lV0CosinePA = theV0.v0cosPA(theCollision.posX(), theCollision.posY(), theCollision.posZ());
    
    fillReconstructedInfoHistograms("beforeCuts/",
                                    theV0,
                                    lV0Tracks,
                                    lV0CosinePA);
                                    
    for (auto &lTrack : theV0Tracks){
      if (!trackPassesCuts(lTrack)) {
        return kFALSE;
      }
    }

    // apply photon cuts
    if (!passesPhotonCuts(theV0, lV0CosinePA)) {
      return kFALSE;
    }

    fillReconstructedInfoHistograms("afterCuts/",
                                    theV0,
                                    lV0Tracks,
                                    lV0CosinePA);
    return kTRUE;
  }

  // todo: maybe combine this with fillV0Histograms
  template <typename TV0, typename TMCGAMMATABLE>
  void fillTruePhotonHistograms(std::string theBAC, TV0 const &theV0, TMCGAMMATABLE const &theTrueGammaTable)
  {
    if (!theTrueGammaTable.size()){
      return;
    }
    auto lTrueGamma = theTrueGammaTable.begin();
    
    // v0 MCInfoOnly histos
    {
      std::string lPath(fPrefixMCInfoNeededHistos + "v0/MCinformationOnly/" + theBAC);
      //~ fillTH1(fRecTrueV0Histos[kTrue], lPath + "hValidatedPtTrue", lTrueGamma.pt());
      fillV0Histograms(kTrue, lPath, lTrueGamma, nullptr);
      //~ fillV0Histograms(kRec, fPrefixReconstructedInfoHistos + "v0/" + theBAC, theV0, theV0CosinePA);

    }
    
    // reconstructed info of MC validated V0s histos
    {
      std::string lPath(fPrefixMCInfoNeededHistos + "v0/reconstructedInfoOfMCvalidated/" + theBAC);
      fillTH1(fV0ReconstructedInfoMCValidatedHistos, lPath + "hValidatedPtRec", theV0.pt());
      fillTH1(fV0ReconstructedInfoMCValidatedHistos, lPath + "hValidatedRRec", theV0.v0radius());
    }
    
    // v0 resolution histos
    {
      TVector3 lConvPointRec(theV0.x(), theV0.y(), theV0.z());
      TVector3 lConvPointTrue(lTrueGamma.x(), lTrueGamma.y(), lTrueGamma.z());

      std::string lPath(fPrefixMCInfoNeededHistos + "v0/resolutions/" + theBAC);
      fillTH1(fV0ResolutionHistos, lPath + "hPtRes", theV0.pt() - lTrueGamma.pt());
      fillTH1(fV0ResolutionHistos, lPath + "hEtaRes", theV0.eta() - lTrueGamma.eta());
      fillTH1(fV0ResolutionHistos, lPath + "hPhiRes", theV0.phi() - lTrueGamma.phi());
      fillTH1(fV0ResolutionHistos, lPath + "hConvPointRRes", theV0.v0radius() - lConvPointTrue.Perp());
      fillTH1(fV0ResolutionHistos, lPath + "hConvPointAbsoluteDistanceRes", TVector3(lConvPointRec - lConvPointTrue).Mag());
    }
    
    

    
  }
/*
template <typename TV0>
void fillV0Histograms(std::string theHistoPath, const TV0& theV0, float theV0CosinePA)
{
  fillTH1(fV0Histos, theHistoPath + "hPtRec", theV0.pt());
  fillTH1(fV0Histos, theHistoPath + "hEtaRec", theV0.eta());
  fillTH1(fV0Histos, theHistoPath + "hPhiRec", theV0.phi());
  fillTH1(fV0Histos, theHistoPath + "hConvPointR", theV0.v0radius());
  fillTH1(fV0Histos, theHistoPath + "hCosPAngle", theV0CosinePA);
  fillTH2(fV0Histos, theHistoPath + "hArmenteros", theV0.alpha(), theV0.qtarm());
  fillTH2(fV0Histos, theHistoPath + "hPsiPtRec", theV0.psipair(), theV0.pt());
}
* */
  void process(aod::Collisions::iterator  const& theCollision,
               aod::V0Datas               const& theV0s,
               aod::GammaConversionTracks const& theV0Tracks,
               aod::GammaConversionsInfoTrue const& theV0sTrue)
  {
    for (auto& lV0 : theV0s) {
      
      auto lConvPhotonMCTable = theV0sTrue.sliceBy(aod::v0data::v0Id, lV0.v0Id());
      fillTruePhotonHistograms("beforeCuts/",
                               lV0,
                               lConvPhotonMCTable);

      if (!processPhoton(theCollision, lV0, theV0Tracks)) {
        continue;
      }
      
      fillTruePhotonHistograms("afterCuts/",
                               lV0,
                               lConvPhotonMCTable);
    }
  }

  template <typename T>
  std::shared_ptr<T> getTH(mapStringHistPtr const& theMap, std::string const& theName)
  {
    auto lPairIt = theMap.find(theName);
    if (lPairIt == theMap.end()) {
      LOGF(warning, "SFS getTH(): No element with key %s in map.", theName.data());
      return std::shared_ptr<T>(); // SFS todo verify this is correct
    }
    if (!std::holds_alternative<std::shared_ptr<T>>(lPairIt->second)) {
      LOGF(warning, "SFS getTH(): No shared_ptr<T> in map for key %s.", theName.data());
      return std::shared_ptr<T>();
    }
    std::shared_ptr<T> lHisto = std::get<std::shared_ptr<T>>(lPairIt->second);
    if (lHisto == nullptr) {
      LOGF(warning, "SFS getTH(): Found shared_ptr<T> for key %s but it is a nullptr.", theName.data());
      return std::shared_ptr<T>();
    }
    return lHisto;
  }

  void fillTH1(mapStringHistPtr const& theMap, std::string const& theName, float theValue)
  {
    std::shared_ptr<TH1> lHisto = getTH<TH1>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValue);
    }
  }

  void fillTH2(mapStringHistPtr const& theMap, std::string const& theName, float theValueX, float theValueY)
  {
    std::shared_ptr<TH2> lHisto = getTH<TH2>(theMap, theName);
    if (lHisto != nullptr) {
      lHisto->Fill(theValueX, theValueY);
    }
  }

  template <typename TTRACKS>
  void fillTrackHistograms(std::string theHistoPath, TTRACKS const &theV0Tracks)
  {
    auto fillTrackHistogramsI = [&](auto const &theTrack) {
      fillTH1(fTrackHistos, theHistoPath + "hTPCFoundOverFindableCls", theTrack.tpcFoundOverFindableCls());
      fillTH1(fTrackHistos, theHistoPath + "hTPCCrossedRowsOverFindableCls", theTrack.tpcCrossedRowsOverFindableCls());
      fillTH2(fTrackHistos, theHistoPath + "hTPCdEdxSigEl", theTrack.p(), theTrack.tpcNSigmaEl());
      fillTH2(fTrackHistos, theHistoPath + "hTPCdEdxSigPi", theTrack.p(), theTrack.tpcNSigmaPi());
      fillTH2(fTrackHistos, theHistoPath + "hTPCdEdx", theTrack.p(), theTrack.tpcSignal());
    };
    
    for (auto &lTrack : theV0Tracks){
      fillTrackHistogramsI(lTrack);
    }
  }

  //todo: turn this into a function with a switch to fill either true or rec
  template <typename TV0>
  void fillV0Histograms(eRecTrueEnum theRecTrue, std::string theHistoPath, TV0 const &theV0, float const *theV0CosinePA)
  {
    
    fillTH1(fRecTrueV0Histos[theRecTrue], theHistoPath + "hEta" + fRecTrueStrings[theRecTrue], theV0.eta());
    fillTH1(fRecTrueV0Histos[theRecTrue], theHistoPath + "hPhi" + fRecTrueStrings[theRecTrue], theV0.phi());
    fillTH1(fRecTrueV0Histos[theRecTrue], theHistoPath + "hPt"  + fRecTrueStrings[theRecTrue], theV0.pt());
    
    fillTH1(fRecTrueV0Histos[theRecTrue], theHistoPath + "hConvPointR" + fRecTrueStrings[theRecTrue], theV0.v0radius());
    if (theV0CosinePA){
      fillTH1(fRecTrueV0Histos[theRecTrue], theHistoPath + "hCosPAngle"  + fRecTrueStrings[theRecTrue], *theV0CosinePA);
    }
    fillTH2(fRecTrueV0Histos[theRecTrue], theHistoPath + "hArmenteros" + fRecTrueStrings[theRecTrue], theV0.alpha(), theV0.qtarm());
    fillTH2(fRecTrueV0Histos[theRecTrue], theHistoPath + "hPsiPt" + fRecTrueStrings[theRecTrue], theV0.psipair(), theV0.pt());
  }

  template <typename T>
  bool trackPassesCuts(const T& theTrack)
  {
    // single track eta cut
    if (TMath::Abs(theTrack.eta()) > fEtaCut) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTrackEta"));
      return kFALSE;
    }

    // single track pt cut
    if (theTrack.pt() < fSinglePtCut) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTrackPt"));
      return kFALSE;
    }

    if (!(selectionPIDTPC_track(theTrack))) {
      return kFALSE;
    }

    if (theTrack.tpcFoundOverFindableCls() < fMinTPCFoundOverFindableCls) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTPCFoundOverFindableCls"));
      return kFALSE;
    }

    if (theTrack.tpcCrossedRowsOverFindableCls() < fMinTPCCrossedRowsOverFindableCls) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kTPCCrossedRowsOverFindableCls"));
      return kFALSE;
    }

    return kTRUE;
  }

  template <typename T>
  bool passesPhotonCuts(const T& theV0, float theV0CosinePA)
  {
    if (theV0.v0radius() < fMinR || theV0.v0radius() > fMaxR) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kV0Radius"));
      return kFALSE;
    }

    if (!ArmenterosQtCut(theV0.alpha(), theV0.qtarm(), theV0.pt())) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kArmenteros"));
      return kFALSE;
    }

    if (TMath::Abs(theV0.psipair()) > fPsiPairCut) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPsiPair"));
      return kFALSE;
    }

    if (theV0CosinePA < fCosPAngleCut) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kCosinePA"));
      return kFALSE;
    }

    return kTRUE;
  }

  template <typename TV0, typename TTRACKS>
  void fillReconstructedInfoHistograms(std::string theBAC, TV0 const &theV0, TTRACKS const &theV0Tracks, float const &theV0CosinePA)
  {
    fillTrackHistograms(fPrefixReconstructedInfoHistos + "track/" + theBAC, theV0Tracks);
    fillV0Histograms(kRec, fPrefixReconstructedInfoHistos + "v0/" + theBAC, theV0, &theV0CosinePA);

    if (theBAC == "beforeCuts/") {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPhotonIn"));
    } else if (theBAC == "afterCuts/") {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPhotonOut"));
    }
  }

  Bool_t ArmenterosQtCut(Double_t theAlpha, Double_t theQt, Double_t thePt)
  {
    // in AliPhysics this is the cut for if fDo2DQt && fDoQtGammaSelection == 2
    Float_t lQtMaxPtDep = fQtPtMax * thePt;
    if (lQtMaxPtDep > fQtMax) {
      lQtMaxPtDep = fQtMax;
    }
    if (!(TMath::Power(theAlpha / fMaxPhotonAsymmetry, 2) + TMath::Power(theQt / lQtMaxPtDep, 2) < 1)) {
      return kFALSE;
    }
    return kTRUE;
  }

  template <typename T>
  bool selectionPIDTPC_track(const T& theTrack)
  {
    // TPC Electron Line
    if (theTrack.tpcNSigmaEl() < fPIDnSigmaBelowElectronLine || theTrack.tpcNSigmaEl() > fPIDnSigmaAboveElectronLine) {
      fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kElectronPID"));
      return kFALSE;
    }

    // TPC Pion Line
    if (theTrack.p() > fPIDMinPnSigmaAbovePionLine) {
      // low pt Pion rej
      if (theTrack.p() < fPIDMaxPnSigmaAbovePionLine) {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLine) {
          fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPionRejLowMom"));
          return kFALSE;
        }
      }
      // High Pt Pion rej
      else {
        if (theTrack.tpcNSigmaEl() > fPIDnSigmaBelowElectronLine && theTrack.tpcNSigmaEl() < fPIDnSigmaAboveElectronLine && theTrack.tpcNSigmaPi() < fPIDnSigmaAbovePionLineHighP) {
          fillTH1(fRecTrueV0Histos[kRec], fFullNameIsPhotonSelectedHisto, getPhotonCutIndex("kPionRejHighMom"));
          return kFALSE;
        }
      }
    }
    return kTRUE;
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{adaptAnalysisTask<GammaConversionsConsumermc>(cfgc)};
}
