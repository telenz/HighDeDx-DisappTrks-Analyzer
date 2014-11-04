#ifndef DECLARATIONSOFCLASSES_H
#define DECLARATIONSOFCLASSES_H
//--------------------------------------------------------------------------------------------------

#include "analyzer.h"
using namespace std;
using namespace evt;

//----------------------------
class Hist
{

 public:
  TH1D *htrackPt;
  TH1D *htrackP;
  TH1D *htrackEta;
  TH1D *htrackd0;
  TH1D *htrackdz;
  TH1D *htrackNValid;
  TH1D *htrackNLostMid;
  TH1D *htrackNLostInner;
  TH1D *htrackNLostOuter;
  TH1D *htrackIsolation;
  TH1D *htrackCaloIsolation;
  TH1D *htrackCaloIsolationSmallRange;
  TH1D *htrackASmi;
  TH1D *htrackASmi_3;
  TH1D *htrackASmi_7;
  TH1D *htrackASmiNP;
  TH1D *htrackASmiNP_3;
  TH1D *htrackASmiNP_7;
  TH1D *htrackDeDxHarm2;
  TH1D *htrackHighPurity;
  TH2D *htrackPtDeDxHarm2;
  TH2D *htrackPtASmi;
  TH2D *htrackPtCaloIso;
  TH2D *htrackPtNLostOuter;


  //TH1D *hRelativePtError;
  TH1D *htrackPdgId;
  TH1D *htrackgenParticle;
  //TH2D *htrackPdgIdMass;
  //TH2D *htrackPMass;
  //TH2D *htrackgenBetaMass;
  //TH2D *htrackgenBetaGammaMass;

  //TH1D *hMass;
  //TH1D *hMass_1Hits;
  //TH1D *hMass_3Hits;
  //TH1D *hMass_7Hits;

  TH1D* h1stjetpt;
  TH2D* htrackpt1stjetpt;
  TH1D* hMet;
  TH1D *hnPFJetsub;

  TH1D *hDeltaPhi;
  TH1D *hDeltaPhiMax;

  TH1D *hgenPtChi;
  TH1D *hgenEtaChi;
  TH1D *hgenPhiChi;
  TH1D *hgenBetaChi;
  TH1D *hgenPChi;

  // chargino plots
  TH1D *htrackPtoverGenPt;
  TH1D *htrackEfficiency;
  TH1D *htrackDeltaRSimRecoTracks;
  TH1D *hSimTrackType;
  
  TH2D *hAllTracksZRho;
  TH2D *hFoundTracksZRho;

 public:Hist(TString histName, outputFile ofile_);

  void FillTrackVariables(std::vector<evt::Track_s> inputCollection,double weight);
  void FillGenParticleHistograms(std::vector<evt::GenParticle_s> inputCollection, double weight);
  void FillCharginoHistograms(std::vector<ChiTrack_s> inputCollection, double weight);
};

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
class Event
{

 public:
  std::vector<Track_s> TrackColl;
  std::vector<Jet_s> JetColl;
  std::vector<Jet_s> subleadingJetColl;
    
  struct GenParticle_s leadJetGenParticle;
  TH1D *countsTrackCriteria;
  TH1D *countsEventCuts;
  bool preselection;
  bool triggerRequirements;
  bool trackPreselection;
  bool qcdSupression;
  bool trackCandidateCut;
  bool trackCandidateSoftCut;
  bool onlyChi;
  bool noChi;


  bool NumOfLostOuterCut;
  bool CaloIsolationCut;
  bool NumOfValidHitsCut;

  bool invertTrackPtRequirement;
  bool invertNumOfValidHitsRequirement;
  bool invertCaloIsolationRequirement;
  bool invertNumOfLostOuterRequirement;
  bool invertDeDxRequirement;

  Hist hist;
  double mass;

 public:
  Event(TString histName, outputFile ofile_);
  int Selection();
  std::vector<Track_s> finalTrackCuts(std::vector<Track_s> trackCollection, TH1D* countsTrackCriteria);

};
//--------------------------------------------------------------------------------------------------
#endif
