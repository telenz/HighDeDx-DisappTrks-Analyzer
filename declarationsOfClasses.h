#ifndef DECLARATIONSOFCLASSES_H
#define DECLARATIONSOFCLASSES_H
//--------------------------------------------------------------------------------------------------

#include "TTree.h"
#include "TRandom3.h"
#include <vector>
using namespace std;
using namespace evt;

bool isSignalC1N1;
bool isSignal;
TRandom3 *randGenerator;
//----------------------------
typedef struct { 

  double                 weight;
  double                 weightReweighting;
  std::vector<Double_t>  trackDeDxASmi;
  std::vector<Double_t>  trackDeDxASmiNP;
  std::vector<Double_t>  trackDeDxASmi_woLastHit;
  std::vector<Double_t>  trackDeDxHarm2;
  std::vector<Double_t>  trackPt;
  std::vector<Double_t>  trackPtError;
  std::vector<Double_t>  trackP;
  std::vector<Double_t>  trackGenPt;
  std::vector<Double_t>  trackGenE;
  std::vector<Double_t>  trackGenEt;
  std::vector<Double_t>  trackEta;
  std::vector<Int_t>     trackNLostOuter;
  std::vector<Int_t>     trackNLostInner;
  std::vector<Int_t>     trackNLostMiddle;
  std::vector<Int_t>     trackNValid;
  std::vector<Int_t>     trackPdgId;
  std::vector<Int_t>     trackStatus;
  std::vector<Double_t>  trackCaloIsolation;
  std::vector<Double_t>  trackMass;
  std::vector<Double_t>  trackIsolation;
  std::vector<Double_t>  trackEndVertexRho;
  std::vector<Double_t>  trackChi2;
  std::vector<Double_t>  trackNdof;
  std::vector<Double_t>  trackd0;
  std::vector<Double_t>  trackdz;
  UInt_t                 event;
  UInt_t                 run;
  UInt_t                 lumiBlock;
  double                 met;
  double                 LeadingJetPt;
  int                    nJets;


  void clearVectors(){

    trackDeDxASmi.clear();
    trackDeDxASmiNP.clear();
    trackDeDxASmi_woLastHit.clear();
    trackDeDxHarm2.clear();
    trackPt.clear();
    trackPtError.clear();
    trackP.clear();
    trackGenPt.clear();
    trackGenE.clear();
    trackGenEt.clear();
    trackEta.clear();
    trackNLostOuter.clear();
    trackNLostInner.clear();
    trackNLostMiddle.clear();
    trackNValid.clear();
    trackPdgId.clear();
    trackStatus.clear();
    trackCaloIsolation.clear();
    trackMass.clear();
    trackIsolation.clear();
    trackEndVertexRho.clear();
    trackChi2.clear();
    trackNdof.clear();
    trackd0.clear();
    trackdz.clear();
  }

} TreeVariables_t;


class Hist
{


 public:

  TTree *tree;
  TreeVariables_t variables;

  TH1D *hPtOfPions;
  TH1D *htrackPt;
  TH1D *htrackPtSmallRange;
  TH1D *htrackPtSmallRangeCoarseBinning;
  TH1D *htrackP;
  TH1D *htrackEta;
  TH1D *htrackAbsEta;
  TH1D *htrackd0;
  TH1D *htrackdz;
  TH1D *htrackNValid;
  TH1D *htrackNValidSmallRange;
  TH1D *htrackNLostMid;
  TH1D *htrackNLostInner;
  TH1D *htrackNLostOuter;
  TH1D *htrackNLostOuterSmallRange;
  TH1D *htrackIsolation;
  TH1D *htrackIsolationSmallRange;
  TH1D *htrackCaloIsolation;
  TH1D *htrackCaloIsolationSmallRange;
  TH1D *htrackASmi;
  TH1D *htrackASmiSmallRange;
  TH1D *htrackASmiSmallRangeFineBinning;
  TH1D *htrackASmi_3;
  TH1D *htrackASmi_7;
  TH1D *htrackASmiNP;
  TH1D *htrackASmiNPSmallRange;
  TH1D *htrackASmiNP_3;
  TH1D *htrackASmiNP_7;
  TH1D *htrackDeltaEDeltaX;
  TH1D *htrackDeDxHarm2;
  TH1D *htrackDeDxHarm2SmallRange;
  TH1D *htrackHighPurity;
  TH1D *htrackMT;

  TH1D *hNumberOfTracks;
  TH2D *htrackPtDeDxHarm2;
  TH2D *htrackPDeDxHarm2;
  TH2D *htrackPtDeDxHarm2LargeRange;
  TH2D *htrackPtDeDxHarm2SmallBinning;
  TH2D *htrackPtASmi;
  TH2D *htrackPtASmiLargeRange;
  TH2D *htrackPtASmiSmallBinning;
  TH2D *htrackPtCaloIso;
  TH2D *htrackPtCaloIsoLargeRange;
  TH2D *htrackPtCaloIsoSmallBinning;
  TH2D *htrackCaloIsoASmi;
  TH2D *htrackCaloIsoASmiLargeRange;
  TH2D *htrackCaloIsoASmiSmallBinning;
  TH2D *htrackPtNLostOuter;

  

  TH1D *htrackPdgId;
  TH1D* htrackdRminTau;
  TH1D* htrackdRminMuon;
  TH1D* htrackdRminElec;
  TH1D* htrackdRminJet;
  TH1D *htrackgenParticle;
  TH1D *htrackgenParticleSmallRange;

  TH1D *htrackMass;
  TH1D* h1stjetpt;
  TH1D* h1stjetptSmallRange;
  TH2D* htrackpt1stjetpt;
  TH1D *hTriggerResults;
  TH1D *hMonoCentralPFJet80_PFMETnoMu95;
  TH1D *hMonoCentralPFJet80_PFMETnoMu105;
  TH1D *hMET120_HBHENoiseCleaned ;
  TH1D *hMonoCentralPFJet80_PFMETnoMu95_prescale;
  TH1D *hMonoCentralPFJet80_PFMETnoMu105_prescale;
  TH1D *hMET120_HBHENoiseCleaned_prescale ;
  TH1D *hLuminosityBlock;
  TH1D* hMet;
  TH1D* hMetSmallRange;
  TH1D *hnPFJetsub;

  TH1D *hDeltaPhi;
  TH1D *hDeltaPhiMax;
  TH1D *hDeltaPhiMaxbeforeCut;
  TH1D *hDeltaPhiJetMetMinbeforeCut;

  TH1D *hgenPtChi;
  TH1D *hgenPtChiCoarseBinning;
  TH1D *hgenPtChiMiddleBinning;
  TH1D *hgenEtaChi;
  TH1D *hgenPhiChi;
  TH1D *hgenBetaChi;
  TH1D *hgenBetaTimesGammaChi;
  TH1D *hgenPChi;

  // chargino plots
  TH1D *htrackPtoverGenPt;
  TH1D *htrackEfficiency;
  TH1D *htrackDeltaRSimRecoTracks;
  TH1D *hSimTrackType;
  TH1D *htrackProperLifetime;
  TH1D *htrackProperLifetimeReweighted;
  
  TH2D *hAllTracksZRho;
  TH2D *hFoundTracksZRho;

 public:Hist(TString histName, outputFile ofile_);
 public:Hist();

  void FillTrackVariables(std::vector<evt::Track_s> inputCollection,double weight);
  void FillGenParticleHistograms(std::vector<evt::GenParticle_s> inputCollection, double weight);
  void FillCharginoHistograms(std::vector<ChiTrack_s> inputCollection, double weight);
  TH1D *iniTH1D(TString histoName,int nBins,double low, double high);
  TH2D *iniTH2D(TString histoName,int nBinsX,double lowX, double highX,int nBinsY,double lowY, double highY);
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
  std::vector<Jet_s> DTsubleadingJetColl;
  std::vector<MuonPFlow_s> MuonColl;
  std::vector<ElectronPFlow_s> ElectronColl;
    
  struct GenParticle_s leadJetGenParticle;
  TH1D *countsTrackCriteria;
  TH1D *countsEventCuts;
  bool triggerRequirements;
  bool trigger;
  bool trackPreselection;
  bool trackGoodQualitySelection;
  bool trackCandidateSelection;
  bool trackMIPSelection;
  bool qcdSupression;
  bool trackCandidateCutFinal;
  bool onlyChi;
  bool noChi;
  bool DTSelection;

  bool isolatedLeptonCut;

  bool TrackPtRequirement;
  bool NumOfLostOuterRequirement;
  bool CaloIsolationRequirement;
  bool DeDxRequirement;
  
  bool invertTrackPtRequirement;
  bool invertCaloIsolationRequirement;
  bool invertNumOfLostOuterRequirement;
  bool invertDeDxRequirement;

  Hist hist;
  double mass;

 public:
  Event();
  Event(TString histName, outputFile ofile_);
  int Selection();
  std::vector<Track_s> finalTrackCuts(std::vector<Track_s> trackCollection, TH1D* countsTrackCriteria);

};
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
class ABCD
{

 public:
  Event CR1;
  Event CR2;
  Event CR3;
  Event SR;

 public:
 ABCD(TString variables,outputFile ofile_,bool isSignal): CR1("chiTracksCR1"+variables,ofile_), CR2("chiTracksCR2"+variables,ofile_), CR3("chiTracksCR3"+variables,ofile_), SR("chiTracksSR"+variables,ofile_) 
    {

      if(isSignal) CR1.onlyChi   = true;
      CR1.triggerRequirements    = true;
      CR1.trackPreselection      = true;
      CR1.qcdSupression          = true;
      CR1.trackCandidateCutFinal = true;
      CR1.NumOfLostOuterRequirement = false;
      CR1.TrackPtRequirement        = false;
      CR1.CaloIsolationRequirement  = false;
      CR1.DeDxRequirement           = false;

      if(isSignal) CR2.onlyChi   = true;
      CR2.triggerRequirements    = true;
      CR2.trackPreselection      = true;
      CR2.qcdSupression          = true;
      CR2.trackCandidateCutFinal = true;
      CR2.NumOfLostOuterRequirement = false;
      CR2.TrackPtRequirement        = false;
      CR2.CaloIsolationRequirement  = false;
      CR2.DeDxRequirement           = false;

      if(isSignal) CR3.onlyChi   = true;
      CR3.triggerRequirements    = true;
      CR3.trackPreselection      = true;
      CR3.qcdSupression          = true;
      CR3.trackCandidateCutFinal = true;
      CR3.NumOfLostOuterRequirement = false;
      CR3.TrackPtRequirement        = false;
      CR3.CaloIsolationRequirement  = false;
      CR3.DeDxRequirement           = false;

      if(isSignal) SR.onlyChi   = true;
      SR.triggerRequirements    = true;
      SR.trackPreselection      = true;
      SR.qcdSupression          = true;
      SR.trackCandidateCutFinal = true;
      SR.NumOfLostOuterRequirement = false;
      SR.TrackPtRequirement        = false;
      SR.CaloIsolationRequirement  = false;
      SR.DeDxRequirement           = false;

    };
    

  
};
//--------------------------------------------------------------------------------------------------
#endif
