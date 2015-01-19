#ifndef DECLARATIONSOFCLASSES_H
#define DECLARATIONSOFCLASSES_H
//--------------------------------------------------------------------------------------------------

#include "analyzer.h"
#include "TTree.h"
#include <vector>
using namespace std;
using namespace evt;

//----------------------------
typedef struct { 

  double                 weight;
  std::vector<Double_t>  trackDeDxASmi; 
  std::vector<Double_t>  trackDeDxHarm2;
  std::vector<Double_t>  trackPt;
  std::vector<Double_t>  trackNLostOuter;
  std::vector<Double_t>  trackNValid;
  std::vector<Double_t>  trackCaloIsolation;
  std::vector<Double_t>  trackMass;
  std::vector<Double_t>  trackIsolation; 


  void clearVectors(){

    trackDeDxASmi.clear();
    trackDeDxHarm2.clear();
    trackPt.clear();
    trackNLostOuter.clear();
    trackNValid.clear();
    trackCaloIsolation.clear();
    trackMass.clear();
    trackIsolation.clear();
  }

} TreeVariables_t;

//TFile f("treeAfterSelection.root","recreate");

class Hist
{


 public:

  TTree *tree;
  TreeVariables_t variables;

  TH1D *htrackPt;
  TH1D *htrackPtSmallRange;
  TH1D *htrackPtEventCount;
  TH1D *htrackPtSmallRangeEventCount;
  TH1D *htrackP;
  TH1D *htrackEta;
  TH1D *htrackd0;
  TH1D *htrackdz;
  TH1D *htrackNValid;
  TH1D *htrackNValidSmallRange;
  TH1D *htrackNValidEventCount;
  TH1D *htrackNValidSmallRangeEventCount;
  TH1D *htrackNLostMid;
  TH1D *htrackNLostInner;
  TH1D *htrackNLostOuter;
  TH1D *htrackNLostOuterSmallRange;
  TH1D *htrackNLostOuterEventCount;
  TH1D *htrackNLostOuterSmallRangeEventCount;
  TH1D *htrackIsolation;
  TH1D *htrackIsolationSmallRange;
  TH1D *htrackIsolationEventCount;
  TH1D *htrackIsolationSmallRangeEventCount;
  TH1D *htrackCaloIsolation;
  TH1D *htrackCaloIsolationSmallRange;
  TH1D *htrackCaloIsolationEventCount;
  TH1D *htrackCaloIsolationSmallRangeEventCount;
  TH1D *htrackASmi;
  TH1D *htrackASmiSmallRange;
  TH1D *htrackASmiEventCount;
  TH1D *htrackASmiSmallRangeEventCount;
  TH1D *htrackASmi_3;
  TH1D *htrackASmi_7;
  TH1D *htrackASmiNP;
  TH1D *htrackASmiNPSmallRange;
  TH1D *htrackASmiNP_3;
  TH1D *htrackASmiNP_7;
  TH1D *htrackDeDxHarm2;
  TH1D *htrackDeDxHarm2SmallRange;
  TH1D *htrackDeDxHarm2EventCount;
  TH1D *htrackDeDxHarm2SmallRangeEventCount;
  TH1D *htrackHighPurity;
  TH1D *hNumberOfTracks;
  TH2D *htrackPtDeDxHarm2;
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
  TH1D *htrackgenParticle;

  TH1D *htrackMass;
  TH1D *htrackMassEventCount;

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
  TH1D *hgenBetaTimesGammaChi;
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
    
  struct GenParticle_s leadJetGenParticle;
  TH1D *countsTrackCriteria;
  TH1D *countsEventCuts;
  bool triggerRequirements;
  bool trackPreselection;
  bool qcdSupression;
  bool trackCandidateCutFinal;
  bool onlyChi;
  bool noChi;

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
