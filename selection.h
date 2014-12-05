#ifndef SELECTION_H
#define SELECTION_H
//-----------------------------------------------------------------------------
#include "analyzer.h"
#include "histogramClass.h"
#include "functions.h"
#include "declarationsOfClasses.h"
#include <iostream>
#include <vector>
using namespace std;
//-----------------------------------------------------------------------------

Event::Event(TString histName, outputFile ofile_):
  hist(histName, ofile_)
    { 

      countsTrackCriteria = new TH1D("countsTrackCriteria","countsTrackCriteria",1,0,1);
      countsTrackCriteria->SetBit(TH1::kCanRebin);
      countsTrackCriteria->SetStats(0);
      countsEventCuts = new TH1D("countsEventCuts","countsEventCuts",1,0,1);
      countsEventCuts->SetBit(TH1::kCanRebin);
      countsEventCuts->SetStats(0);
      
      onlyChi                = false;
      noChi                  = false;
      triggerRequirements    = false;
      trackPreselection      = false;
      qcdSupression          = false;
      trackCandidateCutFinal = false;

      TrackPtRequirement = true;
      NumOfLostOuterCut  = true;
      CaloIsolationCut   = true;
      DeDxRequirement    = true;

      invertTrackPtRequirement         = false;
      invertCaloIsolationRequirement   = false;
      invertNumOfLostOuterRequirement  = false;
      invertDeDxRequirement            = false;
    }

int Event::Selection()
{
    
  TrackColl.clear();
  JetColl.clear();
    
  TrackColl=evt::Track;
  JetColl=evt::Jet;

  double dPhiMax=0;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Blinding!!!!!
  if(edmEventHelper_isRealData){
   
    for(unsigned int i=0; i<TrackColl.size(); i++){
      if(TrackColl[i].pt>=50 && TrackColl[i].ASmi>=0.4 /*&& trackCaloIsolation(&TrackColl[i])<10*/) return 0;
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  countsEventCuts->Fill("noCuts", weight);
    

  subleadingJetColl = getSubleadingJetCollection();
  // Special Collection
  if(onlyChi)
    {
      TrackColl = findChiInRecoTrackCollection(TrackColl,&hist);
    }
  else if(noChi)
    {
      TrackColl = getFakeTracksInTrackCollection(TrackColl);
    }
      
  //.................................................................................//
  //%%%%%%%%% Trigger Requirements %%%%%%%%%%%%%
  if(triggerRequirements){
    // 1.) Trigger Cut
    /* 
    if(edmEventHelper_isRealData){
      if(
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v1   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v6   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v7   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v8   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v9   == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v10  == 1 ||

	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v2  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v6  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v7  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v8  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v9  == 1 ||
	 edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v10 == 1 ||

	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v1                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v2                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v3                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v4                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v5                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v6                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v7                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v8                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v9                   == 1 ||
	 edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v10                  == 1 
	 ){}
      else{ return 0; }
    }
    */
    countsEventCuts->Fill("triggerCut_OnlyData", weight);
	     
    // 2.) MET cut
    if(MET_pt<=100.) return 0;
    countsEventCuts->Fill("metCut", weight);
    
    // 3.) Leading Jet Cut
    if(JetColl.size()==0) return 0;
    if(!leadingJetRequirementsFullfilled(&JetColl[0], countsEventCuts)) return 0;
  }
  //.................................................................................//
  //%%%%%%%%% Track Preselection %%%%%%%%%%%%%
  if(trackPreselection){
    // 4.)
    TrackColl = trackCandidateCuts(TrackColl,countsTrackCriteria);
    countsEventCuts->Fill("trackCandCut", weight);
	  
    // 5.)
    TrackColl = trackParticleMatchingCuts(TrackColl,subleadingJetColl,countsTrackCriteria);
    countsEventCuts->Fill("trackParticleMatchingCut", weight);

    // 6.)
    TrackColl = trackCleaningCuts(TrackColl,countsTrackCriteria);
    countsEventCuts->Fill("trackCleaningCut", weight);
  }
  //.................................................................................//
  //%%%%%%%%% QCD supression %%%%%%%%%%%%%
  for(unsigned int i=0; i<subleadingJetColl.size(); i++){
	  
    for(unsigned int j=i+1; j< subleadingJetColl.size(); j++){
	    
      double dPhi = std::abs(TVector2::Phi_mpi_pi( subleadingJetColl[i].phi- subleadingJetColl[j].phi));
      if(dPhi>dPhiMax) dPhiMax = dPhi;
      hist.hDeltaPhi->Fill(dPhi,weight);     
    }
  }
  if(qcdSupression){
    if(areTwoJetsBackToBack(subleadingJetColl)) return 0;
  }
  countsEventCuts->Fill("DeltaPhiCut", weight);
	
  // DeltaPhi(Jet,MET) cut
  if(qcdSupression)
    {
      if(isMetInJetDirection(subleadingJetColl,MET_phi)) return 0;
    }
  countsEventCuts->Fill("1MetwithDeltaPhiMin2Jetsgt0p5", weight);
   
  //.................................................................................//
  //%%%%%%%%% Final track cuts  BEGIN %%%%%%%%%%%%%
  // Final Track Cuts
  if(trackCandidateCutFinal)
    {
      TrackColl = finalTrackCuts(TrackColl,countsTrackCriteria);
      if(TrackColl.size()==0) return 0;
    }
  countsEventCuts->Fill("finalTrackCuts", weight);
  //%%%%%%%%% Final track cuts  END %%%%%%%%%%%%%%%
  //.................................................................................//
       
  matchTrackToGenParticle(TrackColl);
  hist.FillTrackVariables(TrackColl,weight);
  hist.FillCharginoHistograms(ChiTrack,weight);
  hist.hMet->Fill(evt::MET_pt,weight);
  if(JetColl.size()!=0) hist.h1stjetpt->Fill(JetColl[0].pt,weight);
  if(dPhiMax>0.) hist.hDeltaPhiMax->Fill(dPhiMax,weight);
  hist.hnPFJetsub->Fill(subleadingJetColl.size(),weight);

  for(unsigned int i=0; i<TrackColl.size(); i++){
    hist.htrackpt1stjetpt->Fill(TrackColl[i].pt,JetColl[0].pt,weight);  
  }
  //-----------------------------------------------

  if(chipmGenParticle.size()>0) hist.FillGenParticleHistograms(chipmGenParticle,weight);


  return 0;

};


std::vector<Track_s> Event::finalTrackCuts(std::vector<Track_s> trackCollection, TH1D* countsTrackCriteria){

  std::vector<Track_s> outputCollection;
  if(trackCollection.size()==0) return trackCollection;

  TrackColl.clear();

  for(unsigned int i=0; i<trackCollection.size(); i++){
    //.................................................................................//
    if(TrackPtRequirement){
      if(invertTrackPtRequirement){
	if(trackCollection[i].pt>70.)                                             continue;
      }
      else{
	if(trackCollection[i].pt<=70.)                                            continue;
      }
    }
    countsTrackCriteria->Fill("PtGreater70GeV", weight);
    //.................................................................................//
    if(CaloIsolationCut){
      if(invertCaloIsolationRequirement){
	if(trackCaloIsolation(&trackCollection[i])<=10)                           continue;
      }
      else{
	if(trackCaloIsolation(&trackCollection[i])>10)                            continue;
      }
    }
    countsTrackCriteria->Fill("CaloIsolation0p5", weight);
    //.................................................................................//
    if(NumOfLostOuterCut){
      if(invertNumOfLostOuterRequirement){
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits>=3)           continue;
      }
      else{
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits<3)            continue;
      }
    }
    countsTrackCriteria->Fill("NOfLostHitsOuterGe3", weight);
    //.................................................................................//
    if(DeDxRequirement){
      if(invertDeDxRequirement){
	if(trackCollection[i].ASmi>=0.2){
	  continue;
	}
      }
      else{
	if(trackCollection[i].ASmi<0.2)                                                              continue;
      }
    }
    countsTrackCriteria->Fill("DeDxASmiGe0p2", weight);
    //.................................................................................//

    outputCollection.push_back(trackCollection[i]);
  }


  return outputCollection;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//--------------------------------------------------------------------------------------------------
#endif

