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
//std::vector<Track_s> getCandidateTrackCollection(std::vector<Track_s>& Track, std::vector<Jet_s>& jetColl,  TH1D* countsTrackCriteria, TH1D* countsEventCuts);
//std::vector<Track_s> getCandidateTrackCollection_SoftCuts(std::vector<Track_s>& Track, TH1D* countsTrackCriteria);

Event::Event(TString histName, outputFile ofile_):
  hist(histName, ofile_)
    { 

      countsTrackCriteria = new TH1D("countsTrackCriteria","countsTrackCriteria",1,0,1);
      countsTrackCriteria->SetBit(TH1::kCanRebin);
      countsTrackCriteria->SetStats(0);
      /*
      countsTrackCriteria->Fill("beforeTrackCriteria", 0);   
      countsTrackCriteria->Fill("PtGreater50GeV", 0);                         
      countsTrackCriteria->Fill("EtaLess2p1", 0);
      countsTrackCriteria->Fill("EtaLess1p42Gt1p65", 0);
      countsTrackCriteria->Fill("EtaLess0p15Gt0p35", 0);
      countsTrackCriteria->Fill("EtaLess1p55Gt1p85", 0);
      countsTrackCriteria->Fill("isMatchedDeadEcal", 0);
      countsTrackCriteria->Fill("notWithinECALGap", 0);
      countsTrackCriteria->Fill("isMatchedBadCSC", 0);
      countsTrackCriteria->Fill("d0Less0p2mm", 0);
      countsTrackCriteria->Fill("dZLess5mm", 0);
      countsTrackCriteria->Fill("NOfValidHitsGe7", 0);
      countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", 0);
      countsTrackCriteria->Fill("NOfLostHitsInnerEq0", 0);
      countsTrackCriteria->Fill("TrackIsolationDeltaRLess0p05", 0);
      countsTrackCriteria->Fill("InJetCollectionR0p5", 0);
      countsTrackCriteria->Fill("InTauCollectionR0p15", 0);
      countsTrackCriteria->Fill("InElectronCollectionR0p15", 0);
      countsTrackCriteria->Fill("InMuonCollectionR0p15", 0);
      countsTrackCriteria->Fill("CaloIsolation0p5", 0);
      countsTrackCriteria->Fill("NOfLostHitsOuterGe3", 0);
      countsTrackCriteria->Fill("DeDxHarm2Ge3", 0);
      */
      countsEventCuts = new TH1D("countsEventCuts","countsEventCuts",1,0,1);
      countsEventCuts->SetBit(TH1::kCanRebin);
      countsEventCuts->SetStats(0);
      
      preselection          = false;
      trackCandidateCut     = false;
      trackCandidateSoftCut = false;
      onlyChi               = false;
      noChi                 = false;
      triggerRequirements   = false;
      trackPreselection     = false;
      qcdSupression         = false;

      NumOfLostOuterCut  = false;
      CaloIsolationCut   = false;
      NumOfValidHitsCut  = false;

      invertTrackPtRequirement         = false;
      invertNumOfValidHitsRequirement  = false;
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
      double ASmi = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], 1, template_strip, template_pixel,1);
      if(TrackColl[i].pt>=50 && ASmi>=0.4 /*&& trackCaloIsolation(&TrackColl[i])<10*/) return 0;
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  countsEventCuts->Fill("noCuts", weight);
    

  subleadingJetColl = getSubleadingJetCollection();
  if(preselection)
    {

      //.................................................................................//
      if(triggerRequirements){
	// 1.) Trigger Cut
	  
	/*if(edmEventHelper_isRealData){
	  if(edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5  == 1 || 
	  edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5 == 1 ){}
	  else{  return 0;}
	    
	  if(edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5  == 0 &&
	  edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5 ==1){
	  }
	    
	  if(edmTriggerResultsHelper_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5  == 1 || 
	  edmTriggerResultsHelper_HLT_MET120_HBHENoiseCleaned_v3 == 1 ){}
	  else{  return 0;}
	  }
	*/
	countsEventCuts->Fill("triggerCut", weight);
	     
	// 2.) MET cut
	if(MET_pt<=100.) return 0;
	countsEventCuts->Fill("metCut", weight);
    
	// 3.) Leading Jet Cut
	if(JetColl.size()==0) return 0;
	if(!leadingJetRequirementsFullfilled(&JetColl[0], countsEventCuts)) return 0;
      }
      //.................................................................................//
      if(trackPreselection){
	// 4.)
	TrackColl = trackCandidateCuts(TrackColl,countsTrackCriteria);
	countsEventCuts->Fill("trackCandCut", weight);
	  
	// 5.)
	TrackColl = trackCleaningCuts(TrackColl,countsTrackCriteria);
	countsEventCuts->Fill("trackCleaningCut", weight);

	// 6.)
	TrackColl = trackParticleMatchingCuts(TrackColl,subleadingJetColl,countsTrackCriteria);
	countsEventCuts->Fill("trackParticleMatchingCut", weight);
      }
      //.................................................................................//
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
    }
  //.................................................................................//
  


  // Special Collection
  if(onlyChi)
    {
      TrackColl = findChiInRecoTrackCollection(TrackColl,&hist);
    }
  else if(noChi)
    {
      TrackColl = getFakeTracksInTrackCollection(TrackColl);
    }
        
  // Track Candidate Cut
  if(trackCandidateCut)
    {
      TrackColl = finalTrackCuts(TrackColl,countsTrackCriteria);
      if(TrackColl.size()==0) return 0;
    }
  countsEventCuts->Fill("finalTrackCuts", weight);
    
       
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

  for(unsigned int i=0; i<trackCollection.size(); i++){
    //.................................................................................//
    if(invertTrackPtRequirement){
      if(trackCollection[i].pt>70.)                                             continue;
    }
    else{
      if(trackCollection[i].pt<=70.)                                            continue;
    }
    countsTrackCriteria->Fill("PtGreater50GeV", weight);
    //.................................................................................//
    if(NumOfValidHitsCut){
      if(invertNumOfValidHitsRequirement){
	if(trackCollection[i].numberOfValidHits>=7)                               continue;
      }
      else{
	if(trackCollection[i].numberOfValidHits<7)                                continue;
      }
    }
    countsTrackCriteria->Fill("NOfValidHitsGe7", weight);
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
    double ASmi = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], 1, template_strip, template_pixel,1);
    if(invertDeDxRequirement){
      if(ASmi>=0.1)                                                             continue;
    }
    else{
      if(ASmi<0.1)                                                              continue;
    }
    countsTrackCriteria->Fill("DeDxASmiGe0p4", weight);
    //.................................................................................//
    outputCollection.push_back(trackCollection[i]);
  }


  return outputCollection;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/*
  std::vector<Track_s> getCandidateTrackCollection(std::vector<Track_s>& Track, std::vector<Jet_s>& jetColl, TH1D* countsTrackCriteria, TH1D* countsEventCuts){

  std::vector<Track_s> trackCollection;
  std::vector<Track_s> trackCollectionAux;
  trackCollection.clear();
  trackCollectionAux.clear();


  if(Track.size()==0) return Track;
  //  countsEventCuts->Fill("beforeTrackCriteria", weight);

  trackCollection = Track;
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  countsTrackCriteria->Fill("beforeTrackCriteria", weight);
  if(trackCollection[i].pt<=50.)                                                             continue;
  trackCollectionAux.push_back(trackCollection[i]);
  countsTrackCriteria->Fill("PtGreater50GeV", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("PtGreater50GeV", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(std::abs(trackCollection[i].eta)>2.1)                                                  continue;
  trackCollectionAux.push_back(trackCollection[i]);                                  
  countsTrackCriteria->Fill("EtaLess2p1", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("EtaLess2p1", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(std::abs(trackCollection[i].eta)>1.42 && std::abs(trackCollection[i].eta)<1.65)        continue;
  trackCollectionAux.push_back(trackCollection[i]);                        
  countsTrackCriteria->Fill("EtaLess1p42Gt1p65", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("EtaLess1p42Gt1p65", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(std::abs(trackCollection[i].eta)>0.15 && std::abs(trackCollection[i].eta)<0.35)        continue;
  trackCollectionAux.push_back(trackCollection[i]);                        
  countsTrackCriteria->Fill("EtaLess0p15Gt0p35", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("EtaLess0p15Gt0p35", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(std::abs(trackCollection[i].eta)>1.55 && std::abs(trackCollection[i].eta)<1.85)        continue;
  trackCollectionAux.push_back(trackCollection[i]);               
  countsTrackCriteria->Fill("EtaLess1p55Gt1p85", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("EtaLess1p55Gt1p85", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(getTrkIsMatchedDeadEcal(&trackCollection[i]))                                          continue;
  trackCollectionAux.push_back(trackCollection[i]);               
  countsTrackCriteria->Fill("isMatchedDeadEcal", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("isMatchedDeadEcal", weight);
      
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(isWithinIntermoduleGapsOfECAL(&trackCollection[i]))                                     continue;
  trackCollectionAux.push_back(trackCollection[i]);               
  countsTrackCriteria->Fill("notWithinECALGap", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("notWithinECALGap", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(getTrkIsMatchedBadCSC(&trackCollection[i]))                                            continue;
  trackCollectionAux.push_back(trackCollection[i]);               
  countsTrackCriteria->Fill("isMatchedBadCSC", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("isMatchedBadCSC", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  double _dvx = trackCollection[i].vx - Vertex[0].x;
  double _dvy = trackCollection[i].vy - Vertex[0].y;
  double d0 = ( - _dvx*trackCollection[i].py + _dvy*trackCollection[i].px )/trackCollection[i].pt;
  if(abs(d0)>=0.02)                                                                         continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("d0Less0p2mm", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("d0Less0p2mm", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  double _dvx = trackCollection[i].vx - Vertex[0].x;
  double _dvy = trackCollection[i].vy - Vertex[0].y;
  double _dvz = trackCollection[i].vz - Vertex[0].z;
  double dZ = _dvz - ( _dvx*trackCollection[i].px + _dvy*trackCollection[i].py)/trackCollection[i].pt * (trackCollection[i].pz/trackCollection[i].pt);
  if(abs(dZ)>0.5)                                                                            continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("dZLess5mm", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("dZLess5mm", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){ 
  if(trackCollection[i].numberOfValidHits<7)                                                continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("NOfValidHitsGe7", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("NOfValidHitsGe7", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(trackCollection[i].hitPattern_trackerLayersWithoutMeasurement>0)                       continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("NOfLostHitsMiddleEq0", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(trackCollection[i].trackerExpectedHitsInner_numberOfLostHits>0)                        continue;
  trackCollectionAux.push_back(trackCollection[i]);  
  countsTrackCriteria->Fill("NOfLostHitsInnerEq0", weight);
  }     
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("NOfLostHitsInnerEq", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(trackCollection[i].trackRelIso03>=0.05)                           continue;
    
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("TrackIsolationDeltaRLess0p05", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("TrackIsolationDeltaRLess0p05", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(isTrackReconstructedJet(trackCollection[i], jetColl))                                  continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("InJetCollectionR0p5", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("InJetCollectionR0p5", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(isTrackReconstructedTau(&trackCollection[i]))                                       continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("InTauCollectionR0p15", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("InTauCollectionR0p15", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(isTrackReconstructedElectron(&trackCollection[i]))                                       continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("InElectronCollectionR0p15", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("InElectronCollectionR0p15", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(isTrackReconstructedMuon(&trackCollection[i]))                                           continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("InMuonCollectionR0p15", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("InMuonCollectionR0p15", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(trackCaloIsolation(&trackCollection[i])>10)                                            continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("CaloIsolation0p5", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("CaloIsolation0p5", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits<3)                            continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("NOfLostHitsOuterGe3", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("NOfLostHitsOuterGe3", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  for(unsigned int i=0; i<trackCollection.size(); i++){
  if(trackCollection[i].dEdxHarm2<3)                            continue;
  trackCollectionAux.push_back(trackCollection[i]);       
  countsTrackCriteria->Fill("DeDxHarm2Ge3", weight);
  }
  if(trackCollectionAux.size()==0) return trackCollectionAux;
  countsEventCuts->Fill("DeDxHarm2Ge3", weight);
  trackCollection.clear();
  trackCollection=trackCollectionAux;
  trackCollectionAux.clear();
  //.................................................................................//
  return trackCollection;
  }
*/

//--------------------------------------------------------------------------------------------------
#endif

