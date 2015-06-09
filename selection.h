#ifndef SELECTION_H
#define SELECTION_H
//-----------------------------------------------------------------------------
#include "histogramClass.h"
#include "functions.h"
#include "triggerFunctions.h"
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
  trigger                = false;
  trackPreselection      = false;
  qcdSupression          = false;
  trackCandidateCutFinal = false;

  isolatedLeptonCut      = false;  

  TrackPtRequirement         = false;
  NumOfLostOuterRequirement  = false;
  CaloIsolationRequirement   = false;
  DeDxRequirement            = false;

  invertTrackPtRequirement         = false;
  invertCaloIsolationRequirement   = false;
  invertNumOfLostOuterRequirement  = false;
  invertDeDxRequirement            = false;
}

int Event::Selection()
{
  countsEventCuts      -> LabelsDeflate("X");
  countsTrackCriteria  -> LabelsDeflate("X");

  hist.variables.clearVectors();

  TrackColl.clear();
  JetColl.clear();
    
  TrackColl=evt::Track;
  JetColl=evt::Jet;
  MuonColl=evt::MuonPFlow;
  ElectronColl=evt::ElectronPFlow;

  double dPhiMax=0;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Blinding!!!!!
  //if(edmEventHelper_isRealData){
  //  for(unsigned int i=0; i<TrackColl.size(); i++){
  //    if(TrackColl[i].pt>=35 && TrackColl[i].ASmi>=0.3 && trackCaloIsolation(&TrackColl[i])<5) return 0;
  //  }
  //}
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
  
  //%%%%%%%%% Vertex Requirements - (are already done on ntuple level) %%%%%%%%%%%%%
  //if(!isGoodVertex())  return 0;
  //.................................................................................//

  //%%%%%%%%% Trigger Requirements %%%%%%%%%%%%%
  if(triggerRequirements){
    // 1.) Trigger Cut

    if(trigger){
      int edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v           = getTriggerResult("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v");
      int edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v          = getTriggerResult("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v");
      int edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v                           = getTriggerResult("edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v");
      int Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v  = getTriggerPrescales("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v");
      int Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v = getTriggerPrescales("edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v");
      int Prescale_edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v                  = getTriggerPrescales("edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v");
      hist.hMonoCentralPFJet80_PFMETnoMu95            -> Fill(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v,weight);
      hist.hMonoCentralPFJet80_PFMETnoMu105           -> Fill(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v,weight);
      hist.hMET120_HBHENoiseCleaned                   -> Fill(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v,weight);
      hist.hMonoCentralPFJet80_PFMETnoMu95_prescale   -> Fill(Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v,weight);
      hist.hMonoCentralPFJet80_PFMETnoMu105_prescale  -> Fill(Prescale_edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v,weight);
      hist.hMET120_HBHENoiseCleaned_prescale          -> Fill(Prescale_edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v,weight);
      
      if(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v  == 1 ||
	   edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v == 1 ||
	   edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v                  == 1)
	  {}
      else{ 
	if(edmEventHelper_isRealData) cout<<"something wrong with triggers!"<<endl;
	return 0; 
      }
    }
    countsEventCuts->Fill("triggerCut_OnlyData", weight);
	     
    // 2.) MET cut
    if(MET_pt<=100.) return 0;
    countsEventCuts->Fill("metCut", weight);
    
    // 3.) Leading Jet Cut
    if(!leadingJetRequirementsFullfilled(&JetColl[0], countsEventCuts)) return 0;
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
  hist.hDeltaPhiMaxbeforeCut->Fill(dPhiMax,weight);
  if(qcdSupression){
    if(areTwoJetsBackToBack(subleadingJetColl)) return 0;
  }
  hist.hDeltaPhiMax->Fill(dPhiMax,weight);
  countsEventCuts->Fill("DeltaPhiCut", weight);

  // DeltaPhi(Jet,MET) cut
  if(qcdSupression)
    {
      if(isMetInJetDirection(subleadingJetColl,MET_phi)) return 0;
    }
  countsEventCuts->Fill("1MetwithDeltaPhiMin2Jetsgt0p5", weight);

  //.................................................................................//
  //%%%%%%%%% Isolated Lepton selection %%%%%%%%%%%%%
  if(isolatedLeptonCut){
    
    //MuonColl     = getTightMuonsInEvent();
    //ElectronColl = getTightElectronsInEvent();
    if(evt::ElectronPFlow.size() + evt::MuonPFlow.size() != 1) return 0;

  }
  //.................................................................................//
  //%%%%%%%%% Track Preselection %%%%%%%%%%%%%
  if(trackPreselection){
    // 4.)
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("NonEmptyTrkColl", weight);
    TrackColl = trackCandidateCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackCandCut", weight);
    
    // 5.)
    TrackColl = trackCleaningCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackCleaningCut", weight);

    // 6.)
    TrackColl = trackParticleMatchingCuts(TrackColl,subleadingJetColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackParticleMatchingCut", weight);
  }

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
  matchTrackToSimTrack(TrackColl);     
  matchTrackToGenParticle(TrackColl);
  
  hist.FillTrackVariables(TrackColl,weight);
  hist.FillCharginoHistograms(ChiTrack,weight);
  hist.hMet->Fill(evt::MET_pt,weight);
  hist.hLuminosityBlock->Fill(evt::edmEventHelper_luminosityBlock, weight);
  if(JetColl.size()!=0){
    hist.h1stjetpt->Fill(JetColl[0].pt,weight);
    hist.variables.LeadingJetPt = JetColl[0].pt;
  }
  hist.variables.met       = MET_pt;
  hist.variables.nJets     = subleadingJetColl.size();
  hist.hnPFJetsub->Fill(subleadingJetColl.size(),weight);

  for(unsigned int i=0; i<TrackColl.size(); i++){
    hist.htrackpt1stjetpt->Fill(TrackColl[i].pt,JetColl[0].pt,weight);  
  }
  //-----------------------------------------------

  if(chipmGenParticle.size()>0) hist.FillGenParticleHistograms(chipmGenParticle,weight);

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Blinding!!!!!
  if(edmEventHelper_isRealData){
    for(unsigned int i=0; i<TrackColl.size(); i++){
      if(TrackColl[i].pt>=35 && TrackColl[i].ASmi>=0.3 && trackCaloIsolation(&TrackColl[i])<5) return 0;
    }
  }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  hist.tree->Fill();

  return 0;

};


std::vector<Track_s> Event::finalTrackCuts(std::vector<Track_s> trackCollection, TH1D* countsTrackCriteria){

  std::vector<Track_s> outputCollection;
  if(trackCollection.size()==0) return trackCollection;

  TrackColl.clear();

  for(unsigned int i=0; i<1; i++){

    //.................................................................................//
    if(TrackPtRequirement){
      if(invertTrackPtRequirement){
	if(trackCollection[i].pt>35.)                                             continue;
      }
      else{
	if(trackCollection[i].pt<=35.)                                            continue;
      }
    }
    countsTrackCriteria->Fill("PtGreater70GeV", weight);
    //.................................................................................//
    if(CaloIsolationRequirement){
      if(invertCaloIsolationRequirement){
	if(trackCaloIsolation(&trackCollection[i])<=5)                           continue;
      }
      else{
	if(trackCaloIsolation(&trackCollection[i])>5)                            continue;
      }
    }
    countsTrackCriteria->Fill("CaloIsolation0p5", weight);
    //.................................................................................//
    if(NumOfLostOuterRequirement){
      if(invertNumOfLostOuterRequirement){
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits>=1)           continue;
      }
      else{
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits<1)            continue;
      }
    }
    countsTrackCriteria->Fill("NOfLostHitsOuterGe1", weight);
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

