#ifndef SELECTION_H
#define SELECTION_H
//-----------------------------------------------------------------------------
#include "histogramClass.h"
#include "functions.h"
#include "triggerFunctions.h"
#include "declarationsOfClasses.h"
#include "DTSelection.h"
#include <iostream>
#include <vector>
using namespace std;
//-----------------------------------------------------------------------------
Event::Event():
hist(){};
Event::Event(TString histName, outputFile ofile_):
hist(histName, ofile_)
{ 

  countsTrackCriteria = new TH1D("countsTrackCriteria","countsTrackCriteria",1,0,1);
  countsTrackCriteria->SetBit(TH1::kCanRebin);
  countsTrackCriteria->SetStats(0);
  countsEventCuts = new TH1D("countsEventCuts","countsEventCuts",1,0,1);
  countsEventCuts->SetBit(TH1::kCanRebin);
  countsEventCuts->SetStats(0);
      
  onlyChi                   = false;
  noChi                     = false;
  triggerRequirements       = false;
  trigger                   = false;
  trackPreselection         = false;
  trackGoodQualitySelection = false;
  trackCandidateSelection   = false;
  trackMIPSelection         = false;
  qcdSupression             = false;
  trackCandidateCutFinal    = false;
  DTSelection               = false;

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

  if(trackCandidateCutFinal){
    countsTrackCriteria->Fill("before analysis selection", 0);
    countsTrackCriteria->Fill("trigger fired",0);
    countsTrackCriteria->Fill("p_{T}^{1.jet} > 110GeV", 0);
    countsTrackCriteria->Fill("MET > 100GeV", 0);
    countsTrackCriteria->Fill("#Delta Phi(Jet,Jet) < 2.5", 0);
    countsTrackCriteria->Fill("#Delta Phi(MET,Jet) < 0.5", 0);
    countsTrackCriteria->Fill("#geq 1 recon. trk", 0);
    countsTrackCriteria->Fill("#geq 1 high purity trk", 0);
    countsTrackCriteria->Fill("#geq 1 trk with N_{lost}^{middle} = 0", 0);
    countsTrackCriteria->Fill("#geq 1 trk with N_{lost}^{inner} = 0", 0);
    countsTrackCriteria->Fill("#geq 1 trk with |d0|<0.2mm", 0);
    countsTrackCriteria->Fill("#geq 1 trk with |dz|<5mm", 0);
    countsTrackCriteria->Fill("#geq 1 trk with |#eta|<2.1", 0);
    countsTrackCriteria->Fill("#geq 1 trk with p_{T}>20GeV", 0);  
    countsTrackCriteria->Fill("#geq 1 trk not matched to muon", 0);
    countsTrackCriteria->Fill("#geq 1 trk not matched to electron", 0);
    countsTrackCriteria->Fill("#geq 1 trk not matched to tau", 0);
    countsTrackCriteria->Fill("#geq 1 trk not matched to jet", 0);
    countsTrackCriteria->Fill("#geq 1 trk not matched to dead ECAL cell", 0);
    countsTrackCriteria->Fill("#geq 1 trk not within intermodule gaps", 0);
    countsTrackCriteria->Fill("#geq 1 trk not with 1.42<|#eta|<1.65", 0);
    countsTrackCriteria->Fill("#geq 1 trk not matched to bad CSC cell", 0);
    countsTrackCriteria->Fill("#geq 1 isolated trk", 0);
    countsTrackCriteria->Fill("#geq 1 trk with E_{calo}<5GeV", 0);
  }

  countsEventCuts      -> LabelsDeflate("X");
  countsTrackCriteria  -> LabelsDeflate("X");

  hist.variables.clearVectors();

  TrackColl.clear();
  JetColl.clear();
  PionsFromDecay.clear();
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
  countsEventCuts->Fill("before analysis selection", weight);
  countsTrackCriteria->Fill("before analysis selection", weight);
  
  
  subleadingJetColl = getSubleadingJetCollection(20.);
  DTsubleadingJetColl = getSubleadingJetCollection(30.);
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
      

      if(isSignal){
	double num = randGenerator->Rndm();
	int fired = 0;
	// fraction of luminosity collected in run A+B is 0.268
	if(num<0.2684672792811) fired = edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v;
	else                    fired = edmEventHelperExtra_emulated_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95;

	if(fired                                                         == 1 ||
	   edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v   == 1)
	  {}
	else{ 
	  return 0; 
	}
      }
      else{
	if(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v  == 1 ||
	   edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v == 1 ||
	   edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v                  == 1
	   )
	  {}
	else{ 
	  if(edmEventHelper_isRealData) cout<<"something wrong with triggers!"<<endl;
	  return 0; 
	}
      }
    }
    countsEventCuts->Fill("Trigger fired", weight);
    countsTrackCriteria->Fill("trigger fired",weight);    
    // 3.) Leading Jet Cut
    if(!leadingJetRequirementsFullfilled(&JetColl[0], countsEventCuts)) return 0;
    countsTrackCriteria->Fill("p_{T}^{1.jet} > 110GeV", weight);
    // 2.) MET cut
    if(MET_pt<=0.) return 0;
    countsEventCuts->Fill("MET > 100GeV", weight);
    countsTrackCriteria->Fill("MET > 100GeV", weight);
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
  countsEventCuts->Fill("#Delta Phi(Jet,Jet) < 2.5", weight);
  countsTrackCriteria->Fill("#Delta Phi(Jet,Jet) < 2.5", weight);

  // DeltaPhi(Jet,MET) cut

  double dPhiMin=100000;
  for(unsigned int i=0; i<subleadingJetColl.size(); i++){
    
    if(i==2) break;  
    double dPhi = std::abs(TVector2::Phi_mpi_pi(subleadingJetColl[i].phi-MET_phi));
    if(dPhi<dPhiMin) dPhiMin = dPhi;
    
  }
  hist.hDeltaPhiJetMetMinbeforeCut->Fill(dPhiMin,weight);
  if(qcdSupression){
    if(isMetInJetDirection(subleadingJetColl,MET_phi)) return 0;
  }
 
  countsEventCuts->Fill("#Delta Phi(MET,Jet) < 0.5", weight);
  countsTrackCriteria->Fill("#Delta Phi(MET,Jet) < 0.5", weight);

  //.................................................................................//
  //%%%%%%%%% Isolated Lepton selection %%%%%%%%%%%%%
  if(isolatedLeptonCut){
    
    //MuonColl     = getTightMuonsInEvent();
    //ElectronColl = getTightElectronsInEvent();
    if(evt::ElectronPFlow.size() + evt::MuonPFlow.size() != 1) return 0;

  }
  //.................................................................................//
  //%%%%%%%%% Track MIP selection %%%%%%%%%%%%%
  if(trackMIPSelection){
    // 4.)
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("> 0 reco. track", weight);
    TrackColl = trackMIPCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackMIPCut", weight);
  }
  //%%%%%%%%% Track Good quality selection %%%%%%%%%%%%%
  if(trackGoodQualitySelection){
    // 4.)
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("> 0 reco. track", weight);
    countsTrackCriteria->Fill("> 0 reco. track", weight);
    TrackColl = trackGoodQuality(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackQualityCuts", weight);
  }
 //%%%%%%%%% Track Preselection %%%%%%%%%%%%%
  if(trackCandidateSelection){
  
    // 4.)
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("> 0 reco. track", weight);
    TrackColl = trackGoodQuality(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("> 0 quality track", weight);

    // 5.)
    TrackColl = trackParticleMatchingCuts(TrackColl,subleadingJetColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("> 0 not SM-like track", weight);
    
    // 6.)
    TrackColl = trackCleaningCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    //countsEventCuts->Fill("trackCleaningCut", weight);
  }
  //%%%%%%%%% Track Preselection %%%%%%%%%%%%%
  if(trackPreselection){
  
    // 4.)
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("> 0 reco. track", weight);
    TrackColl = trackGoodQuality(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackQualityCuts", weight);
    
    // 5.)
    TrackColl = trackParticleMatchingCuts(TrackColl,subleadingJetColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackParticleMatchingCut", weight);

    // 6.)
    TrackColl = trackCleaningCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackCleaningCut", weight);

    // 7.)
    TrackColl = trackAnalysisCuts(TrackColl,countsTrackCriteria);
    if(TrackColl.size()==0) return 0;
    countsEventCuts->Fill("trackAnalysisCuts", weight);


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
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DT selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(DTSelection)
    {

      for(unsigned int i=0; i<DTsubleadingJetColl.size(); i++){
	for(unsigned int j=i+1; j<DTsubleadingJetColl.size(); j++){
	  
	  double dPhi = std::abs(TVector2::Phi_mpi_pi( DTsubleadingJetColl[i].phi- DTsubleadingJetColl[j].phi));
	  if(dPhi>dPhiMax) dPhiMax = dPhi;
	  hist.hDeltaPhi->Fill(dPhi,weight);     
	}
      }
      hist.hDeltaPhiMaxbeforeCut->Fill(dPhiMax,weight);
      if(areTwoJetsBackToBack(DTsubleadingJetColl)) return 0;
      hist.hDeltaPhiMax->Fill(dPhiMax,weight);
      countsEventCuts->Fill("#Delta Phi(Jet,Jet) < 2.5", weight);
      countsTrackCriteria->Fill("#Delta Phi(Jet,Jet) < 2.5", weight);
      // DeltaPhi(Jet,MET) cut
      double dPhiMin=100000;
      for(unsigned int i=0; i<DTsubleadingJetColl.size(); i++){
	if(i==2) break;  
	double dPhi = std::abs(TVector2::Phi_mpi_pi(DTsubleadingJetColl[i].phi-MET_phi));
	if(dPhi<dPhiMin) dPhiMin = dPhi;
    
      }
      hist.hDeltaPhiJetMetMinbeforeCut->Fill(dPhiMin,weight);
      if(isMetInJetDirection(DTsubleadingJetColl,MET_phi)) return 0;
      countsEventCuts->Fill("#Delta Phi(MET,Jet) < 0.5", weight);
      countsTrackCriteria->Fill("#Delta Phi(MET,Jet) < 0.5", weight);

      TrackColl = DTSelectionCuts(TrackColl,DTsubleadingJetColl,countsTrackCriteria);
      if(TrackColl.size()==0) return 0;
    }
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //.................................................................................//
  // Fill delta R variables
  
  for(int i=0; i<TrackColl.size();i++){
    isTrackReconstructedJet(&TrackColl[i], subleadingJetColl);
    isTrackReconstructedMuon(&TrackColl[i]);
    isTrackReconstructedTau(&TrackColl[i]);
    isTrackReconstructedElectron(&TrackColl[i]);

    hist.htrackpt1stjetpt->Fill(TrackColl[i].pt,JetColl[0].pt,weight);  
  }
  //.................................................................................//
  matchTrackToSimTrack(TrackColl);     
  matchTrackToGenParticle(TrackColl);
  
  hist.FillTrackVariables(TrackColl,weight);
  hist.FillCharginoHistograms(ChiTrack,weight);
  hist.hMet->Fill(evt::MET_pt,weight);
  hist.hMetSmallRange->Fill(evt::MET_pt,weight);
  hist.hLuminosityBlock->Fill(evt::edmEventHelper_luminosityBlock, weight);
  if(JetColl.size()!=0){
    hist.h1stjetpt->Fill(JetColl[0].pt,weight);
    hist.h1stjetptSmallRange->Fill(JetColl[0].pt,weight);
    hist.variables.LeadingJetPt = JetColl[0].pt;
  }
  hist.variables.met       = MET_pt;
  hist.variables.nJets     = subleadingJetColl.size();
  hist.hnPFJetsub->Fill(subleadingJetColl.size(),weight);

  //-----------------------------------------------
  /*
  for(unsigned int i; i<PionsFromDecay.size(); i++){
    TLorentzVector v4;
    v4.SetPtEtaPhiE(PionsFromDecay[i].momentum_pt,PionsFromDecay[i].momentum_eta,PionsFromDecay[i].momentum_phi,PionsFromDecay[i].momentum_energy);
  
    hist.hPtOfPions->Fill(PionsFromDecay[i].momentum_pt,weight);
    //hist.hPtOfPions->Fill(v4.P(),weight);
  }
  */
  if(chipmGenParticle.size()>0) hist.FillGenParticleHistograms(chipmGenParticle,weight);
  hist.tree->Fill();

  return 0;

};


std::vector<Track_s> Event::finalTrackCuts(std::vector<Track_s> trackCollection, TH1D* countsTrackCriteria){

  std::vector<Track_s> outputCollection;
  if(trackCollection.size()==0) return trackCollection;

  TrackColl.clear();

  bool firstTrack1  = true;
  bool firstTrack2  = true;
  bool firstTrack3  = true;
  bool firstTrack4  = true;

  for(unsigned int i=0; i<trackCollection.size(); i++){

    //.................................................................................//
    if(CaloIsolationRequirement){
      if(invertCaloIsolationRequirement){
	if(trackCaloIsolation(&trackCollection[i])<=5)                           continue;
      }
      else{
	if(trackCaloIsolation(&trackCollection[i])>5)                            continue;
      }
      if(firstTrack2){
	countsTrackCriteria->Fill("#geq 1 trk with E_{calo}<5GeV", weight);
	firstTrack2 = false;
      }
    }
    //.................................................................................//
    if(TrackPtRequirement){
      if(invertTrackPtRequirement){
	if(trackCollection[i].pt>30.)                                             continue;
      }
      else{
	if(trackCollection[i].pt<=30.)                                            continue;
      }
      if(firstTrack1){
	countsTrackCriteria->Fill("PtGreater30GeV", weight);
	firstTrack1 = false;
      }
    }
    //.................................................................................//
    if(DeDxRequirement){
      if(invertDeDxRequirement){
	if(trackCollection[i].ASmi>=0.05){
	  continue;
	}
      }
      else{
	if(trackCollection[i].ASmi<0.05)                                                              continue;
      }
      if(firstTrack4){
	countsTrackCriteria->Fill("#geq 1 trk with I_{as}>0.05", weight);
	firstTrack4 = false;
      }
    }    
    //.................................................................................//
    if(NumOfLostOuterRequirement){
      if(invertNumOfLostOuterRequirement){
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits>=1)           continue;
      }
      else{
	if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits<1)            continue;
      }
      if(firstTrack3){
	countsTrackCriteria->Fill("NOfLostHitsOuterGe1", weight);
	firstTrack3 = false;
      }
    }    
    //.................................................................................//

    outputCollection.push_back(trackCollection[i]);
  }


  return outputCollection;
}
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//--------------------------------------------------------------------------------------------------
#endif

