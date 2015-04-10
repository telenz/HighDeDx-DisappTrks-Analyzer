#ifndef HISTOGRAMCLASS_H
#define HISTOGRAMCLASS_H

#include "analyzer.h"
#include "functions.h"
#include "hitInformation.h"
#include "declarationsOfClasses.h"
#include "TH1.h"
#include "TVector3.h"
#include <vector>

using namespace std;

const double K   = 2.529; // checked again at 5th december (same for data and MC -> ask Loic about this)
const double C   = 2.772;


Hist::Hist(TString histName, outputFile ofile_)
{

  ofile_.file_->mkdir(histName);
  ofile_.file_->cd(histName);
  tree=new TTree("Variables","a Tree with all relevant variables after selection");
  tree->Branch("weight",&variables.weight);
  tree->Branch("event",&variables.event);
  tree->Branch("run",&variables.run);
  tree->Branch("lumiBlock",&variables.lumiBlock);
  tree->Branch("MET",&variables.met);
  tree->Branch("LeadingJetPt",&variables.LeadingJetPt);
  tree->Branch("trackDeDxASmi",&variables.trackDeDxASmi);
  tree->Branch("trackDeDxASmi_woLastHit",&variables.trackDeDxASmi_woLastHit);
  tree->Branch("trackDeDxHarm2",&variables.trackDeDxHarm2);
  tree->Branch("trackPt",&variables.trackPt);
  tree->Branch("trackP",&variables.trackP);
  tree->Branch("trackGenPt",&variables.trackGenPt);
  tree->Branch("trackGenE",&variables.trackGenE);
  tree->Branch("trackGenEt",&variables.trackGenEt);
  tree->Branch("trackEta",&variables.trackEta);
  tree->Branch("trackNLostOuter",&variables.trackNLostOuter);
  tree->Branch("trackNLostInner",&variables.trackNLostInner);
  tree->Branch("trackNLostMiddle",&variables.trackNLostMiddle);
  tree->Branch("trackNValid",&variables.trackNValid);
  tree->Branch("trackPdgId",&variables.trackPdgId);
  tree->Branch("trackCaloIsolation",&variables.trackCaloIsolation);
  tree->Branch("trackHCALRp5Isolation",&variables.trackHCALRp5Isolation);
  tree->Branch("trackECALRp5Isolation",&variables.trackECALRp5Isolation);
  tree->Branch("trackHCALRp4Isolation",&variables.trackHCALRp4Isolation);
  tree->Branch("trackECALRp4Isolation",&variables.trackECALRp4Isolation);
  tree->Branch("trackHCALRp3Isolation",&variables.trackHCALRp3Isolation);
  tree->Branch("trackECALRp3Isolation",&variables.trackECALRp3Isolation);
  tree->Branch("trackMass",&variables.trackMass);
  tree->Branch("trackIsolation",&variables.trackIsolation);
  tree->Branch("trackEndVertexRho",&variables.trackEndVertexRho);
  tree->Branch("trackChi2",&variables.trackChi2);
  tree->Branch("trackNdof",&variables.trackNdof);
  tree->Branch("trackDeDx1",&variables.trackDeDx1);
  tree->Branch("trackDeDx2",&variables.trackDeDx2);
  tree->Branch("trackDeDx3",&variables.trackDeDx3);
  tree->Branch("trackDeDx4",&variables.trackDeDx4);
  tree->Branch("trackDx1",&variables.trackDx1);
  tree->Branch("trackDx2",&variables.trackDx2);
  tree->Branch("trackDx3",&variables.trackDx3);
  tree->Branch("trackDx4",&variables.trackDx4);
  tree->Branch("trackMeasSize",&variables.trackMeasSize);


  htrackPt            = iniTH1D("htrackPt",100,0,2000);
  htrackPtSmallRange  = iniTH1D("htrackPtSmallRange",40,0,400);
  htrackPtEventCount           = iniTH1D("htrackPtEventCount",100,0,2000);
  htrackPtSmallRangeEventCount = iniTH1D("htrackPtSmallRangeEventCount",40,0,400);
  htrackP             = iniTH1D("htrackP",100,0,2000);
  htrackEta           = iniTH1D("htrackEta",100,-5,5);
  htrackd0            = iniTH1D("htrackd0",100,-1,1);
  htrackdz            = iniTH1D("htrackdz",100,-5,5);
  htrackNValid        = iniTH1D("htrackNValid",40,0,40);
  htrackNValidSmallRange = iniTH1D("htrackNValidSmallRange",20,0,20);
  htrackNValidEventCount           = iniTH1D("htrackNValidEventCount",40,0,40);
  htrackNValidSmallRangeEventCount = iniTH1D("htrackNValidSmallRangeEventCount",20,0,20);
  htrackNLostMid      = iniTH1D("htrackNLostMid",10,0,10);
  htrackNLostInner    = iniTH1D("htrackNLostInner",20,0,20);
  htrackNLostOuter                     = iniTH1D("htrackNLostOuter",15,0,15);
  htrackNLostOuterSmallRange           = iniTH1D("htrackNLostOuterSmallRange",10,0,10);
  htrackNLostOuterEventCount           = iniTH1D("htrackNLostOuterEventCount",15,0,15);
  htrackNLostOuterSmallRangeEventCount = iniTH1D("htrackNLostOuterSmallRangeEventCount",10,0,10);
  htrackIsolation                     = iniTH1D("htrackIsolation",100,0,5);
  htrackIsolationSmallRange           = iniTH1D("htrackIsolationSmallRange",20,0,0.2);
  htrackIsolationEventCount           = iniTH1D("htrackIsolationEventCount",100,0,5);
  htrackIsolationSmallRangeEventCount = iniTH1D("htrackIsolationSmallRangeEventCount",20,0,0.2);
  htrackCaloIsolation                     = iniTH1D("htrackCaloIsolation",100,0,500);
  htrackCaloIsolationSmallRange           = iniTH1D("htrackCaloIsolationSmallRange",20,0,20);
  htrackCaloIsolationEventCount           = iniTH1D("htrackCaloIsolationEventCount",100,0,500);
  htrackCaloIsolationSmallRangeEventCount = iniTH1D("htrackCaloIsolationSmallRangeEventCount",20,0,20);
  htrackMass             = iniTH1D("htrackMass",100,0,1000);   
  htrackMassEventCount   = iniTH1D("htrackMassEventCount",100,0,1000);   
  htrackDeDxHarm2                     = iniTH1D("htrackDeDxHarm2",100,0,40);
  htrackDeDxHarm2SmallRange           = iniTH1D("htrackDeDxHarm2SmallRange",50,0,10);
  htrackDeDxHarm2EventCount           = iniTH1D("htrackDeDxHarm2EventCount",100,0,40);
  htrackDeDxHarm2SmallRangeEventCount = iniTH1D("htrackDeDxHarm2SmallRangeEventCount",50,0,10);
  htrackASmi                     = iniTH1D("htrackASmi",50,-1.1,1.1);
  htrackASmiSmallRange           = iniTH1D("htrackASmiSmallRange",20,0,1);
  htrackASmiEventCount           = iniTH1D("htrackASmiEventCount",50,-1.1,1.1);
  htrackASmiSmallRangeEventCount = iniTH1D("htrackASmiSmallRangeEventCount",20,0,1);
  htrackASmi_3           = iniTH1D("htrackASmi_3",50,-1.1,1.1);
  htrackASmi_7           = iniTH1D("htrackASmi_7",50,-1.1,1.1);
  htrackASmiNP           = iniTH1D("htrackASmiNP",50,-1.1,1.1);
  htrackASmiNPSmallRange = iniTH1D("htrackASmiNPSmallRange",20,0,1);
  htrackASmiNP_3         = iniTH1D("htrackASmiNP_3",50,-1.1,1.1);
  htrackASmiNP_7         = iniTH1D("htrackASmiNP_7",50,-1.1,1.1);
  htrackHighPurity       = iniTH1D("htrackHighPurity",2,0,2);
  htrackPdgId       = iniTH1D("htrackPdgId",500,0,500);
  htrackgenParticle = iniTH1D("htrackgenParticle",1,0,1);
  htrackgenParticle->Fill("unmatched", 0);
  htrackgenParticle->Fill("d", 0);
  htrackgenParticle->Fill("u", 0);
  htrackgenParticle->Fill("s", 0);
  htrackgenParticle->Fill("c", 0);
  htrackgenParticle->Fill("b", 0);
  htrackgenParticle->Fill("t", 0);
  htrackgenParticle->Fill("e", 0);
  htrackgenParticle->Fill("mu", 0);
  htrackgenParticle->Fill("#tau", 0);
  htrackgenParticle->Fill("g", 0);
  htrackgenParticle->Fill("#gamma", 0);
  htrackgenParticle->Fill("pi", 0);
  htrackgenParticle->Fill("mesons", 0);
  htrackgenParticle->Fill("baryons", 0);
  htrackgenParticle->Fill("others", 0);
  htrackgenParticleSmallRange = iniTH1D("htrackgenParticleSmallRange",1,0,1);
  htrackgenParticleSmallRange -> Fill("unmatched", 0);
  htrackgenParticleSmallRange -> Fill("e", 0);
  htrackgenParticleSmallRange -> Fill("mu", 0);
  htrackgenParticleSmallRange -> Fill("pi", 0);
  htrackgenParticleSmallRange -> Fill("others", 0);
  

  hNumberOfTracks = iniTH1D("hNumberOfTracks",10,0,10);
  htrackMT           = iniTH1D("htrackMT",250,0,500);

  hnPFJetsub         = iniTH1D("hnPFJetsub",20,0,20);
  hDeltaPhi          = iniTH1D("hDeltaPhi",32,0,3.2);
  hDeltaPhiMax       = iniTH1D("hDeltaPhiMax",32,0,3.2);
  hDeltaPhiMaxbeforeCut = iniTH1D("hDeltaPhiMaxbeforeCut",32,0,3.2);
  h1stjetpt          = iniTH1D("h1stjetpt",200,0,2000);
  htrackpt1stjetpt              = iniTH2D("htrackpt1stjetpt",100,0,1000,200,0,2000);
  htrackPtDeDxHarm2             = iniTH2D("htrackPtDeDxHarm2",50,0,500,50,0,50);
  htrackPtDeDxHarm2LargeRange   = iniTH2D("htrackPtDeDxHarm2LargeRange",2000,0,2000,50,0,50);
  htrackPtDeDxHarm2SmallBinning = iniTH2D("htrackPtDeDxHarm2SmallBinning",100,0,500,100,0,50);
  htrackPtASmi                  = iniTH2D("htrackPtASmi",50,00,500,20,0,1);
  htrackPtASmiLargeRange        = iniTH2D("htrackPtASmiLargeRange",2000,00,2000,20,0,1);
  htrackPtASmiSmallBinning      = iniTH2D("htrackPtASmiSmallBinning",100,00,500,50,0,1);
  htrackPtCaloIso               = iniTH2D("htrackPtCaloIso",50,0,500,20,0,100);
  htrackPtCaloIsoLargeRange     = iniTH2D("htrackPtCaloIsoLargeRange",2000,0,2000,20,0,100);
  htrackPtCaloIsoSmallBinning   = iniTH2D("htrackPtCaloIsoSmallBinning",100,0,500,50,0,100);
  htrackPtNLostOuter            = iniTH2D("htrackPtNLostOuter",50,0,500,20,0,20);
  htrackCaloIsoASmi             = iniTH2D("htrackCaloIsoASmi",20,0,100,20,0,1);
  htrackCaloIsoASmiLargeRange   = iniTH2D("htrackCaloIsoASmiLargeRange",40,0,200,10,0,1);
  htrackCaloIsoASmiSmallBinning = iniTH2D("htrackCaloIsoASmiSmallBinning",50,0,100,50,0,1);
  
  hMonoCentralPFJet80_PFMETnoMu95            = iniTH1D("hMonoCentralPFJet80_PFMETnoMu95",3,-1,2);
  hMonoCentralPFJet80_PFMETnoMu105           = iniTH1D("hMonoCentralPFJet80_PFMETnoMu105",3,-1,2);
  hMET120_HBHENoiseCleaned                   = iniTH1D("hMET120_HBHENoiseCleaned",3,-1,2);
  hMonoCentralPFJet80_PFMETnoMu95_prescale   = iniTH1D("hMonoCentralPFJet80_PFMETnoMu95_prescale",101,-1,100);
  hMonoCentralPFJet80_PFMETnoMu105_prescale  = iniTH1D("hMonoCentralPFJet80_PFMETnoMu105_prescale",101,-1,100);
  hMET120_HBHENoiseCleaned_prescale          = iniTH1D("hMET120_HBHENoiseCleaned_prescale",101,-1,100);
  hLuminosityBlock                           = iniTH1D("hLuminosityBlock",2500,0,2500);
  hMet                          = iniTH1D("hMet",150,0,1500);

  hgenPtChi              = iniTH1D("hgenPtChi",150,0,1500);
  hgenPChi               = iniTH1D("hgenPChi",300,0,3000);
  hgenEtaChi             = iniTH1D("hgenEtaChi",200,-5,5);
  hgenPhiChi             = iniTH1D("hgenPhiChi",100,0,3.142);
  hgenBetaChi            = iniTH1D("hgenBetaChi",100,0,1);
  hgenBetaTimesGammaChi  = iniTH1D("hgenBetaTimesGammaChi",20,0,200);
      
  //chargino plots
  htrackPtoverGenPt         = iniTH1D("htrackPtoverGenPt",80,0,4);
  htrackDeltaRSimRecoTracks = iniTH1D("htrackDeltaRSimRecoTracks",100,0,1);
  hSimTrackType             = iniTH1D("hSimTrackType",100,0,10000000);
  htrackEfficiency          = iniTH1D("htrackEfficiency",2,0,2);
  hAllTracksZRho            = iniTH2D("hAllTracksZRho",700,0,1400,400,0,800.);
  hFoundTracksZRho          = iniTH2D("hFoundTracksZRho",700,0,1400,400,0,800.);
  htrackProperLifetime      = iniTH1D("htrackProperLifetime",1000,0,1000);
     
};


void Hist::FillTrackVariables(std::vector<evt::Track_s> trkCollection,double weight)
{  

  variables.weight=weight; 

  hNumberOfTracks    ->   Fill(trkCollection.size(), weight);

  if(trkCollection.size()>1){

    htrackPtEventCount                      -> Fill(trkCollection[0].pt, weight);
    htrackPtSmallRangeEventCount            -> Fill(trkCollection[0].pt, weight);
    htrackASmiEventCount                    -> Fill(trkCollection[0].ASmi, weight);
    htrackASmiSmallRangeEventCount          -> Fill(trkCollection[0].ASmi, weight);
    htrackDeDxHarm2EventCount               -> Fill(trkCollection[0].dEdxHarm2, weight);
    htrackDeDxHarm2SmallRangeEventCount     -> Fill(trkCollection[0].dEdxHarm2, weight);
    htrackCaloIsolationEventCount           -> Fill(trackCaloIsolation(&trkCollection[0]), weight);
    htrackCaloIsolationSmallRangeEventCount -> Fill(trackCaloIsolation(&trkCollection[0]), weight);
    htrackNLostOuterEventCount              -> Fill(trkCollection[0].trackerExpectedHitsOuter_numberOfHits, weight);
    htrackNLostOuterSmallRangeEventCount    -> Fill(trkCollection[0].trackerExpectedHitsOuter_numberOfHits, weight);
    htrackIsolationEventCount               -> Fill(trkCollection[0].trackRelIso03, weight);
    htrackIsolationSmallRangeEventCount     -> Fill(trkCollection[0].trackRelIso03, weight);
    htrackNValidEventCount                  -> Fill(trkCollection[0].numberOfValidHits, weight); 
    htrackNValidSmallRangeEventCount        -> Fill(trkCollection[0].numberOfValidHits, weight); 
    double p    = std::sqrt(std::pow(trkCollection[0].pt,2) + std::pow(trkCollection[0].pz,2));
    double mass = sqrt(pow(p,2)/K*(trkCollection[0].dEdxHarm2-C));
    htrackMassEventCount                    -> Fill(mass,weight);

  }

  for(unsigned int i=0; i<trkCollection.size(); i++){
    TVector3 trackVector;
    trackVector.SetPtEtaPhi(trkCollection[i].pt,trkCollection[i].eta,trkCollection[i].phi);

    double p     = std::sqrt(std::pow(trkCollection[i].pt,2) + std::pow(trkCollection[i].pz,2));
    double _dvx  = trkCollection[i].vx - Vertex[0].x;
    double _dvy  = trkCollection[i].vy - Vertex[0].y;
    double d0    = ( - _dvx*trkCollection[i].py + _dvy*trkCollection[i].px )/trkCollection[i].pt;
    double _dvz  = trkCollection[i].vz - Vertex[0].z;
    double dZ    = _dvz - ( _dvx*trkCollection[i].px + _dvy*trkCollection[i].py)/trkCollection[i].pt * (trkCollection[i].pz/trkCollection[i].pt);
    
    htrackP                  ->Fill(p, weight);
    htrackPt                 ->Fill(trkCollection[i].pt, weight);
    htrackPtSmallRange       ->Fill(trkCollection[i].pt, weight);
      
    htrackEta                ->Fill(trkCollection[i].eta, weight);
    htrackNValid             ->Fill(trkCollection[i].numberOfValidHits, weight);
    htrackNValidSmallRange   ->Fill(trkCollection[i].numberOfValidHits, weight);
    htrackNLostMid           ->Fill(trkCollection[i].hitPattern_trackerLayersWithoutMeasurement, weight);
    htrackNLostInner         ->Fill(trkCollection[i].trackerExpectedHitsInner_numberOfLostHits, weight);
    htrackIsolation          ->Fill(trkCollection[i].trackRelIso03, weight);
    htrackIsolationSmallRange->Fill(trkCollection[i].trackRelIso03, weight);
    htrackCaloIsolation           ->Fill(trackCaloIsolation(&trkCollection[i]), weight);
    htrackCaloIsolationSmallRange ->Fill(trackCaloIsolation(&trkCollection[i]), weight);
    htrackNLostOuter              ->Fill(trkCollection[i].trackerExpectedHitsOuter_numberOfHits, weight);
    htrackNLostOuterSmallRange    ->Fill(trkCollection[i].trackerExpectedHitsOuter_numberOfHits, weight);
    double MT = sqrt(2*trkCollection[i].pt*MET_pt*(1.-TMath::Cos(TVector2::Phi_mpi_pi(trkCollection[i].phi-MET_phi))));
    htrackMT                      ->Fill(MT, weight); 

    double mass = sqrt(pow(p,2)/K*(trkCollection[i].dEdxHarm2-C));

    htrackMass               ->Fill(mass,weight);
    htrackDeDxHarm2          ->Fill(trkCollection[i].dEdxHarm2, weight);
    htrackDeDxHarm2SmallRange->Fill(trkCollection[i].dEdxHarm2, weight);
    htrackASmi               ->Fill(trkCollection[i].ASmi, weight);
    htrackASmiSmallRange     ->Fill(trkCollection[i].ASmi, weight);
    htrackASmi_3             ->Fill(trkCollection[i].ASmi_3, weight);
    htrackASmi_7             ->Fill(trkCollection[i].ASmi_7, weight);
    htrackASmiNP             ->Fill(trkCollection[i].ASmiNP, weight);
    htrackASmiNPSmallRange   ->Fill(trkCollection[i].ASmiNP, weight);
    htrackASmiNP_3           ->Fill(trkCollection[i].ASmiNP_3, weight);
    htrackASmiNP_7           ->Fill(trkCollection[i].ASmiNP_7, weight);
    htrackHighPurity              ->Fill(trkCollection[i].trackHighPurity, weight);
    htrackPtDeDxHarm2             ->Fill(trkCollection[i].pt,trkCollection[i].dEdxHarm2, weight);
    htrackPtDeDxHarm2LargeRange   ->Fill(trkCollection[i].pt,trkCollection[i].dEdxHarm2, weight);
    htrackPtDeDxHarm2SmallBinning ->Fill(trkCollection[i].pt,trkCollection[i].dEdxHarm2, weight);
    htrackPtASmi                  ->Fill(trkCollection[i].pt,trkCollection[i].ASmi, weight);
    htrackPtASmiLargeRange        ->Fill(trkCollection[i].pt,trkCollection[i].ASmi, weight);
    htrackPtASmiSmallBinning      ->Fill(trkCollection[i].pt,trkCollection[i].ASmi, weight);
    htrackPtCaloIso               ->Fill(trkCollection[i].pt,trackCaloIsolation(&trkCollection[i]), weight);
    htrackPtCaloIsoLargeRange     ->Fill(trkCollection[i].pt,trackCaloIsolation(&trkCollection[i]), weight);
    htrackPtCaloIsoSmallBinning   ->Fill(trkCollection[i].pt,trackCaloIsolation(&trkCollection[i]), weight);
    htrackCaloIsoASmi             ->Fill(trackCaloIsolation(&trkCollection[i]), trkCollection[i].ASmi, weight);
    htrackCaloIsoASmiLargeRange   ->Fill(trackCaloIsolation(&trkCollection[i]), trkCollection[i].ASmi, weight);
    htrackCaloIsoASmiSmallBinning ->Fill(trackCaloIsolation(&trkCollection[i]), trkCollection[i].ASmi, weight);
    htrackPtNLostOuter            ->Fill(trkCollection[i].pt,trkCollection[i].trackerExpectedHitsOuter_numberOfHits, weight);


    htrackd0                 ->Fill(d0, weight);
    htrackdz                 ->Fill(dZ, weight);
 
    if(abs(trkCollection[i].pdgId==0))        htrackgenParticle->Fill("unmatched", weight);
    else if(abs(trkCollection[i].pdgId)==1)   htrackgenParticle->Fill("d", weight);
    else if(abs(trkCollection[i].pdgId)==2)   htrackgenParticle->Fill("u", weight);
    else if(abs(trkCollection[i].pdgId)==3)   htrackgenParticle->Fill("s", weight);
    else if(abs(trkCollection[i].pdgId)==4)   htrackgenParticle->Fill("c", weight);
    else if(abs(trkCollection[i].pdgId)==5)   htrackgenParticle->Fill("b", weight);
    else if(abs(trkCollection[i].pdgId)==6)   htrackgenParticle->Fill("t", weight);
    else if(abs(trkCollection[i].pdgId)==11)  htrackgenParticle->Fill("e", weight);
    else if(abs(trkCollection[i].pdgId)==13)  htrackgenParticle->Fill("mu", weight);
    else if(abs(trkCollection[i].pdgId)==15)  htrackgenParticle->Fill("#tau", weight);
    else if(abs(trkCollection[i].pdgId)==21)  htrackgenParticle->Fill("g", weight);
    else if(abs(trkCollection[i].pdgId)==22)  htrackgenParticle->Fill("#gamma", weight);
    else if(abs(trkCollection[i].pdgId)==211) htrackgenParticle->Fill("pi", weight);
    else if(abs(trkCollection[i].pdgId)<1000 &&abs(trkCollection[i].pdgId)!=211 && abs(trkCollection[i].pdgId)>23)   htrackgenParticle->Fill("mesons", weight);
    else if(abs(trkCollection[i].pdgId)>1000 && abs(trkCollection[i].pdgId)<10000)   htrackgenParticle->Fill("baryons", weight);
    else                                        htrackgenParticle->Fill("others", weight);

    if(abs(trkCollection[i].pdgId==0))        htrackgenParticleSmallRange->Fill("unmatched", weight);
    else if(abs(trkCollection[i].pdgId)==11)  htrackgenParticleSmallRange->Fill("e", weight);
    else if(abs(trkCollection[i].pdgId)==13)  htrackgenParticleSmallRange->Fill("mu", weight);
    else if(abs(trkCollection[i].pdgId)==211) htrackgenParticleSmallRange->Fill("pi", weight);
    else                                      htrackgenParticleSmallRange->Fill("others", weight);


    if(trkCollection[i].beta<10){
      htrackPdgId->Fill(abs(trkCollection[i].pdgId), weight);
      //htrackPdgIdMass->Fill(abs(trkCollection[i].pdgId),mass weight);
      //htrackgenBetaMass-> Fill(trkCollection[i].beta,mass weight);
      //if(trkCollection[i].beta>0.9999999){
      //  htrackgenBetaGammaMass-> Fill(2500,mass weight);
      //}
      //else{
      //  htrackgenBetaGammaMass-> Fill(trkCollection[i].beta*1/sqrt(1.-pow(trkCollection[i].beta,2)),mass weight);
      //	}
    }

    
    // Fill tree variables
    variables.trackDeDxASmi.push_back(trkCollection[i].ASmi);
    variables.trackDeDxASmi_woLastHit.push_back(trkCollection[i].ASmi_woLastHit);
    variables.trackDeDxHarm2.push_back(trkCollection[i].dEdxHarm2);
    variables.trackPt.push_back(trkCollection[i].pt);
    variables.trackP.push_back(p);
    variables.trackGenPt.push_back(trkCollection[i].genPt);
    variables.trackGenE.push_back(trkCollection[i].genE);
    variables.trackGenEt.push_back(trkCollection[i].genEt);
    variables.trackEta.push_back(trkCollection[i].eta);
    variables.trackNLostOuter.push_back(trkCollection[i].trackerExpectedHitsOuter_numberOfHits);
    variables.trackNLostInner.push_back(trkCollection[i].trackerExpectedHitsInner_numberOfLostHits);
    variables.trackNLostMiddle.push_back(trkCollection[i].hitPattern_trackerLayersWithoutMeasurement);
    variables.trackNValid.push_back(trkCollection[i].numberOfValidHits);
    variables.trackPdgId.push_back(trkCollection[i].pdgId);
    variables.trackCaloIsolation.push_back(trackCaloIsolation(&trkCollection[i]));
    variables.trackHCALRp5Isolation.push_back(trkCollection[i].caloHadDeltaRp5);
    variables.trackECALRp5Isolation.push_back(trkCollection[i].caloEMDeltaRp5);
    variables.trackHCALRp4Isolation.push_back(trkCollection[i].caloHadDeltaRp4);
    variables.trackECALRp4Isolation.push_back(trkCollection[i].caloEMDeltaRp4);
    variables.trackHCALRp3Isolation.push_back(trkCollection[i].caloHadDeltaRp3);
    variables.trackECALRp3Isolation.push_back(trkCollection[i].caloEMDeltaRp3);
    variables.trackMass.push_back(mass);
    variables.trackIsolation.push_back(trkCollection[i].trackRelIso03);
    variables.trackEndVertexRho.push_back(trkCollection[i].simEndVertexRho);
    variables.trackChi2.push_back(trkCollection[i].chi2);
    variables.trackNdof.push_back(trkCollection[i].ndof);
    variables.trackDeDx1.push_back(trkCollection[i].DeDx1);
    variables.trackDeDx2.push_back(trkCollection[i].DeDx2);
    variables.trackDeDx3.push_back(trkCollection[i].DeDx3);
    variables.trackDeDx4.push_back(trkCollection[i].DeDx4);
    variables.trackDx1.push_back(trkCollection[i].Dx1);
    variables.trackDx2.push_back(trkCollection[i].Dx2);
    variables.trackDx3.push_back(trkCollection[i].Dx3);
    variables.trackDx4.push_back(trkCollection[i].Dx4);
    variables.trackMeasSize.push_back(trkCollection[i].MeasSize);

  }

  variables.weight    = weight;
  variables.event     = edmEventHelper_event;
  variables.run       = edmEventHelper_run;
  variables.lumiBlock = edmEventHelper_luminosityBlock;
  
};

void Hist::FillGenParticleHistograms(std::vector<evt::GenParticle_s> genCollection, double weight)
{
  for(unsigned int i=0; i<genCollection.size(); i++){

    double betaAux  = genCollection[i].p/genCollection[i].energy;
    double gammaAux = 1./TMath::Sqrt(1.-pow(betaAux,2));

    hgenBetaChi             -> Fill(betaAux, weight);
    hgenBetaTimesGammaChi   -> Fill(betaAux*gammaAux, weight);
    hgenPtChi               -> Fill(genCollection[i].pt, weight);
    hgenPChi                -> Fill(genCollection[i].p, weight);
    hgenEtaChi              -> Fill(genCollection[i].eta, weight);
    hgenPhiChi              -> Fill(genCollection[i].phi, weight);
  }
}


void Hist::FillCharginoHistograms(std::vector<ChiTrack_s> chitrkCollection, double weight)
{

  for(unsigned int i=0; i<chitrkCollection.size(); i++){

    htrackEfficiency ->Fill(chitrkCollection[i].matched, weight);
    if(chitrkCollection[i].matched){
      htrackPtoverGenPt ->Fill(chitrkCollection[i].genpt/chitrkCollection[i].pt, weight);
    }

    double rho = sqrt(pow(chitrkCollection[i].SimVertexposition_x - Vertex[0].x,2)+pow(chitrkCollection[i].SimVertexposition_y - Vertex[0].y,2));
    double z   = abs(chitrkCollection[i].SimVertexposition_z - Vertex[0].z);

    double distance = sqrt( pow(rho,2) + pow(z,2) );
    double beta  = chitrkCollection[i].genp/chitrkCollection[i].genenergy;
    double gamma = 1./sqrt(1.-pow(chitrkCollection[i].genp/chitrkCollection[i].genenergy,2));
    double properlifetime = distance/(beta*gamma);
    htrackProperLifetime -> Fill(properlifetime,weight);

    if(chitrkCollection[i].SimVertexFound){
      hAllTracksZRho->Fill(z,rho);
      if(chitrkCollection[i].matched) hFoundTracksZRho->Fill(z,rho);
    }
    
  }
}

TH1D* Hist::iniTH1D(TString histoName,int nBins,double low, double high){

  TH1D* histo = new TH1D(histoName,histoName,nBins,low,high);
  histo->Sumw2();

  return histo;

}

TH2D* Hist::iniTH2D(TString histoName,int nBinsX,double lowX, double highX,int nBinsY,double lowY, double highY){

  TH2D* histo = new TH2D(histoName,histoName,nBinsX,lowX,highX,nBinsY,lowY,highY);
  histo->Sumw2();

  return histo;

}

//--------------------------------------------------------------------------------------------------
#endif
