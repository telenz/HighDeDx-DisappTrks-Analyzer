#ifndef HISTOGRAMCLASS_H
#define HISTOGRAMCLASS_H

#include "analyzer.h"
#include "functions.h"
#include "hitInformation.h"
#include "TH1.h"
#include "TVector3.h"

using namespace std;

const double K = 2.529;
const double C = 2.772;

class Hist
{

 public:
  TH1D *htrackPt;
  TH1D *htrackP;
  TH1D *htrackEta;
  //TH1D *htrackd0;
  //TH1D *htrackdz;
  TH1D *htrackNValid;
  TH1D *htrackNLostMid;
  TH1D *htrackNLostInner;
  TH1D *htrackNLostOuter;
  TH1D *htrackIsolation;
  TH1D *htrackCaloIsolation;
  TH1D *htrackASmi;
  TH1D *htrackASmiNP;
  TH1D *htrackDeDxHarm2;
  TH2D *htrackPtDeDxHarm2;
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


 public:Hist(TString histName, outputFile ofile_)
    {
      ofile_.file_->mkdir(histName);
      ofile_.file_->cd(histName);
      htrackPt = new TH1D("htrackPt","htrackPt",100,0,2000);
      htrackP = new TH1D("htrackP","htrackP",100,0,2000);
      htrackEta = new TH1D("htrackEta","htrackEta",100,-5,5);
      //htrackd0 = new TH1D("htrackd0","htrackd0",100,-5,5);
      //htrackdz = new TH1D("htrackdz","htrackdz",100,-50,50);
      htrackNValid = new TH1D("htrackNValid","htrackNValid",40,0,40);
      htrackNLostMid = new TH1D("htrackNLostMid","htrackNLostMid",10,0,10);
      htrackNLostInner = new TH1D("htrackNLostInner","htrackNLostInner",20,0,20);
      htrackNLostOuter = new TH1D("htrackNLostOuter","htrackNLostOuter",20,0,20);
      //hRelativePtError = new TH1D("hRelativePtError" ,"hRelativePtError",300,0,3);
      htrackIsolation = new TH1D("htrackIsolation" ,"htrackIsolation",1000,0,50);
      htrackCaloIsolation = new TH1D("htrackCaloIsolation" ,"htrackCaloIsolation",1000,0,3000);
      htrackDeDxHarm2 = new TH1D("htrackDeDxHarm2" ,"htrackDeDxHarm2",100,0,50);
      htrackASmi = new TH1D("htrackASmi" ,"htrackASmi",1000,-1.1,1.1);
      htrackASmiNP = new TH1D("htrackASmiNP" ,"htrackASmiNP",1000,-1.1,1.1);

      htrackPdgId = new TH1D("htrackPdgId","htrackPdgId",500,0,500);
      htrackgenParticle = new TH1D("htrackgenParticle","htrackgenParticle",1,0,1);
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
      //htrackPdgIdMass = new TH2D("htrackPdgIdMass","htrackPdgIdMass",500,0,500,75,0,1500);
      //htrackPMass = new TH2D("htrackPMass","htrackPMass",150,0,1500,75,0,1500);
      //htrackgenBetaMass = new TH2D("htrackgenBetaMass","htrackgenBetaMass",100,0,1,75,0,1500);
      //htrackgenBetaGammaMass = new TH2D("htrackgenBetaGammaMass","htrackgenBetaGammaMass",30000,0,3000,75,0,1500);

      //hMass = new TH1D("hMass","hMass",75,0,1500);
      //hMass_7Hits = new TH1D("hMass_7Hits","hMass_7Hits",75,0,1500);
      //hMass_3Hits = new TH1D("hMass_3Hits","hMass_3Hits",75,0,1500);
      //hMass_1Hits = new TH1D("hMass_1Hits","hMass_1Hits",75,0,1500);

      hnPFJetsub    = new TH1D("hnPFJetsub","hnPFJetsub",20,0,20);
      hDeltaPhi     = new TH1D("hDeltaPhi","hDeltaPhi",32,0,3.2);
      hDeltaPhiMax  = new TH1D("hDeltaPhiMax","hDeltaPhiMax",32,0,3.2);
      h1stjetpt     = new TH1D("h1stjetpt","h1stjetpt",200,0,2000);
      htrackpt1stjetpt  = new TH2D("htrackpt1stjetpt","htrackpt1stjetpt",100,0,1000,200,0,2000);
      htrackPtDeDxHarm2 = new TH2D("htrackPtDeDxHarm2","htrackPtDeDxHarm2",100,0,1000,100,0,50);
      htrackPtCaloIso   = new TH2D("htrackPtCaloIso","htrackPtCaloIso",100,0,1000,100,0,100);
      htrackPtNLostOuter = new TH2D("htrackPtNLostOuter","htrackPtNLostOuter",100,0,1000,20,0,20);
      hMet          = new TH1D("hMet","hMet",150,0,1500);

      hgenPtChi    = new TH1D("hgenPtChi","hgenPtChi",300,0,3000);
      hgenPChi     = new TH1D("hgenPChi","hgenPChi",300,0,3000);
      hgenEtaChi   = new TH1D("hgenEtaChi","hgenEtaChi",200,-5,5);
      hgenPhiChi   = new TH1D("hgenPhiChi","hgenPhiChi",100,0,3.142);
      hgenBetaChi  = new TH1D("hgenBetaChi","hgenBetaChi",100,0,1);
      
      //chargino plots
      htrackPtoverGenPt = new TH1D("htrackPtoverGenPt","htrackPtoverGenPt",80,0,4);
      htrackEfficiency  = new TH1D("htrackEfficiency","htrackEfficiency",2,0,2);
      
     
    };

  void FillTrackVariables(std::vector<evt::Track_s> inputCollection,double weight)
  {  

    for(unsigned int i=0; i<inputCollection.size(); i++){
      TVector3 trackVector;
      trackVector.SetPtEtaPhi(inputCollection[i].pt,inputCollection[i].eta,inputCollection[i].phi);
      double p = std::sqrt(std::pow(inputCollection[i].pt,2) + std::pow(trackVector.Pz(),2));
      double ASmi   = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], 1, template_strip, template_pixel,1);
      double ASmiNP = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], 1, template_strip, template_strip,0);
      htrackP                  ->Fill(p, weight);
      htrackPt                 ->Fill(inputCollection[i].pt, weight);
      
      htrackEta                ->Fill(inputCollection[i].eta, weight);
      htrackNValid             ->Fill(inputCollection[i].numberOfValidHits, weight);
      htrackNLostMid           ->Fill(inputCollection[i].hitPattern_trackerLayersWithoutMeasurement, weight);
      htrackNLostInner         ->Fill(inputCollection[i].trackerExpectedHitsInner_numberOfLostHits, weight);
      htrackIsolation          ->Fill(inputCollection[i].trackRelIso03, weight);
      htrackCaloIsolation      ->Fill(trackCaloIsolation(&inputCollection[i]), weight);
      htrackNLostOuter         ->Fill(inputCollection[i].trackerExpectedHitsOuter_numberOfHits, weight);
      htrackDeDxHarm2          ->Fill(inputCollection[i].dEdxHarm2, weight);
      htrackASmi               ->Fill(ASmi, weight);
      htrackASmiNP             ->Fill(ASmiNP, weight);
      htrackPtDeDxHarm2        ->Fill(inputCollection[i].pt,inputCollection[i].dEdxHarm2, weight);
      htrackPtCaloIso          ->Fill(inputCollection[i].pt,trackCaloIsolation(&inputCollection[i]), weight);
      htrackPtNLostOuter       ->Fill(inputCollection[i].pt,inputCollection[i].trackerExpectedHitsOuter_numberOfHits, weight);

      //htrackd0->Fill(inputCollection[i].d0, weight);
      //htrackdz->Fill(inputCollection[i].dz, weight);
      //htrackPMass-> Fill(p,mass, weight);
      //double mass =  sqrt((inputCollection[i].dEdxHitsNPHarm2_1000 - C)*pow(p,2)/K, weight);
      //double mass_7Hits =  sqrt((inputCollection[i].dEdxHitsNPHarm2_7 - C)*pow(p,2)/K, weight);
      //double mass_3Hits =  sqrt((inputCollection[i].dEdxHitsNPHarm2_3 - C)*pow(p,2)/K, weight);
      //double mass_1Hits =  sqrt((inputCollection[i].dEdxHitsNPHarm2_1 - C)*pow(p,2)/K, weight);
      //hMass->Fill(mass, weight);
      //hMass_7Hits->Fill(mass_7Hits, weight);
      //hMass_3Hits->Fill(mass_3Hits, weight);
      //hMass_1Hits->Fill(mass_1Hits, weight);


      //if(inputCollection[i].beta*1/sqrt(1.-pow(inputCollection[i].beta,2))>3000){
     
      if(abs(inputCollection[i].pdgId==0))        htrackgenParticle->Fill("unmatched", weight);
      else if(abs(inputCollection[i].pdgId)==1)   htrackgenParticle->Fill("d", weight);
      else if(abs(inputCollection[i].pdgId)==2)   htrackgenParticle->Fill("u", weight);
      else if(abs(inputCollection[i].pdgId)==3)   htrackgenParticle->Fill("s", weight);
      else if(abs(inputCollection[i].pdgId)==4)   htrackgenParticle->Fill("c", weight);
      else if(abs(inputCollection[i].pdgId)==5)   htrackgenParticle->Fill("b", weight);
      else if(abs(inputCollection[i].pdgId)==6)   htrackgenParticle->Fill("t", weight);
      else if(abs(inputCollection[i].pdgId)==11)  htrackgenParticle->Fill("e", weight);
      else if(abs(inputCollection[i].pdgId)==13)  htrackgenParticle->Fill("mu", weight);
      else if(abs(inputCollection[i].pdgId)==15)  htrackgenParticle->Fill("#tau", weight);
      else if(abs(inputCollection[i].pdgId)==21)  htrackgenParticle->Fill("g", weight);
      else if(abs(inputCollection[i].pdgId)==22)  htrackgenParticle->Fill("#gamma", weight);
      else if(abs(inputCollection[i].pdgId)==211) htrackgenParticle->Fill("pi", weight);
      else if(abs(inputCollection[i].pdgId)<1000 &&abs(inputCollection[i].pdgId)!=211 && abs(inputCollection[i].pdgId)>23)   htrackgenParticle->Fill("mesons", weight);
      else if(abs(inputCollection[i].pdgId)>1000 && abs(inputCollection[i].pdgId)<10000)   htrackgenParticle->Fill("baryons", weight);
      else                                        htrackgenParticle->Fill("others", weight);

      if(inputCollection[i].beta<10){
	htrackPdgId->Fill(abs(inputCollection[i].pdgId), weight);
	//htrackPdgIdMass->Fill(abs(inputCollection[i].pdgId),mass weight);
	//htrackgenBetaMass-> Fill(inputCollection[i].beta,mass weight);
	//if(inputCollection[i].beta>0.9999999){
	//  htrackgenBetaGammaMass-> Fill(2500,mass weight);
	//}
	//else{
	//  htrackgenBetaGammaMass-> Fill(inputCollection[i].beta*1/sqrt(1.-pow(inputCollection[i].beta,2)),mass weight);
	//	}
      }
    }
  };

  void FillGenParticleHistograms(std::vector<evt::GenParticle_s> inputCollection, double weight)
  {
    for(unsigned int i=0; i<inputCollection.size(); i++){

      hgenBetaChi       -> Fill(inputCollection[i].p/inputCollection[i].energy, weight);
      hgenPtChi         -> Fill(inputCollection[i].pt, weight);
      hgenPChi          -> Fill(inputCollection[i].p, weight);
      hgenEtaChi        -> Fill(inputCollection[i].eta, weight);
      hgenPhiChi        -> Fill(inputCollection[i].phi, weight);
    }
  }


  void FillCharginoHistograms(std::vector<ChiTrack_s> inputCollection, double weight)
  {

    for(unsigned int i=0; i<inputCollection.size(); i++){
      
      htrackEfficiency ->Fill(inputCollection[i].matched, weight);
      if(inputCollection[i].matched){
	htrackPtoverGenPt ->Fill(inputCollection[i].genpt/inputCollection[i].pt, weight);
      }
      
    }
  }

};


#endif
