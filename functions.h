#ifndef FUNCTIONS_H
#define FUNCTIONS_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "chiClass.h"
#include "declarationsOfClasses.h"
#include "TH3.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------

struct evt::GenParticle_s chipGenParticle;
struct evt::GenParticle_s chimGenParticle;

std::vector<evt::GenParticle_s> chipmGenParticle;


std::vector<double> etaCSC, phiCSC, etaEcal, phiEcal;  

bool zeroChip = true;
bool zeroChim = true;
bool zeroJet  = true;

int nChi0 = 0;


int mismatchedGenChiToTrack = 0;
int matchedGenChiToTrack    = 0;

int nChiInSimTrack          = 0;
int nChiInSimVertex         = 0;

//--------------------------------------------------------------------------------------------------

double trackCaloIsolation(struct Track_s *track){


  double caloTowerIso05RhoCorr = 
    std::max(0.,track->caloHadDeltaRp5 + track->caloEMDeltaRp5 - sdouble_value*TMath::Pi()*0.5*0.5); 

  return caloTowerIso05RhoCorr;
}

//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedTau(struct Track_s* track){


  for(unsigned int i=0; i<Tau.size(); i++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0; 

    if(Tau[i].pt<20)                                    continue;
    if(abs(Tau[i].eta)>2.3)                             continue;
    if(Tau[i].byLooseCombinedIsolationDeltaBetaCorr<=0) continue;
    if(Tau[i].decayModeFinding<=0)                      continue;
    if(Tau[i].againstElectronLoose<=0)                  continue;
    if(Tau[i].againstMuonTight<=0)                      continue;

    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Tau[i].phi));
    dEta = std::abs(track->eta - Tau[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedElectron(struct Track_s* track){


  for(unsigned int i=0; i<Electron.size(); i++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0; 

    if(Electron[i].mvaNonTrigV0<0) continue;
    if(Electron[i].pt<10)          continue;
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Electron[i].phi));
    dEta = std::abs(track->eta - Electron[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedMuon(struct Track_s* track){

  for(unsigned int i=0; i<Muon.size(); i++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0;

    if(Muon[i].pt<10.)   continue;
    dPhi = std::abs(TVector2::Phi_mpi_pi(track->phi-Muon[i].phi));
    dEta = std::abs(track->eta - Muon[i].eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta); 

    if(dR<0.15) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isTrackReconstructedJet(struct evt::Track_s track, std::vector<evt::Jet_s>& jetColl){

  for(unsigned int i=0; i<jetColl.size(); i++){
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(track.phi-jetColl[i].phi));
    double dEta = std::abs(track.eta - jetColl[i].eta);
    double dR   = std::sqrt(pow(dPhi,2) + pow(dEta,2)); 

    if(dR<0.5) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
/*****************************************
  Find chargino in GenParticle collection 
*****************************************/
void findChiInGenParticleCollection(){

  zeroChip = true;
  zeroChim = true;

  chipmGenParticle.clear();
  ChiTrack.clear();
  ChiTrack.resize(2);

  int idx=0;

  for(int i=0; i<evt::nGenParticle; i++){

    if(abs(evt::GenParticle[i].pdgId)==1000024){

      if(evt::GenParticle[i].pdgId>0 && zeroChip){

	chipGenParticle = evt::GenParticle[i];
	chipmGenParticle.push_back(evt::GenParticle[i]);
	fillChiTrackWithGenParticleVariables(&ChiTrack[idx], &evt::GenParticle[i]);
	idx+=1;
	zeroChip = false;
      }
      else if(evt::GenParticle[i].pdgId<0 && zeroChim){

	chimGenParticle = evt::GenParticle[i];
	chipmGenParticle.push_back(evt::GenParticle[i]);
	fillChiTrackWithGenParticleVariables(&ChiTrack[idx], &evt::GenParticle[i]);
	idx+=1;
	zeroChim = false;
      }	      
    }
    if(!zeroChip && !zeroChim) break;
  }

  //if(zeroChip || zeroChim) cout<<"To few charginos in GenParticle collection!"<<endl;
}
//--------------------------------------------------------------------------------------------------
/*
// Get only the Chi from the full Track Collection
std::vector<Track_s> findChiInRecoTrackCollection(std::vector<Track_s>& inputCollection){

std::vector<Track_s> chiTrackCollection;
chiTrackCollection.clear();

double dPhi = 0;
double dEta = 0;
double dR   = 0;
int idxMin  = -100;
 
for(unsigned int j=0; j<ChiTrack.size();j++){

double dRmin=10000.;

for(unsigned int i=0; i<inputCollection.size(); i++){
dPhi = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - ChiTrack[j].genphi));
dEta = std::abs(inputCollection[i].eta - ChiTrack[j].geneta);
dR   = std::sqrt( dPhi*dPhi + dEta*dEta );

if(dR<dRmin){
dRmin=dR;
idxMin=i;
}
}
    
if(dRmin<0.001){
chiTrackCollection.push_back(inputCollection[idxMin]);
fillChiTrackWithRecoTrackVariables(&ChiTrack[j], &inputCollection[idxMin]);
ChiTrack[j].matched = true;
matchedGenChiToTrack += 1;
}


}

return chiTrackCollection;
}
*/
//--------------------------------------------------------------------------------------------------
// Get only the Chi from the full Track Collection
std::vector<Track_s> findChiInRecoTrackCollection(std::vector<Track_s>& trkCollection, Hist* hist){

  std::vector<Track_s> chiTrackCollection;
  chiTrackCollection.clear();

  int idxMin  = -100;
 
  for(unsigned int j=0; j<ChiTrack.size();j++){

    double dPhi = 0;
    double dEta = 0;
    double dR   = 0;
    double dRmin=10000.;

    for(unsigned int i=0; i<trkCollection.size(); i++){
      dPhi = std::abs(TVector2::Phi_mpi_pi(trkCollection[i].phi - ChiTrack[j].genphi));
      dEta = std::abs(trkCollection[i].eta - ChiTrack[j].geneta);
      dR   = std::sqrt( dPhi*dPhi + dEta*dEta );

      if(dR<dRmin){
	dRmin=dR;
	idxMin=i;
      }
    }
    
    if(dRmin<0.01){
      chiTrackCollection.push_back(trkCollection[idxMin]);
      fillChiTrackWithRecoTrackVariables(&ChiTrack[j], &trkCollection[idxMin]);
      ChiTrack[j].matched = true;
      matchedGenChiToTrack += 1;
    }
    hist->htrackDeltaRSimRecoTracks -> Fill(dRmin,weight);
    hist->hSimTrackType             -> Fill(std::abs(ChiTrack[j].SimTracktype),weight);

  }

  return chiTrackCollection;
}

//--------------------------------------------------------------------------------------------------
void findChiInSimTrackCollection(){

  for(unsigned int j=0; j<ChiTrack.size();j++){

    for(int i=0; i<nSimTrack; i++){

      if(SimTrack[i].type==ChiTrack[j].genpdgId){

	fillChiTrackWithSimTrackVariables(&ChiTrack[j], &SimTrack[i]);
	
	nChiInSimTrack += 1;

      }
    }
  }
}
//--------------------------------------------------------------------------------------------------
void findChiDecayVertex(){

  for(unsigned int j=0; j<ChiTrack.size();j++){

    ChiTrack[j].SimVertexFound = false;

    for(int i=0; i<nSimVertex; i++){

      bool decayed = false;

      if((unsigned int) SimVertex[i].parentIndex==ChiTrack[j].SimTracktrackId){


	for(int l=0; l<nSimTrack; l++){

	  if(std::abs(SimTrack[l].type)==1000022){
	    if((unsigned int) SimTrack[l].vertIndex == SimVertex[i].vertexId){
	      decayed = true;
	      break;
	    }
	  }
	}

	if(decayed){
	  fillChiTrackWithSimVertexVariables(&ChiTrack[j], &SimVertex[i]);
	  nChiInSimVertex += 1;
	  ChiTrack[j].SimVertexFound = true;
	}
      }
    }
  }
}
//--------------------------------------------------------------------------------------------------
   
struct evt::GenParticle_s  findLeadingJetInGenParticleCollectionWithPtGt30(){

  struct evt::GenParticle_s leadingJetGenParticle;
  zeroJet = true;
  for(int i=0; i<evt::nGenParticle; i++){

    if(abs(evt::GenParticle[i].pdgId)<=6){
      if(evt::GenParticle[i].pt<30) continue;
      leadingJetGenParticle = evt::GenParticle[i];
      zeroJet = false;
      break;
    }
	       
  }
  return leadingJetGenParticle;
}

//--------------------------------------------------------------------------------------------------
bool areTwoJetsBackToBack(std::vector<evt::Jet_s>& jetColl){

  for(unsigned int i=0; i<jetColl.size(); i++){
    for(unsigned int j=i+1; j<jetColl.size(); j++){

      double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-jetColl[j].phi));
      if(abs(dPhi)>=2.5) return true;     
      
    }
  }

  return false;

}

//--------------------------------------------------------------------------------------------------
bool isMetInJetDirection(std::vector<evt::Jet_s> jetColl, double metPhi){

  for(unsigned int i=0; i<jetColl.size(); i++){

    if(i>1) break;
    
    double dPhi = std::abs(TVector2::Phi_mpi_pi(jetColl[i].phi-metPhi));
    if(dPhi<=0.5) return true;     
    
  }

  return false;

}


//--------------------------------------------------------------------------------------------------
bool leadingJetRequirementsFullfilled(struct evt::Jet_s* leadingJet, TH1D* countsEventCuts){

  // jetId: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
  if(leadingJet==0)                                return false;
  if(leadingJet->pt<=110.)                         return false;
  countsEventCuts->Fill("leadingJetPtGt110GeV", evt::weight);
  if(std::abs(leadingJet->eta)>=2.4)               return false;
  countsEventCuts->Fill("absLeadJetEtaLt2p4", evt::weight);
  if(leadingJet->chargedHadronEnergyFraction<=0.2) return false;
  countsEventCuts->Fill("CHEFgt0p2", evt::weight);
  if(leadingJet->chargedEmEnergyFraction>=0.5)     return false;
  countsEventCuts->Fill("CHEmEFle0p5", evt::weight);
  if(leadingJet->neutralHadronEnergyFraction>=0.7) return false;
  countsEventCuts->Fill("NHEFle0p7", evt::weight);
  if(leadingJet->neutralEmEnergyFraction>=0.7)     return false;
  countsEventCuts->Fill("NEmEFle0p7", evt::weight);

  return true;
}
//--------------------------------------------------------------------------------------------------

std::vector<evt::Jet_s>  getSubleadingJetCollection(){

  std::vector<evt::Jet_s> jetCollection;
  jetCollection.clear();
  for(unsigned int i=0; i<evt::Jet.size(); i++){

    //bool isCharginoCandidate = false;

    if(evt::Jet[i].pt<=20.)                          continue;
    if(std::abs(evt::Jet[i].eta)>=4.5)               continue;
    //if(evt::Jet[i].neutralHadronEnergyFraction>=0.7) continue;
    //if(evt::Jet[i].chargedEmEnergyFraction>=0.5)     continue;

    jetCollection.push_back(evt::Jet[i]);
  }

  return jetCollection;
}
//--------------------------------------------------------------------------------------------------
  
bool isGoodVertex(){
    
  if(evt::Vertex[0].z > 24.)                                                     return false;
  if(std::sqrt(std::pow(evt::Vertex[0].x,2) + std::pow(evt::Vertex[0].y,2)) > 2) return false;
  if(evt::Vertex[0].ndof < 4)                                                    return false;
  return true;
}
//--------------------------------------------------------------------------------------------------
bool getTrkIsMatchedDeadEcal(struct evt::Track_s *track){

  for(unsigned int i=0; i<etaEcal.size(); i++){

    double dPhi  = 0;
    double dEta  = 0;
    double dR    = 0;  
        
    dPhi = std::abs(TVector2::Phi_mpi_pi(phiEcal[i] - track->phi));
    dEta = std::abs(etaEcal[i] - track->eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.05) return true;
  }  
  return false;
}
//--------------------------------------------------------------------------------------------------
bool getTrkIsMatchedBadCSC(struct evt::Track_s *track){

  for(unsigned int i=0; i<etaCSC.size(); i++){

    double dPhi  = 0;
    double dEta  = 0;
    double dR    = 0;  
          
    dPhi = std::abs(TVector2::Phi_mpi_pi(phiCSC[i] - track->phi));
    dEta = std::abs(etaCSC[i] - track->eta);
    dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

    if(dR<0.25) return true;
  }

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isWithinIntermoduleGapsOfECAL(struct evt::Track_s *track){

  if(track->eta<-1.14018   && track->eta>-1.1439)     return true;
  if(track->eta<-0.791884  && track->eta>-0.796051)   return true;
  if(track->eta<-0.44356   && track->eta>-0.447911)   return true;
  if(track->eta<0.00238527 && track->eta>-0.00330793) return true;
  if(track->eta<0.446183   && track->eta>0.441949)    return true;
  if(track->eta<0.793955   && track->eta>0.789963)    return true;
  if(track->eta<1.14164    && track->eta>1.13812)     return true;
  
  return false;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackCuts(std::vector<evt::Track_s> inputColl, bool keepCriteria)
{

  std::vector<evt::Track_s> outputColl;

  for(unsigned int i=0; i<inputColl.size(); i++){
    if(!keepCriteria) continue;
    outputColl.push_back(inputColl[i]);
  }

  return outputColl;

}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackCandidateCuts(std::vector<evt::Track_s> trackCollection, TH1D* countsTrackCriteria)
{

  std::vector<evt::Track_s> outputColl;
  bool firstTrack1 = true;
  bool firstTrack2 = true;
  bool firstTrack3 = true;
  bool firstTrack4 = true;

  for(unsigned int i=0; i<trackCollection.size(); i++){
    if(firstTrack1){
      countsTrackCriteria->Fill("beforeTrackCriteria", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    if(trackCollection[i].pt<=10.)                                            continue;
    if(firstTrack2){
      countsTrackCriteria->Fill("PtGreater10GeV", weight);  
      firstTrack2 = false;
    }
    //.................................................................................//
    if(std::abs(trackCollection[i].eta)>2.4)                                  continue;
    if(firstTrack3){
      countsTrackCriteria->Fill("EtaLess2p4", weight);
      firstTrack3 = false;
    }
    //.................................................................................//
    if(!trackCollection[i].trackHighPurity)                                   continue;
    if(firstTrack4){
      countsTrackCriteria->Fill("highPurity", weight);
      firstTrack4 = false;
    }

    //.................................................................................//
    outputColl.push_back(trackCollection[i]);
  }
  return outputColl;
  
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackCleaningCuts(std::vector<evt::Track_s> trackCollection, TH1D* countsTrackCriteria)
{

  std::vector<evt::Track_s> outputColl;

  bool firstTrack1  = true;
  bool firstTrack2  = true;
  bool firstTrack3  = true;
  bool firstTrack4  = true;
  bool firstTrack5  = true;
  bool firstTrack6  = true;
  bool firstTrack7  = true;
  bool firstTrack8  = true;
  bool firstTrack9  = true;
  bool firstTrack10 = true;
  bool firstTrack11 = true;
  bool firstTrack12 = true;


  for(unsigned int i=0; i<trackCollection.size(); i++){

    //.................................................................................//
    //if(std::abs(trackCollection[i].eta)>1.42 && std::abs(trackCollection[i].eta)<1.65)        continue;
    if(firstTrack1){
      countsTrackCriteria->Fill("EtaLess1p42Gt1p65", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    //if(std::abs(trackCollection[i].eta)>0.15 && std::abs(trackCollection[i].eta)<0.35)        continue;
    if(firstTrack1){
      countsTrackCriteria->Fill("EtaLess0p15Gt0p35", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    //if(std::abs(trackCollection[i].eta)>1.55 && std::abs(trackCollection[i].eta)<1.85)        continue;
    if(firstTrack2){
      countsTrackCriteria->Fill("EtaLess1p55Gt1p85", weight);
      firstTrack2 = false;
    }
    //.................................................................................//
    if(getTrkIsMatchedDeadEcal(&trackCollection[i]))                                          continue;
    if(firstTrack3){
      countsTrackCriteria->Fill("isMatchedDeadEcal", weight);
      firstTrack3 = false;
    }
    //.................................................................................//
    if(isWithinIntermoduleGapsOfECAL(&trackCollection[i]))                                    continue;
    if(firstTrack4){
      countsTrackCriteria->Fill("notWithinECALGap", weight);
      firstTrack4 = false;
    }
    //.................................................................................//
    if(getTrkIsMatchedBadCSC(&trackCollection[i]))                                            continue;
    if(firstTrack5){
      countsTrackCriteria->Fill("isMatchedBadCSC", weight);
      firstTrack5 = false;
    }
    //.................................................................................//
    double _dvx = trackCollection[i].vx - Vertex[0].x;
    double _dvy = trackCollection[i].vy - Vertex[0].y;
    double d0 = abs( - _dvx*trackCollection[i].py + _dvy*trackCollection[i].px)/trackCollection[i].pt;
    if(abs(d0)>0.02)                                                                          continue;
    if(firstTrack6){
      countsTrackCriteria->Fill("d0Less0p2mm", weight);
      firstTrack6 = false;
    }
    //.................................................................................//
    double _dvz = trackCollection[i].vz - Vertex[0].z;
    double dZ = _dvz - ( _dvx*trackCollection[i].px + _dvy*trackCollection[i].py)/trackCollection[i].pt * (trackCollection[i].pz/trackCollection[i].pt);
    if(abs(dZ)>0.5)                                                                           continue;
    if(firstTrack7){
      countsTrackCriteria->Fill("dZLess5mm", weight);
      firstTrack7 = false;
    }
    //.................................................................................//
    //if(trackCollection[i].numberOfValidHits<7)                                                continue;
    if(firstTrack8){
      countsTrackCriteria->Fill("NOfValidHitsGreater7", weight);
      firstTrack8 = false;
    }
    //.................................................................................//
    if(trackCollection[i].hitPattern_trackerLayersWithoutMeasurement>0)                       continue;
    if(firstTrack9){
      countsTrackCriteria->Fill("NOfLostHitsMiddleEq0", weight);
      firstTrack9 = false;
    }
    //.................................................................................//
    if(trackCollection[i].trackerExpectedHitsInner_numberOfLostHits>0)                        continue;
    if(firstTrack10){
      countsTrackCriteria->Fill("NOfLostHitsInnerEq0", weight);
      firstTrack10 = false;
    }
    //.................................................................................//
    //if(trackCollection[i].hitPattern_trackerLayersWithoutMeasurement==0 && trackCollection[i].trackerExpectedHitsInner_numberOfLostHits==0)            continue;
    if(firstTrack11){
      countsTrackCriteria->Fill("missingMiddleAndInnerHits", weight);
      firstTrack11 = false;
    }
    //.................................................................................//
    if(trackCollection[i].trackRelIso03>=0.1)                                                 continue;
    if(firstTrack12){
      countsTrackCriteria->Fill("TrackIsolationDeltaR0p3Less0p1", weight);
      firstTrack12 = false;
    }
    //.................................................................................//

    outputColl.push_back(trackCollection[i]);
    //.................................................................................//
  
  }
  
  
  return outputColl;
  
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> trackParticleMatchingCuts(std::vector<evt::Track_s> trackCollection, std::vector<Jet_s>& jetColl, TH1D* countsTrackCriteria)
{

  std::vector<evt::Track_s> outputColl;

  bool firstTrack1  = true;
  bool firstTrack2  = true;
  bool firstTrack3  = true;
  bool firstTrack4  = true;

  for(unsigned int i=0; i<trackCollection.size(); i++){

    //.................................................................................//
    if(isTrackReconstructedJet(trackCollection[i], jetColl))                                  continue;
    if(firstTrack1){
      countsTrackCriteria->Fill("InJetCollectionR0p5", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedTau(&trackCollection[i]))                                          continue;
    if(firstTrack2){
      countsTrackCriteria->Fill("InTauCollectionR0p15", weight);
      firstTrack2 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedElectron(&trackCollection[i]))                                     continue;
    if(firstTrack3){
      countsTrackCriteria->Fill("InElectronCollectionR0p15", weight);
      firstTrack3 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedMuon(&trackCollection[i]))                                         continue;
    if(firstTrack4){
      countsTrackCriteria->Fill("InMuonCollectionR0p15", weight);
      firstTrack4 = false;
    }
    //.................................................................................//
    outputColl.push_back(trackCollection[i]);
    //.................................................................................//

  }
  
  
  return outputColl;
  
}
//--------------------------------------------------------------------------------------------------
TH3* loadDeDxTemplate(string path){
  TFile* InputFile = new TFile(path.c_str());
  TH3*  DeDxMap_;
  InputFile->GetObject("Charge_Vs_Path",DeDxMap_);
  if(!DeDxMap_){printf("dEdx templates in file %s can't be open\n", path.c_str()); exit(0);}

  TH3* Prob_ChargePath  = (TH3*)(DeDxMap_->Clone("Prob_ChargePath"));
  Prob_ChargePath->Reset();
  Prob_ChargePath->SetDirectory(0);

  for(int i=0;i<=Prob_ChargePath->GetXaxis()->GetNbins()+1;i++){        // pt dimension
    for(int j=0;j<=Prob_ChargePath->GetYaxis()->GetNbins()+1;j++){      // path dimension
      double Ni = 0;
      for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){Ni+=DeDxMap_->GetBinContent(i,j,k);}  
      for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){  // charge/path dimension 
	double tmp = 0;
	for(int l=0;l<=k;l++){ tmp+=DeDxMap_->GetBinContent(i,j,l);} // sum over all bins with dEdx equal or smaller than the current one
	if(Ni>0){
	  Prob_ChargePath->SetBinContent (i, j, k, tmp/Ni);
	}else{
	  Prob_ChargePath->SetBinContent (i, j, k, 0);
	}
      }
    }
  }
  InputFile->Close();
  return Prob_ChargePath;
}
//--------------------------------------------------------------------------------------------------
/* This is a sorting algorithm, which sorts the dE/dx entries of vector_prob according the radius 
   where they are deposited. The sorted vector is filled into vect_prob_sorted.
   The new vector vect_prob_sorted is sorted such that the dE/dx value deposited closest to the 
   beam pipe is the first entry. */
void sortVectorAccordingRadius(vector<double> *vect_prob_sorted, vector<double> *vect_prob, vector<double> *vect_R, vector<double> *vect_dx, double* dedx1=0, double* dedx2=0, double* dedx3=0, double* dedx4=0, double* dx1=0, double* dx2=0, double* dx3=0, double* dx4=0)
{
  int *idx = new int[vect_prob->size()];
  vector<double> aux_R (*vect_R);
  vector<double> aux_dx_sorted;

  for(unsigned int j=0; j<vect_R->size(); j++)
    {
      double aux = 10000;
      for(unsigned int i=0; i<aux_R.size(); i++)
	{ 
	  if(aux_R[i]<aux)
	    {
	      aux=aux_R[i];
	      idx[j] = i;
	    }
	}
      aux_R[idx[j]]=1000000.;
    }

  for(unsigned int j=0; j<vect_prob->size(); j++){
    vect_prob_sorted->push_back((*vect_prob)[idx[j]]);
    aux_dx_sorted.push_back((*vect_dx)[idx[j]]);
  }

  if(dedx1 && vect_R->size()>0){
    *dedx1=(*vect_prob_sorted)[0];
    *dx1  =aux_dx_sorted[0];
  }
  if(dedx2 && vect_R->size()>0){
    *dedx2=(*vect_prob_sorted)[1];
    *dx2  =aux_dx_sorted[1];
  }
  if(dedx3 && vect_R->size()>0){
    *dedx3=(*vect_prob_sorted)[2];
    *dx3  =aux_dx_sorted[2];
  }
  if(dedx4 && vect_R->size()>0){
    *dedx4=(*vect_prob_sorted)[3];
    *dx4  =aux_dx_sorted[3];
  }

  delete[] idx;
}
//--------------------------------------------------------------------------------------------------
double dEdxOnTheFly(std::vector<double> *HitsDeDx, std::vector<int> *HitsShapetest, std::vector<double> *HitsPathlength, std::vector<int> *HitsSubdetid, std::vector<double> *HitsTransverse, bool isData, TH3 *templateHistoStrip, TH3 *templateHistoPixel, bool usePixel, int nHits=0, double* dedx1=0, double* dedx2=0, double* dedx3=0, double* dedx4=0, double* dx1=0, double* dx2=0, double* dx3=0, double* dx4=0, int* measSize=0)
{

  

  const double globalPixelMC    = 3.50843;
  const double globalPixelData  = 3.5090;
  const double globalStrip      = 3.3027;
  
  int numStripMeas = 0;
  int numPixelMeas = 0;
  
  std::vector<double> vect_probs_strip;
  std::vector<double> vect_probs_pixel;
  std::vector<double> vect_probs;
  std::vector<double> vect_R;
  std::vector<double> vect_pathlength;
    
  double scaleFactorStrip = 1.;
  double scaleFactorPixel = 0.;
  if(isData) scaleFactorPixel = globalStrip/globalPixelData;
  else       scaleFactorPixel = globalStrip/globalPixelMC;

  scaleFactorPixel = scaleFactorPixel/265.;
  
  for(unsigned int j=0;j<(*HitsDeDx).size();j++){
  
    if((*HitsSubdetid)[j]>2 && !(*HitsShapetest)[j]){continue;}

    if(!usePixel && (*HitsSubdetid)[j]<=2) continue;

    double ProbStrip = 1.;
    double ProbPixel = 1.;
    
    // Strip
    if((*HitsSubdetid)[j]>2){
  
      if(templateHistoStrip){
	numStripMeas += 1;
	int    BinX   = templateHistoStrip->GetXaxis()->FindBin(50.0); // momentum bin -> not important
	int    BinY   = templateHistoStrip->GetYaxis()->FindBin((*HitsPathlength)[j]);
	int    BinZ   = templateHistoStrip->GetZaxis()->FindBin(scaleFactorStrip*(*HitsDeDx)[j]/(*HitsPathlength)[j]);
	ProbStrip     = templateHistoStrip->GetBinContent(BinX,BinY,BinZ);
      }
  
      vect_probs_strip.push_back(ProbStrip);
      vect_probs.push_back(ProbStrip);
      vect_R.push_back((*HitsTransverse)[j]);
      vect_pathlength.push_back((*HitsPathlength)[j]);
    }
    else{

      if(templateHistoPixel){
	numPixelMeas += 1;
	int    BinX   = templateHistoPixel->GetXaxis()->FindBin(50.0); // momentum bin -> not important
	int    BinY   = templateHistoPixel->GetYaxis()->FindBin((*HitsPathlength)[j]);
	int    BinZ   = templateHistoPixel->GetZaxis()->FindBin(scaleFactorPixel*(*HitsDeDx)[j]/(*HitsPathlength)[j]);
	ProbPixel     = templateHistoPixel->GetBinContent(BinX,BinY,BinZ);
      }
      vect_probs_pixel.push_back(ProbPixel);
      vect_probs.push_back(ProbPixel);
      vect_R.push_back((*HitsTransverse)[j]);
      vect_pathlength.push_back((*HitsPathlength)[j]);
    }
  }
  
  int size      = vect_probs.size();

  if(size<=5 && dedx1){
    vector<double> vector_prob_sorted;
    sortVectorAccordingRadius(&vector_prob_sorted, &vect_probs, &vect_R, &vect_pathlength, dedx1, dedx2, dedx3, dedx4, dx1, dx2, dx3, dx4);
    if(measSize) *measSize=size;
  }

  //vect_prob_sorted.erase(vect_prob_sorted.begin()+vect_prob_sorted.size()-1,vect_prob_sorted.end());

  if(nHits !=0  && size+nHits>=0){
    if(nHits<0)   size=size+nHits;
    else{
      if(size>nHits) size=nHits;
    }
    vect_probs.erase(vect_probs.begin()+size,vect_probs.end());
  }
 
  //Ias pixel + strip
  double P = 1.0/(12*size);

  // Why do I need to sort here? This makes a difference -> Yes it does
  std::sort(vect_probs.begin(), vect_probs.end(), std::less<double>() );
  for(int i=1;i<=size;i++){
    P += vect_probs[i-1] * pow(vect_probs[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
    //P += vector_prob_sorted[i-1] * pow(vector_prob_sorted[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
  }

  P *= (3.0/size);
  if(size<=0) P=-1;
  
  return P;
}
//--------------------------------------------------------------------------------------------------
std::vector<Track_s> getFakeTracksInTrackCollection(const std::vector<Track_s>& inputCollection){

  std::vector<Track_s> outputCollection;
  outputCollection.clear();
  double dRchip = 0;
  double dRchim = 0;

  for(unsigned int i=0; i<inputCollection.size(); i++){

    dRchip=10000;
    dRchim=10000;

    if(!zeroChip){
      double dPhichip = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chipGenParticle.phi));
      double dEtachip = std::abs(inputCollection[i].eta - chipGenParticle.eta);
      dRchip   = std::sqrt(pow(dPhichip,2) + pow(dEtachip,2));
    }
    if(!zeroChim){
      double dPhichim = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - chimGenParticle.phi));
      double dEtachim = std::abs(inputCollection[i].eta - chimGenParticle.eta);
      dRchim   = std::sqrt(pow(dPhichim,2) + pow(dEtachim,2));
    }
    if(dRchim>0.01 && dRchip>0.01) outputCollection.push_back(inputCollection[i]);
  }
  return outputCollection;

}
//--------------------------------------------------------------------------------------------------
// Match Tracks to a generator particle in the GenCollection
void matchTrackToGenParticle(std::vector<Track_s>& inputCollection){

  double dPhi  = 0;
  double dEta  = 0;
  double dR    = 0;
  


  for(unsigned int i=0; i<inputCollection.size(); i++){

    inputCollection[i].pdgId=0;
    inputCollection[i].status=-1;
    inputCollection[i].beta=10;
    double dRmin = 0.01;
    int    idx   = -1;

    
    for(unsigned int j=0; j<GenParticle.size(); j++){

      if(GenParticle[j].status==2) continue;

      dEta = std::abs(inputCollection[i].eta - GenParticle[j].eta);
      dPhi = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - GenParticle[j].phi));
      dR   = std::sqrt( dPhi*dPhi + dEta*dEta );
      if(dR<dRmin){
	dRmin = dR;
	idx = j;
      }
    }
    if(dRmin<0.01){
	inputCollection[i].pdgId = GenParticle[idx].pdgId;
	inputCollection[i].status= GenParticle[idx].status;
	inputCollection[i].beta  = GenParticle[idx].p/GenParticle[idx].energy;
	inputCollection[i].genPt = GenParticle[idx].pt;
	inputCollection[i].genE  = GenParticle[idx].energy;
	inputCollection[i].genEt = GenParticle[idx].et;
    }
  }
}
//--------------------------------------------------------------------------------------------------
// Match Tracks to a SimTrack
void matchTrackToSimTrack(std::vector<Track_s>& inputCollection){

  double dPhi    = 0.0;
  double dEta    = 0.0;
  double dR      = 0.0;
 

  for(unsigned int i=0; i<inputCollection.size(); i++){

    double dRmin = 0.01;
    int    idx   = -1;

    for(int j=0; j<nSimTrack; j++){

      dEta = std::abs(inputCollection[i].eta - SimTrack[j].momentum_eta);
      dPhi = std::abs(TVector2::Phi_mpi_pi(inputCollection[i].phi - SimTrack[j].momentum_phi));
      dR   = std::sqrt( dPhi*dPhi + dEta*dEta );
      
      if(dR<dRmin){
	// Match save index j
	idx=j;
	dRmin=dR;
      }
    }
    if(idx==-1){
      inputCollection[i].simEndVertexRho=-1;
      //inputCollection[i].pdgId = 0;
      continue;
    }
    //else inputCollection[i].pdgId = SimTrack[idx].type;

    // Find corresponding SimVertex where the track ends.
    double rho = -1;
    for(int j=0; j<nSimVertex; j++){

      if((unsigned int) SimVertex[j].parentIndex==SimTrack[idx].trackId){
	rho = sqrt(pow(SimVertex[j].position_x,2)+pow(SimVertex[j].position_y,2));
      }
    }
    inputCollection[i].simEndVertexRho = rho;
  }
}
//--------------------------------------------------------------------------------------------------
double getLifetime(string filename){

  unsigned found         = filename.find("ctau");
  string TargetLifetime  = filename.substr(found + 5);
  found                  = TargetLifetime.find("cm");
  TargetLifetime         = TargetLifetime.substr(0,found);
  //cout<<"Lifetime  = "<<TargetLifetime<<endl;

  return atoi(TargetLifetime.c_str());

}
/*
//--------------------------------------------------------------------------------------------------
bool isMuonTight(struct evt::MuonPFlow_s* muon){
  
  if( !muon->isGlobalMuon )                                          return false;
  if( !muon->isPFMuon )                                              return false;
  if( muon->globalTrack_chi2/muon->globalTrack_ndof>10 )             return false;
  if( muon->globalTrack_hitPattern_numberOfValidMuonHits<1 )         return false;
  if( muon->numberOfMatchedStations<2 )                              return false;
  if( muon->innerTrack_hitPattern_trackerLayersWithMeasurement<6 )   return false;
  if( muon->innerTrack_hitPattern_numberOfValidPixelHits<1 )         return false;
  if( muon->dB>=0.2 )                                                return false;
  if( std::abs( muon->vertex_z - evt::Vertex[0].z ) >= 0.5 )         return false;

  return true;
}
//--------------------------------------------------------------------------------------------------
float relPFIso(struct evt::MuonPFlow_s* muon){

  float relIso = 
    ( muon->chargedHadronIso + 
      std::max( 0.0, muon->neutralHadronIso + muon->photonIso - 0.5*muon->puChargedHadronIso ) ) / muon->pt; 

  return relIso;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::MuonPFlow_s>  getTightMuonsInEvent(){

  std::vector<evt::MuonPFlow_s> muonCollection;
  muonCollection.clear();
  
  for(unsigned int i=0; i<evt::MuonPFlow.size(); i++){

    if( evt::MuonPFlow[i].pt<=25 )                     continue;
    if( std::abs(evt::MuonPFlow[i].eta)>2.5 )          continue;
    if( !isMuonTight(&evt::MuonPFlow[i]) )             continue;
    if( relPFIso(&MuonPFlow[i])>0.12 )                   continue;

    muonCollection.push_back(evt::MuonPFlow[i]);
  }

  return muonCollection;
}
//--------------------------------------------------------------------------------------------------
float isGoodMVAElectron(struct evt::ElectronPFlow_s* electron){

  double sceta = fabs(electron->superCluster_eta);

  if(electron->pt >= 20.0) {

    if((0.0 <= sceta) && (sceta <= 0.8)) {
      if(electron->mvaTrigV0 > 0.94) {
	return true;
      }
    } else if((0.8 < sceta) && (sceta <= 1.479)) {
      if(electron->mvaTrigV0 > 0.85) {
	return true;
      }
    } else if((1.479 < sceta) && (sceta <= 2.5) ) {
      if(electron->mvaTrigV0 > 0.92) {
	return true;
      }
    }
	  
  }else if(electron->pt > 10.0) {
    if((0.0 <= sceta) && (sceta <= 0.8)) {
      if(electron->mvaTrigV0 > 0.00) {
	return true;
      }
    } else if((0.8 < sceta) && (sceta <= 1.479)) {
      if(electron->mvaTrigV0 > 0.10) {
	return true;
      }
    } else if((1.479 < sceta) && (sceta <= 2.5) ) {
      if(electron->mvaTrigV0 > 0.62) {
	return true;
      }
    } 
  }
	

  return false;
}
//--------------------------------------------------------------------------------------------------
bool isElectronTight(struct evt::ElectronPFlow_s* electron){
  
  if( !electron->passConversionVeto )                                          return false;
  if( electron->gsfTrack_trackerExpectedHitsInner_numberOfLostHits >0 )        return false;
  if( !isGoodMVAElectron(electron) )                                           return false;

  return true;
}
//--------------------------------------------------------------------------------------------------
float relPFIsoRho(struct evt::ElectronPFlow_s* electron){

  float relIsoRho =
    ( electron->chargedHadronIso + 
      std::max( 0.0, electron->neutralHadronIso + electron->photonIso - sdoublePF_value*electron->Aeff04 ) ) 
    / electron->pt; 

  return relIsoRho;
}
//--------------------------------------------------------------------------------------------------
std::vector<evt::ElectronPFlow_s>  getTightElectronsInEvent(){

  std::vector<evt::ElectronPFlow_s> electronCollection;
  electronCollection.clear();
  
  for(unsigned int i=0; i<evt::ElectronPFlow.size(); i++){

    if( evt::ElectronPFlow[i].pt<=25 )                     continue;
    if( std::abs(evt::ElectronPFlow[i].eta)>2.5 )          continue;
    if( !isElectronTight(&evt::ElectronPFlow[i]) )         continue;
    if( relPFIsoRho(&ElectronPFlow[i])>0.15 )              continue;

    electronCollection.push_back(evt::ElectronPFlow[i]);
  }

  return electronCollection;
}
*/
//--------------------------------------------------------------------------------------------------
TLorentzVector lorentzVector(float pt, float eta, float phi, float energy){
  
  TLorentzVector v4;
  v4.SetPtEtaPhiE(pt,eta,phi,energy);
  
  return v4;
};
//--------------------------------------------------------------------------------------------------
#endif

