#ifndef DTSelection_H
#define DTSelection_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
#include "TLorentzVector.h"
#include "chiClass.h"
#include "declarationsOfClasses.h"
#include "functions.h"
#include "TH3.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
std::vector<evt::Track_s> DTSelectionCuts(std::vector<evt::Track_s> trackCollection, std::vector<Jet_s>& jetColl, TH1D* countsTrackCriteria )
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
  bool firstTrack10  = true;
  bool firstTrack11  = true;
  bool firstTrack12  = true;
  bool firstTrack13  = true;
  bool firstTrack14  = true;
  bool firstTrack15  = true;
  bool firstTrack16  = true;
  bool firstTrack17  = true;
  bool firstTrack18  = true;
  bool firstTrack19  = true;
  bool firstTrack20  = true;
  bool firstTrack21  = true;


  for(unsigned int i=0; i<trackCollection.size(); i++){
    if(firstTrack1){
      countsTrackCriteria->Fill("#geq 1 recon. trk", weight);
      firstTrack1 = false;
    }
    //.................................................................................//
    if(trackCollection[i].pt<=50.)                                            continue;
    if(firstTrack2){
      countsTrackCriteria->Fill("#geq 1 trk with pt > 50GeV", weight);  
      firstTrack2 = false;
    }
    //.................................................................................//
    if(std::abs(trackCollection[i].eta)>2.1)                                  continue;
    if(firstTrack3){
      countsTrackCriteria->Fill("#geq 1 trk with |eta| > 2.1", weight);
      firstTrack3 = false;
    }
    //.................................................................................//
    if(std::abs(trackCollection[i].eta)>1.42 && std::abs(trackCollection[i].eta)<1.65)        continue;
    if(firstTrack4){
      countsTrackCriteria->Fill("#geq 1 trk not with 1.42<|#eta|<1.65", weight);
      firstTrack4 = false;
    }
    //.................................................................................//
    if(std::abs(trackCollection[i].eta)>0.15 && std::abs(trackCollection[i].eta)<0.35)        continue;
    if(firstTrack5){
      countsTrackCriteria->Fill("#geq 1 trk not with 0.15<|#eta|<0.35", weight);
      firstTrack5 = false;
    }
    //.................................................................................//
    if(std::abs(trackCollection[i].eta)>1.55 && std::abs(trackCollection[i].eta)<1.85)        continue;
    if(firstTrack6){
      countsTrackCriteria->Fill("#geq 1 trk not with 1.55<|#eta|<1.85", weight);
      firstTrack6 = false;
    }
    //.................................................................................//
    if(getTrkIsMatchedDeadEcal(&trackCollection[i]))                                          continue;
    if(firstTrack7){
      countsTrackCriteria->Fill("#geq 1 trk not matched to dead ECAL cell", weight);
      firstTrack7 = false;
    }
    //.................................................................................//
    if(isWithinIntermoduleGapsOfECAL(&trackCollection[i]))                                    continue;
    if(firstTrack8){
      countsTrackCriteria->Fill("#geq 1 trk not within intermodule gaps", weight);
      firstTrack8 = false;
    }
   
    //.................................................................................//
    if(getTrkIsMatchedBadCSC(&trackCollection[i]))                                            continue;
    if(firstTrack9){
      countsTrackCriteria->Fill("#geq 1 trk not matched to bad CSC cell", weight);
      firstTrack9 = false;
    }
    //.................................................................................//
    double _dvx = trackCollection[i].vx - Vertex[0].x;
    double _dvy = trackCollection[i].vy - Vertex[0].y;
    double d0 = abs( - _dvx*trackCollection[i].py + _dvy*trackCollection[i].px)/trackCollection[i].pt;
    if(abs(d0)>0.02)                                                                          continue;
    if(firstTrack10){
      countsTrackCriteria->Fill("#geq 1 trk with |d0|<0.2mm", weight);
      firstTrack10 = false;
    }
    //.................................................................................//
    double _dvz = trackCollection[i].vz - Vertex[0].z;
    double dZ = _dvz - ( _dvx*trackCollection[i].px + _dvy*trackCollection[i].py)/trackCollection[i].pt * (trackCollection[i].pz/trackCollection[i].pt);
    if(abs(dZ)>0.5)                                                                           continue;
    if(firstTrack11){
      countsTrackCriteria->Fill("#geq 1 trk with |dz|<5mm", weight);
      firstTrack11 = false;
    }
    //.................................................................................//
    if(trackCollection[i].numberOfValidHits<7)                       continue;
    if(firstTrack12){
      countsTrackCriteria->Fill("#geq 1 trk with N_{hits} > 6", weight);
      firstTrack12 = false;
    }
    //.................................................................................//
    if(trackCollection[i].hitPattern_trackerLayersWithoutMeasurement>0)                       continue;
    if(firstTrack13){
      countsTrackCriteria->Fill("#geq 1 trk with N_{lost}^{middle} = 0", weight);
      firstTrack13 = false;
    }
    //.................................................................................//
    if(trackCollection[i].trackerExpectedHitsInner_numberOfLostHits>0)                        continue;
    if(firstTrack14){
      countsTrackCriteria->Fill("#geq 1 trk with N_{lost}^{inner} = 0", weight);
      firstTrack14 = false;
    }
    //.................................................................................//
    if(trackCollection[i].trackRelIso03>=0.05)                                   continue;
    if(firstTrack15){
      countsTrackCriteria->Fill("#geq 1 isolated trk", weight);
      firstTrack15 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedJet(&trackCollection[i], jetColl))                                  continue;
    if(firstTrack16){
      countsTrackCriteria->Fill("#geq 1 trk not matched to jet", weight);
      firstTrack16 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedTau(&trackCollection[i]))                                          continue;
    if(firstTrack17){
      countsTrackCriteria->Fill("#geq 1 trk not matched to tau", weight);
      firstTrack17 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedElectron(&trackCollection[i]))                                     continue;
    if(firstTrack18){
      countsTrackCriteria->Fill("#geq 1 trk not matched to electron", weight);
      firstTrack18 = false;
    }
    //.................................................................................//
    if(isTrackReconstructedMuon(&trackCollection[i]))                                         continue;
    if(firstTrack19){
      countsTrackCriteria->Fill("#geq 1 trk not matched to muon", weight);
      firstTrack19 = false;
    }
    //.................................................................................//
    if(trackCaloIsolation(&trackCollection[i])>10)                            continue;
    if(firstTrack20){
      countsTrackCriteria->Fill("#geq 1 trk with E_{calo}<5GeV", weight);
      firstTrack20 = false;
    }
    //.................................................................................//
    if(trackCollection[i].trackerExpectedHitsOuter_numberOfHits<3)            continue;
    if(firstTrack21){
      countsTrackCriteria->Fill("NOfLostHitsOuterGeq3", weight);
      firstTrack21 = false;
    }
    //.................................................................................//
    outputColl.push_back(trackCollection[i]);
    //.................................................................................//
  
  }
    
  return outputColl;
   
}
//--------------------------------------------------------------------------------------------------
#endif

