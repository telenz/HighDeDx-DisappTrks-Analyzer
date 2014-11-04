#ifndef CHICLASS_H
#define CHICLASS_H
//-----------------------------------------------------------------------------
// File:        analyzer.h
// Description: Analyzer header for ntuples created by TheNtupleMaker
// Created:     Wed Sep  3 16:58:58 2014 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
// -- System

#include <stdio.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cassert>
#include <cmath>

#include "analyzerutil.h"
#include "treestream.h"
#include "pdg.h"

// -- Root

#include "TROOT.h"
#include "TApplication.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TH1F.h"
#include "TH2F.h"
//-----------------------------------------------------------------------------
struct ChiTrack_s
{
  // RecoTrack information
  double	pt;
  double	px;
  double	py;
  double	pz;
  double	phi;
  double	eta;
  double	vx;
  double	vy;
  double	vz;
  unsigned short	numberOfValidHits;
  unsigned short	hitPattern_trackerLayersWithoutMeasurement;
  unsigned short	trackerExpectedHitsInner_numberOfLostHits;
  unsigned short	trackerExpectedHitsOuter_numberOfHits;
  int	trackHighPurity;
  double	trackRelIso03;
  double	caloEMDeltaRp3;
  double	caloHadDeltaRp3;
  double	caloEMDeltaRp4;
  double	caloHadDeltaRp4;
  double	caloEMDeltaRp5;
  double	caloHadDeltaRp5;
  double	dEdxNPHarm2;
  double	dEdxNPTru40;
  unsigned int	dEdxNPNoM;
  double	dEdxHarm2;
  double	dEdxTru40;
  unsigned int	dEdxNoM;
  double	dEdxHitsHarm2_1000;
  double	dEdxHitsHarm2_7;
  double	dEdxHitsHarm2_5;
  double	dEdxHitsHarm2_3;
  double	dEdxHitsHarm2_2;
  double	dEdxHitsHarm2_1;
  double	dEdxHitsTrun40_1000;
  double	dEdxHitsTrun40_7;
  double	dEdxHitsTrun40_5;
  double	dEdxHitsTrun40_3;
  double	dEdxHitsTrun40_2;
  double	dEdxHitsTrun40_1;
  double	dEdxHitsMedian_1000;
  double	dEdxHitsMedian_7;
  double	dEdxHitsMedian_5;
  double	dEdxHitsMedian_3;
  double	dEdxHitsMedian_2;
  double	dEdxHitsMedian_1;
  double        beta;
  bool          matched;
  // GenParticle Information
  int	        gencharge;
  double	genp;
  double	genenergy;
  double	genet;
  double	genpz;
  double	genpt;
  double	genphi;
  double	geneta;
  double	genmass;
  int	        genpdgId;
  // SimTrack information
  float	        SimTrackcharge;
  int	        SimTrackgenpartIndex;
  double	SimTrackmomentum_energy;
  double	SimTrackmomentum_eta;
  double	SimTrackmomentum_phi;
  double	SimTrackmomentum_pt;
  int	        SimTracknoGenpart;
  int	        SimTracknoVertex;
  unsigned int	SimTracktrackId;
  int	        SimTracktype;
  int	        SimTrackvertIndex;
  // Decay Vertex
  int	        SimVertexparentIndex;
  int	        SimVertexnoParent;
  unsigned int	SimVertexvertexId;
  double	SimVertexposition_x;
  double	SimVertexposition_y;
  double	SimVertexposition_z;
  double	SimVertexposition_t;
  bool          SimVertexFound;
};
std::vector<ChiTrack_s> ChiTrack(2);

//-----------------------------------------------------------------------------
inline void fillChiTrackWithRecoTrackVariables(struct ChiTrack_s* chi, struct evt::Track_s* track)
{
  chi->pt	= track->pt;
  chi->px	= track->px;
  chi->py	= track->py;
  chi->pz	= track->pz;
  chi->phi	= track->phi;
  chi->eta	= track->eta;
  chi->vx	= track->vx;
  chi->vy	= track->vy;
  chi->vz	= track->vz;
  chi->numberOfValidHits	                        = track->numberOfValidHits;
  chi->hitPattern_trackerLayersWithoutMeasurement	= track->hitPattern_trackerLayersWithoutMeasurement;
  chi->trackerExpectedHitsInner_numberOfLostHits	= track->trackerExpectedHitsInner_numberOfLostHits;
  chi->trackerExpectedHitsOuter_numberOfHits	        = track->trackerExpectedHitsOuter_numberOfHits;
  chi->trackHighPurity	= track->trackHighPurity;
  chi->trackRelIso03	= track->trackRelIso03;
  chi->caloEMDeltaRp3	= track->caloEMDeltaRp3;
  chi->caloHadDeltaRp3	= track->caloHadDeltaRp3;
  chi->caloEMDeltaRp4	= track->caloEMDeltaRp4;
  chi->caloHadDeltaRp4	= track->caloHadDeltaRp4;
  chi->caloEMDeltaRp5	= track->caloEMDeltaRp5;
  chi->caloHadDeltaRp5	= track->caloHadDeltaRp5;
  chi->dEdxNPHarm2  	= track->dEdxNPHarm2;
  chi->dEdxNPTru40	        = track->dEdxNPTru40;
  chi->dEdxNPNoM    	        = track->dEdxNPNoM;
  chi->dEdxHarm2    	        = track->dEdxHarm2;
  chi->dEdxTru40    	        = track->dEdxTru40;
  chi->dEdxNoM      	        = track->dEdxNoM;
  chi->dEdxHitsHarm2_1000	= track->dEdxHitsHarm2_1000;
  chi->dEdxHitsHarm2_7	        = track->dEdxHitsHarm2_7;
  chi->dEdxHitsHarm2_5	        = track->dEdxHitsHarm2_5;
  chi->dEdxHitsHarm2_3	        = track->dEdxHitsHarm2_3;
  chi->dEdxHitsHarm2_2	        = track->dEdxHitsHarm2_2;
  chi->dEdxHitsHarm2_1	        = track->dEdxHitsHarm2_1;
  chi->dEdxHitsTrun40_1000	= track->dEdxHitsTrun40_1000;
  chi->dEdxHitsTrun40_7	        = track->dEdxHitsTrun40_7;
  chi->dEdxHitsTrun40_5	        = track->dEdxHitsTrun40_5;
  chi->dEdxHitsTrun40_3	        = track->dEdxHitsTrun40_3;
  chi->dEdxHitsTrun40_2	        = track->dEdxHitsTrun40_2;
  chi->dEdxHitsTrun40_1	        = track->dEdxHitsTrun40_1;
  chi->dEdxHitsMedian_1000	= track->dEdxHitsMedian_1000;
  chi->dEdxHitsMedian_7	        = track->dEdxHitsMedian_7;
  chi->dEdxHitsMedian_5	        = track->dEdxHitsMedian_5;
  chi->dEdxHitsMedian_3	        = track->dEdxHitsMedian_3;
  chi->dEdxHitsMedian_2	        = track->dEdxHitsMedian_2;
  chi->dEdxHitsMedian_1	        = track->dEdxHitsMedian_1;
}

inline void fillChiTrackWithSimTrackVariables(struct ChiTrack_s* chi, struct evt::SimTrack_s* track)
{
  chi->SimTrackcharge	        = track->charge;
  chi->SimTrackvertIndex        = track->vertIndex;
  chi->SimTracknoVertex	        = track->noVertex;
  chi->SimTrackgenpartIndex     = track->genpartIndex;
  chi->SimTracknoGenpart	= track->noGenpart;
  chi->SimTracktype	        = track->type;
  chi->SimTracktrackId	        = track->trackId;
  chi->SimTrackmomentum_pt	= track->momentum_pt;
  chi->SimTrackmomentum_phi	= track->momentum_phi;
  chi->SimTrackmomentum_eta	= track->momentum_eta;
  chi->SimTrackmomentum_energy	= track->momentum_energy;
}

inline void fillChiTrackWithSimVertexVariables(struct ChiTrack_s* chi, struct evt::SimVertex_s* vertex)
{
  chi->SimVertexparentIndex	= vertex->parentIndex;
  chi->SimVertexnoParent	= vertex->noParent;
  chi->SimVertexvertexId	= vertex->vertexId;
  chi->SimVertexposition_x	= vertex->position_x;
  chi->SimVertexposition_y	= vertex->position_y;
  chi->SimVertexposition_z	= vertex->position_z;
  chi->SimVertexposition_t	= vertex->position_t;
}

inline void fillChiTrackWithGenParticleVariables(struct ChiTrack_s* chi, struct evt::GenParticle_s* genParticle)
{
  chi->gencharge    = genParticle->charge;
  chi->genp	    = genParticle->p;
  chi->genenergy    = genParticle->energy;
  chi->genet        = genParticle->et;
  chi->genpz	    = genParticle->pz;
  chi->genpt	    = genParticle->pt;
  chi->genphi  	    = genParticle->phi;
  chi->geneta  	    = genParticle->eta;
  chi->genmass 	    = genParticle->mass;
  chi->genpdgId	    = genParticle->pdgId;
  chi->matched	    = false;
 
}
#endif
