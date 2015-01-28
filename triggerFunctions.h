#ifndef TRIGGERFUNCTIONS_H
#define TRIGGERFUNCTIONS_H
//-----------------------------------------------------------------------------
#include <iostream>
#include "TVector2.h"
#include "chiClass.h"
#include "declarationsOfClasses.h"
#include "TH3.h"
using namespace std;
using namespace evt;
//-----------------------------------------------------------------------------
int getTriggerResult(TString trigger){

  int result = 0;

  if(trigger == "edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v"){

    if( 
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v1    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v6    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v7    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v8    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v9    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v10   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v11   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v12   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v13   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v14   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v15   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v16   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v17   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v18   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v19   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v20   == 1 )
      {result     =  1; }
    else if(  
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v1    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v6    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v7    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v8    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v9    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v10   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v11   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v12   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v13   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v14   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v15   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v16   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v17   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v18   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v19   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v20   == -1 )
      {     result     =  -1; }
  }

  if(trigger == "edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v"){
    
    if( 
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v2    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v6    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v7    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v8    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v9    == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v10   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v11   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v12   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v13   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v14   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v15   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v16   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v17   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v18   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v19   == 1 ||
       edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v20   == 1 )
      {result     =  1; }
    else if(  
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v2    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v6    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v7    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v8    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v9    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v10   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v11   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v12   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v13   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v14   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v15   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v16   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v17   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v18   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v19   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v20   == -1 )
      {     result     =  -1; }
  }

  if(trigger == "edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v"){

    if(  
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v1                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v2                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v3                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v4                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v5                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v6                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v7                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v8                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v9                    == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v10                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v11                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v12                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v13                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v14                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v15                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v16                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v17                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v18                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v19                   == 1 ||
       edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v20                   == 1 )
      {result                     =  1; }
    else if(  
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v1                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v2                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v3                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v4                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v5                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v6                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v7                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v8                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v9                    == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v10                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v11                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v12                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v13                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v14                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v15                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v16                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v17                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v18                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v19                   == -1 &&
	    edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v20                   == -1 )
      {    result                    =  -1; }
  }


  return result;


}

//-----------------------------------------------------------------------------
int getTriggerPrescales(TString trigger){

  vector<int> MonoCentralPFJet80_PFMETnoMu95_NHEF0p95;
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v1);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v6);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v7);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v8);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v9);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v10);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v11);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v12);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v13);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v14);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v15);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v16);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v17);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v18);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v19);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v20);
  vector<int> MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale;
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v1);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v2);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v3);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v4);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v5);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v6);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v7);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v8);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v9);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v10);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v11);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v12);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v13);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v14);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v15);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v16);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v17);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v18);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v19);
  MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v20);

  vector<int> MonoCentralPFJet80_PFMETnoMu105_NHEF0p95;
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v2);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v6);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v7);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v8);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v9);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v10);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v11);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v12);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v13);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v14);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v15);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v16);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v17);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v18);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v19);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.push_back(edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v20);
  vector<int> MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale;
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v1);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v2);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v3);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v4);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v5);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v6);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v7);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v8);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v9);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v10);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v11);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v12);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v13);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v14);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v15);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v16);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v17);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v18);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v19);
  MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v20);

  vector<int> MET120_HBHENoiseCleaned;
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v1);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v2);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v3);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v4);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v5);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v6);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v7);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v8);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v9);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v10);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v11);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v12);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v13);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v14);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v15);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v16);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v17);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v18);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v19);
  MET120_HBHENoiseCleaned.push_back(edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v20);
  vector<int> MET120_HBHENoiseCleaned_prescale;
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v1);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v2);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v3);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v4);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v5);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v6);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v7);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v8);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v9);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v10);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v11);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v12);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v13);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v14);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v15);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v16);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v17);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v18);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v19);
  MET120_HBHENoiseCleaned_prescale.push_back(edmTriggerResultsHelper_prescale_HLT_MET120_HBHENoiseCleaned_v20);



  if(trigger == "edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_v"){

    for(unsigned int i=0; i<MonoCentralPFJet80_PFMETnoMu95_NHEF0p95.size(); i++){ if(MonoCentralPFJet80_PFMETnoMu95_NHEF0p95[i]==1){ return MonoCentralPFJet80_PFMETnoMu95_NHEF0p95_prescale[i];}}

  }
  if(trigger == "edmTriggerResultsHelper_value_HLT_MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_v"){

    for(unsigned int i=0; i<MonoCentralPFJet80_PFMETnoMu105_NHEF0p95.size(); i++){ if(MonoCentralPFJet80_PFMETnoMu105_NHEF0p95[i]==1){ return MonoCentralPFJet80_PFMETnoMu105_NHEF0p95_prescale[i];}}

  }


  if(trigger == "edmTriggerResultsHelper_value_HLT_MET120_HBHENoiseCleaned_v"){

    for(unsigned int i=0; i<MET120_HBHENoiseCleaned.size(); i++){ 
      if(MET120_HBHENoiseCleaned[i]==1) return MET120_HBHENoiseCleaned_prescale[i];
    }

  }

  return 0;

}

//--------------------------------------------------------------------------------------------------
#endif

