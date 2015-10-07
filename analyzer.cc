//-----------------------------------------------------------------------------
// File:        analyzer.cc
// Description: Analyzer for ntuples created by TheNtupleMaker
// Created:     Tue Sep 24 14:30:05 2013 by mkanalyzer.py
// Author:      Teresa Lenz
//-----------------------------------------------------------------------------
#include "config.h"
#include "analyzer.h"
#include "functions.h"
#include "selection.h"
#include "hitInformation.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TStyle.h"
#include <time.h>
#include <fstream>
#include <iomanip>
using namespace std;
using namespace evt;

//-----------------------------------------------------------------------------
int main(int argc, char** argv)
{

  // Get file list and histogram filename from command line
  clock_t t;
  t = clock();
  commandLine cmdline;
  decodeCommandLine(argc, argv, cmdline);

  // Get names of ntuple files to be processed and open chain of ntuples
  isSignal      = false;
  isSignalC1N1  = false;
  vector<string> filenames = getFilenames(cmdline.filelist);
  if(filenames[0].find("pMSSM") != std::string::npos || filenames[0].find("RECO_RAW2DIGI_L1Reco_RECO") != std::string::npos) isSignal = true;
  if(filenames[0].find("signal_C1N1") != std::string::npos && filenames[0].find("RECO_RAW2DIGI_L1Reco_RECO") != std::string::npos) isSignalC1N1 = true;
  string outputfilename = cmdline.outputfilename;

  double TargetLifetime  = 0;
  double CurrentLifetime = 0;
  if(isSignal){
    TargetLifetime = getLifetime(outputfilename);
    cout<<endl<<"TargetLifetime = "<<TargetLifetime<<endl<<endl;
  }
  bool isData = false;
  if(filenames[0].find("MET_Run2012") != std::string::npos) isData = true;
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  //cout << "Number of events: " << nevents << endl;

  // *************************************************************************************************************************************************************
  // ************************************************************** Systematic uncertainties *********************************************************************
  bool up       = false;
  bool down     = false;
  // *************************************************************************************************************************************************************
  bool ISRunc     = false;
  bool PUunc      = false;
  bool JECunc     = false;
  bool JERunc     = false;
  bool TrigEffunc = false;
  bool PDFunc     = false;
  // ************************************************************** Systematic uncertainties *********************************************************************
  // *************************************************************************************************************************************************************
  // Select variables to be read
  selectVariables(stream);

  // The root application is needed to make canvases visible during
  // program execution. If this is not needed, just comment out the
  // following line

  // TApplication app("analyzer", &argc, argv);

  /**
     Notes 1
     -------
     1. Use
     ofile = outputFile(cmdline.outputfile, stream)

     to skim events to output file in addition to writing out histograms.

     2. Use
     ofile.addEvent(event-weight)

     to specify that the current event is to be added to the output file.
     If omitted, the event-weight is defaulted to 1.

     3. Use
     ofile.count(cut-name, event-weight)

     to keep track, in the count histogram, of the number of events
     passing a given cut. If omitted, the event-weight is taken to be 1.
     If you want the counts in the count histogram to appear in a given
     order, specify the order, before entering the event loop, as in
     the example below

     ofile.count("NoCuts", 0)
     ofile.count("GoodEvent", 0)
     ofile.count("Vertex", 0)
     ofile.count("MET", 0)

     Notes 2
     -------
     By default all variables are saved. Before the event loop, you can use
  
     select(objectname)
	  
     e.g.,
	
     select("GenParticle")
  
     to declare that you intend to select objects of this type. The
     selection is done using

     select(objectname, index)
	  
     e.g.,
	  
     select("GenParticle", 3),
  
     which is called within the event loop. Call saveSelectedObjects()
     before a call to addEvent if you wish to save the selected objects.
     All other objects are saved by default.
	 
     NB: If you declare your intention to select objects of a given type
     by calling select(objectname), but subsequently fail to select
     them using select(objectname, index) then none will be saved!
  */
  

  outputFile ofile(cmdline.outputfilename);

  //-------------------------------------------------------------------------------------
  // Declare pileup histograms
  //-------------------------------------------------------------------------------------
  gStyle->SetOptStat(111111);

  TH1D *nVertices_0            = new TH1D("nVertices_0","nVertices_0",100,0,100);
  TH1D *hPU_NumInteractions_0  = new TH1D("hPU_NumInteractions_0","hPU_NumInteractions_0",100,0,100);
  TH1D *hTrueNumInteractions_0 = new TH1D("hTrueNumInteractions_0","hTrueNumInteractions_0",100,0,100);
  //-------------------------------------------------------------------------------------
  // PDF stuff
  //-------------------------------------------------------------------------------------
  unsigned int  nbins_pdfweight = 0;
  vector<double> pdf_weight_sum(50);
  Event PDFselection[45] = {Event("chiTrackspreselectionTriggerPDF_0",ofile ), Event("chiTrackspreselectionTriggerPDF_1",ofile ), Event("chiTrackspreselectionTriggerPDF_2",ofile ), Event("chiTrackspreselectionTriggerPDF_3",ofile ), Event("chiTrackspreselectionTriggerPDF_4",ofile ), Event("chiTrackspreselectionTriggerPDF_5",ofile ), Event("chiTrackspreselectionTriggerPDF_6",ofile ), Event("chiTrackspreselectionTriggerPDF_7",ofile ), Event("chiTrackspreselectionTriggerPDF_8",ofile ), Event("chiTrackspreselectionTriggerPDF_9",ofile ), Event("chiTrackspreselectionTriggerPDF_10",ofile ), Event("chiTrackspreselectionTriggerPDF_11",ofile ), Event("chiTrackspreselectionTriggerPDF_12",ofile ), Event("chiTrackspreselectionTriggerPDF_13",ofile ), Event("chiTrackspreselectionTriggerPDF_14",ofile ), Event("chiTrackspreselectionTriggerPDF_15",ofile ), Event("chiTrackspreselectionTriggerPDF_16",ofile ), Event("chiTrackspreselectionTriggerPDF_17",ofile ), Event("chiTrackspreselectionTriggerPDF_18",ofile ), Event("chiTrackspreselectionTriggerPDF_19",ofile ), Event("chiTrackspreselectionTriggerPDF_20",ofile ), Event("chiTrackspreselectionTriggerPDF_21",ofile ), Event("chiTrackspreselectionTriggerPDF_22",ofile ), Event("chiTrackspreselectionTriggerPDF_23",ofile ), Event("chiTrackspreselectionTriggerPDF_24",ofile ), Event("chiTrackspreselectionTriggerPDF_25",ofile ), Event("chiTrackspreselectionTriggerPDF_26",ofile ), Event("chiTrackspreselectionTriggerPDF_27",ofile ), Event("chiTrackspreselectionTriggerPDF_28",ofile ), Event("chiTrackspreselectionTriggerPDF_29",ofile ), Event("chiTrackspreselectionTriggerPDF_30",ofile ), Event("chiTrackspreselectionTriggerPDF_31",ofile ), Event("chiTrackspreselectionTriggerPDF_32",ofile ), Event("chiTrackspreselectionTriggerPDF_33",ofile ), Event("chiTrackspreselectionTriggerPDF_34",ofile ), Event("chiTrackspreselectionTriggerPDF_35",ofile ), Event("chiTrackspreselectionTriggerPDF_36",ofile ), Event("chiTrackspreselectionTriggerPDF_37",ofile ), Event("chiTrackspreselectionTriggerPDF_38",ofile ), Event("chiTrackspreselectionTriggerPDF_39",ofile ), Event("chiTrackspreselectionTriggerPDF_40",ofile ), Event("chiTrackspreselectionTriggerPDF_41",ofile ), Event("chiTrackspreselectionTriggerPDF_42",ofile ), Event("chiTrackspreselectionTriggerPDF_43",ofile ), Event("chiTrackspreselectionTriggerPDF_44",ofile )} ;
  if(PDFunc && isSignal){

    nbins_pdfweight = 45;        //pdfWeights_cteq66.size();
    
    for (unsigned int i = 0; i < nbins_pdfweight; i++){ 
    
      if(isSignal) PDFselection[i].onlyChi= true;
      PDFselection[i].triggerRequirements = true;
      PDFselection[i].trigger             = true;
      PDFselection[i].trackPreselection   = true;
      PDFselection[i].qcdSupression       = true;
    }

    // 1) Run over all events to get the sum of pdf weights (no overall xsection uncertainty needed)
    for(int entry=0; entry <nevents; ++entry)
      {
	// Read event into memory
	stream.read(entry);
	stream._chain->GetEntry(entry);
	initialize();
	fillObjects();

	for (unsigned int i = 0; i < nbins_pdfweight; i++){ 
	  pdf_weight_sum[i] += (pdfWeights_cteq66[i]/pdfWeights_cteq66[0]);
	}
      }
  }

  //-------------------------------------------------------------------------------------
  // Declaration of Variables
  //-------------------------------------------------------------------------------------
  cout<<endl<<endl<<"------------ Is it signal : "<<isSignal<<" --------------"<<endl<<endl;
  cout<<endl<<endl<<"------------ Is it C1N1   : "<<isSignalC1N1<<" --------------"<<endl<<endl;
  cout<<endl<<endl<<"------------ Is it data   : "<<isData<<" --------------"<<endl<<endl;
  /*
  class Event noSelection("noSelection",ofile);
  class Event triggerRequirements("triggerRequirements",ofile);
  triggerRequirements.triggerRequirements = true;
  class Event triggerRequirementsTrigger("triggerRequirementsTrigger",ofile);
  triggerRequirementsTrigger.triggerRequirements = true;
  triggerRequirementsTrigger.trigger             = true;
  class Event preselection("preselection",ofile);
  preselection.triggerRequirements        = true;
  preselection.trackPreselection          = true;
  preselection.qcdSupression              = true;
  class Event fullSelection("fullSelection",ofile);
  fullSelection.triggerRequirements       = true;
  fullSelection.trackPreselection         = true;
  fullSelection.qcdSupression             = true;
  fullSelection.trackCandidateCutFinal    = true;
  fullSelection.TrackPtRequirement        = true;
  fullSelection.NumOfLostOuterRequirement = true;
  fullSelection.CaloIsolationRequirement  = true;
  fullSelection.DeDxRequirement           = true;
  */
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Only Chi
  class Event chiTracksnoSelection("chiTracksnoSelection",ofile);
  if(isSignal) chiTracksnoSelection.onlyChi = true;

  class Event chiTrackstriggerRequirements("chiTrackstriggerRequirements",ofile);
  if(isSignal) chiTrackstriggerRequirements.onlyChi = true;
  chiTrackstriggerRequirements.triggerRequirements  = true;

  class Event chiTrackstriggerRequirementsTrigger("chiTrackstriggerRequirementsTrigger",ofile);
  if(isSignal) chiTrackstriggerRequirementsTrigger.onlyChi = true;
  chiTrackstriggerRequirementsTrigger.triggerRequirements  = true;
  chiTrackstriggerRequirementsTrigger.trigger              = true;

  class Event chiTracksQCDsupression("chiTracksQCDsupression",ofile);
  if(isSignal) chiTracksQCDsupression.onlyChi = true;
  chiTracksQCDsupression.triggerRequirements  = true;
  chiTracksQCDsupression.qcdSupression        = true;

  class Event chiTracksQCDsupressionTrigger("chiTracksQCDsupressionTrigger",ofile);
  if(isSignal) chiTracksQCDsupressionTrigger.onlyChi = true;
  chiTracksQCDsupressionTrigger.triggerRequirements  = true;
  chiTracksQCDsupressionTrigger.trigger              = true;
  chiTracksQCDsupressionTrigger.qcdSupression        = true;

  class Event chiTrackspreselection("chiTrackspreselection",ofile);
  if(isSignal) chiTrackspreselection.onlyChi= true;
  chiTrackspreselection.triggerRequirements = true;
  chiTrackspreselection.trackPreselection   = true;
  chiTrackspreselection.qcdSupression       = true;

  class Event chiTrackspreselectionTrigger("chiTrackspreselectionTrigger",ofile);
  if(isSignal) chiTrackspreselectionTrigger.onlyChi= true;
  chiTrackspreselectionTrigger.triggerRequirements = true;
  chiTrackspreselectionTrigger.trigger             = true;
  chiTrackspreselectionTrigger.trackPreselection   = true;
  chiTrackspreselectionTrigger.qcdSupression       = true;

  class Event chiTrackspreselectionNoTrigger("chiTrackspreselectionNoTrigger",ofile);
  if(isSignal) chiTrackspreselectionNoTrigger.onlyChi= true;
  chiTrackspreselectionNoTrigger.trackPreselection   = true;
  chiTrackspreselectionNoTrigger.qcdSupression       = true;

  class Event chiTrackspreselectionNoQCDCuts("chiTrackspreselectionNoQCDCuts",ofile);
  if(isSignal) chiTrackspreselectionNoQCDCuts.onlyChi= true;
  chiTrackspreselectionNoQCDCuts.triggerRequirements = true;
  chiTrackspreselectionNoQCDCuts.trackPreselection   = true;

  class Event chiTrackspreselectionTriggerNoQCDCuts("chiTrackspreselectionTriggerNoQCDCuts",ofile);
  if(isSignal) chiTrackspreselectionTriggerNoQCDCuts.onlyChi= true;
  chiTrackspreselectionTriggerNoQCDCuts.triggerRequirements = true;
  chiTrackspreselectionTriggerNoQCDCuts.trigger             = true;
  chiTrackspreselectionTriggerNoQCDCuts.trackPreselection   = true;

   class Event chiTrackspreselectionNoQCDCutsNoTrigger("chiTrackspreselectionNoQCDCutsNoTrigger",ofile);
  if(isSignal) chiTrackspreselectionNoQCDCutsNoTrigger.onlyChi= true;
  chiTrackspreselectionNoQCDCutsNoTrigger.trackPreselection   = true;

  class Event chiTrackspreselectionWjets("chiTrackspreselectionWjets",ofile);
  if(isSignal) chiTrackspreselectionWjets.onlyChi= true;
  chiTrackspreselectionWjets.triggerRequirements = true;
  chiTrackspreselectionWjets.trackPreselection   = true;
  chiTrackspreselectionWjets.qcdSupression       = true;
  chiTrackspreselectionWjets.isolatedLeptonCut   = true;

  class Event chiTrackspreselectionPlusNLostGt1("chiTrackspreselectionPlusNLostGt1",ofile);
  if(isSignal) chiTrackspreselectionPlusNLostGt1.onlyChi      = true;
  chiTrackspreselectionPlusNLostGt1.triggerRequirements       = true;
  chiTrackspreselectionPlusNLostGt1.trackPreselection         = true;
  chiTrackspreselectionPlusNLostGt1.qcdSupression             = true;
  chiTrackspreselectionPlusNLostGt1.trackCandidateCutFinal    = true;
  chiTrackspreselectionPlusNLostGt1.TrackPtRequirement        = false;
  chiTrackspreselectionPlusNLostGt1.NumOfLostOuterRequirement = true;
  chiTrackspreselectionPlusNLostGt1.CaloIsolationRequirement  = false;
  chiTrackspreselectionPlusNLostGt1.DeDxRequirement           = false;

  class Event chiTrackspreselection1LostPlusIsoCut("chiTrackspreselection1LostPlusIsoCut",ofile);
  if(isSignal) chiTrackspreselection1LostPlusIsoCut.onlyChi      = true;
  chiTrackspreselection1LostPlusIsoCut.triggerRequirements       = true;
  chiTrackspreselection1LostPlusIsoCut.trackPreselection         = true;
  chiTrackspreselection1LostPlusIsoCut.qcdSupression             = true;
  chiTrackspreselection1LostPlusIsoCut.trackCandidateCutFinal    = true;
  chiTrackspreselection1LostPlusIsoCut.TrackPtRequirement        = false;
  chiTrackspreselection1LostPlusIsoCut.NumOfLostOuterRequirement = true;
  chiTrackspreselection1LostPlusIsoCut.CaloIsolationRequirement  = true;
  chiTrackspreselection1LostPlusIsoCut.DeDxRequirement           = false;

  class Event chiTrackspreselection1LostPlusPtCut("chiTrackspreselection1LostPlusPtCut",ofile);
  if(isSignal) chiTrackspreselection1LostPlusPtCut.onlyChi      = true;
  chiTrackspreselection1LostPlusPtCut.triggerRequirements       = true;
  chiTrackspreselection1LostPlusPtCut.trackPreselection         = true;
  chiTrackspreselection1LostPlusPtCut.qcdSupression             = true;
  chiTrackspreselection1LostPlusPtCut.trackCandidateCutFinal    = true;
  chiTrackspreselection1LostPlusPtCut.TrackPtRequirement        = true;
  chiTrackspreselection1LostPlusPtCut.NumOfLostOuterRequirement = true;
  chiTrackspreselection1LostPlusPtCut.CaloIsolationRequirement  = false;
  chiTrackspreselection1LostPlusPtCut.DeDxRequirement           = false;

  class Event chiTrackspreselection1LostPlusIsoAndPtCut("chiTrackspreselection1LostPlusIsoAndPtCut",ofile);
  if(isSignal) chiTrackspreselection1LostPlusIsoAndPtCut.onlyChi      = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.triggerRequirements       = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.trackPreselection         = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.qcdSupression             = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.trackCandidateCutFinal    = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.TrackPtRequirement        = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.NumOfLostOuterRequirement = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.CaloIsolationRequirement  = true;
  chiTrackspreselection1LostPlusIsoAndPtCut.DeDxRequirement           = false;

  class Event chiTrackspreselection1LostPlusIsoAndDeDxCut("chiTrackspreselection1LostPlusIsoAndDeDxCut",ofile);
  if(isSignal) chiTrackspreselection1LostPlusIsoAndDeDxCut.onlyChi      = true;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.triggerRequirements       = true;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.trackPreselection         = true;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.qcdSupression             = true;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.trackCandidateCutFinal    = true;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.TrackPtRequirement        = false;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.NumOfLostOuterRequirement = true;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.CaloIsolationRequirement  = true;
  chiTrackspreselection1LostPlusIsoAndDeDxCut.DeDxRequirement           = true;

  class Event chiTrackspreselection1LostPlusPtAndDeDxCut("chiTrackspreselection1LostPlusPtAndDeDxCut",ofile);
  if(isSignal) chiTrackspreselection1LostPlusPtAndDeDxCut.onlyChi      = true;
  chiTrackspreselection1LostPlusPtAndDeDxCut.triggerRequirements       = true;
  chiTrackspreselection1LostPlusPtAndDeDxCut.trackPreselection         = true;
  chiTrackspreselection1LostPlusPtAndDeDxCut.qcdSupression             = true;
  chiTrackspreselection1LostPlusPtAndDeDxCut.trackCandidateCutFinal    = true;
  chiTrackspreselection1LostPlusPtAndDeDxCut.TrackPtRequirement        = true;
  chiTrackspreselection1LostPlusPtAndDeDxCut.NumOfLostOuterRequirement = true;
  chiTrackspreselection1LostPlusPtAndDeDxCut.CaloIsolationRequirement  = false;
  chiTrackspreselection1LostPlusPtAndDeDxCut.DeDxRequirement           = true;

  class Event chiTracksfullSelectionTrigger("chiTracksfullSelectionTrigger",ofile);
  if(isSignal) chiTracksfullSelectionTrigger.onlyChi   = true;
  chiTracksfullSelectionTrigger.triggerRequirements    = true;
  chiTracksfullSelectionTrigger.trigger             = true;
  chiTracksfullSelectionTrigger.trackPreselection      = true;
  chiTracksfullSelectionTrigger.qcdSupression          = true;
  chiTracksfullSelectionTrigger.trackCandidateCutFinal = true;
  chiTracksfullSelectionTrigger.TrackPtRequirement        = true;
  chiTracksfullSelectionTrigger.NumOfLostOuterRequirement = false;
  chiTracksfullSelectionTrigger.CaloIsolationRequirement  = true;
  chiTracksfullSelectionTrigger.DeDxRequirement           = true;

  class Event chiTracksfullSelection("chiTracksfullSelection",ofile);
  if(isSignal) chiTracksfullSelection.onlyChi   = true;
  chiTracksfullSelection.triggerRequirements    = true;
  chiTracksfullSelection.trackPreselection      = true;
  chiTracksfullSelection.qcdSupression          = true;
  chiTracksfullSelection.trackCandidateCutFinal = true;
  chiTracksfullSelection.TrackPtRequirement        = true;
  chiTracksfullSelection.NumOfLostOuterRequirement = false;
  chiTracksfullSelection.CaloIsolationRequirement  = true;
  chiTracksfullSelection.DeDxRequirement           = true;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // SM model control regions
  class Event chiTracksSMControlCalo("chiTracksSMControlCalo",ofile);
  if(isSignal) chiTracksSMControlCalo.onlyChi   = true;
  chiTracksSMControlCalo.triggerRequirements    = true;
  chiTracksSMControlCalo.trackPreselection      = true;
  chiTracksSMControlCalo.qcdSupression          = true;
  chiTracksSMControlCalo.trackCandidateCutFinal = true;
  chiTracksSMControlCalo.TrackPtRequirement             = false;
  chiTracksSMControlCalo.NumOfLostOuterRequirement      = false;
  chiTracksSMControlCalo.CaloIsolationRequirement       = true;
  chiTracksSMControlCalo.DeDxRequirement                = false;
  chiTracksSMControlCalo.invertCaloIsolationRequirement = true;
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ABCD method with track Pt and DeDx
  /*
  ABCD trackPt_DeDx("_trackPt_DeDx",ofile,isSignal);
 
  trackPt_DeDx.CR1.TrackPtRequirement = true;
  trackPt_DeDx.CR1.DeDxRequirement    = true;
  trackPt_DeDx.CR1.invertTrackPtRequirement = false;
  trackPt_DeDx.CR1.invertDeDxRequirement    = true;

  trackPt_DeDx.CR2.TrackPtRequirement = true;
  trackPt_DeDx.CR2.DeDxRequirement    = true;
  trackPt_DeDx.CR2.invertTrackPtRequirement = true;
  trackPt_DeDx.CR2.invertDeDxRequirement    = true;

  trackPt_DeDx.CR3.TrackPtRequirement = true;
  trackPt_DeDx.CR3.DeDxRequirement    = true;
  trackPt_DeDx.CR3.invertTrackPtRequirement = true;
  trackPt_DeDx.CR3.invertDeDxRequirement    = false;

  trackPt_DeDx.SR.TrackPtRequirement = true;
  trackPt_DeDx.SR.DeDxRequirement    = true;
  trackPt_DeDx.SR.invertTrackPtRequirement = false;
  trackPt_DeDx.SR.invertDeDxRequirement    = false;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ABCD method with track Pt and CaloIsolation
  ABCD trackPt_CaloIso("_trackPt_CaloIso",ofile,isSignal);
 
  trackPt_CaloIso.CR1.TrackPtRequirement       = true;
  trackPt_CaloIso.CR1.CaloIsolationRequirement = true;
  trackPt_CaloIso.CR1.invertTrackPtRequirement       = false;
  trackPt_CaloIso.CR1.invertCaloIsolationRequirement = true;

  trackPt_CaloIso.CR2.TrackPtRequirement       = true;
  trackPt_CaloIso.CR2.CaloIsolationRequirement = true;
  trackPt_CaloIso.CR2.invertTrackPtRequirement       = true;
  trackPt_CaloIso.CR2.invertCaloIsolationRequirement = true;

  trackPt_CaloIso.CR3.TrackPtRequirement       = true;
  trackPt_CaloIso.CR3.CaloIsolationRequirement = true;
  trackPt_CaloIso.CR3.invertTrackPtRequirement       = true;
  trackPt_CaloIso.CR3.invertCaloIsolationRequirement = false;

  trackPt_CaloIso.SR.TrackPtRequirement       = true;
  trackPt_CaloIso.SR.CaloIsolationRequirement = true;
  trackPt_CaloIso.SR.invertTrackPtRequirement       = false;
  trackPt_CaloIso.SR.invertCaloIsolationRequirement = false;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ABCD method with track Pt and DeDx
  ABCD DeDx_CaloIso("_DeDx_CaloIso",ofile,isSignal);
 
  DeDx_CaloIso.CR1.DeDxRequirement          = true;
  DeDx_CaloIso.CR1.CaloIsolationRequirement = true;
  DeDx_CaloIso.CR1.invertDeDxRequirement          = false;
  DeDx_CaloIso.CR1.invertCaloIsolationRequirement = true;

  DeDx_CaloIso.CR2.DeDxRequirement          = true;
  DeDx_CaloIso.CR2.CaloIsolationRequirement = true;
  DeDx_CaloIso.CR2.invertDeDxRequirement          = true;
  DeDx_CaloIso.CR2.invertCaloIsolationRequirement = true;

  DeDx_CaloIso.CR3.DeDxRequirement          = true;
  DeDx_CaloIso.CR3.CaloIsolationRequirement = true;
  DeDx_CaloIso.CR3.invertDeDxRequirement          = true;
  DeDx_CaloIso.CR3.invertCaloIsolationRequirement = false;

  DeDx_CaloIso.SR.DeDxRequirement          = true;
  DeDx_CaloIso.SR.CaloIsolationRequirement = true;
  DeDx_CaloIso.SR.invertDeDxRequirement          = false;
  DeDx_CaloIso.SR.invertCaloIsolationRequirement = false;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ABCD method with track Pt and DeDx
  ABCD trackPt_lostOuterHits("_trackPt_lostOuterHits",ofile,isSignal);
 
  trackPt_lostOuterHits.CR1.TrackPtRequirement        = true;
  trackPt_lostOuterHits.CR1.NumOfLostOuterRequirement = true;
  trackPt_lostOuterHits.CR1.invertTrackPtRequirement        = false;
  trackPt_lostOuterHits.CR1.invertNumOfLostOuterRequirement = true;

  trackPt_lostOuterHits.CR2.TrackPtRequirement        = true;
  trackPt_lostOuterHits.CR2.NumOfLostOuterRequirement = true;
  trackPt_lostOuterHits.CR2.invertTrackPtRequirement        = true;
  trackPt_lostOuterHits.CR2.invertNumOfLostOuterRequirement = true;

  trackPt_lostOuterHits.CR3.TrackPtRequirement        = true;
  trackPt_lostOuterHits.CR3.NumOfLostOuterRequirement = true;
  trackPt_lostOuterHits.CR3.invertTrackPtRequirement        = true;
  trackPt_lostOuterHits.CR3.invertNumOfLostOuterRequirement = false;

  trackPt_lostOuterHits.SR.TrackPtRequirement        = true;
  trackPt_lostOuterHits.SR.NumOfLostOuterRequirement = true;
  trackPt_lostOuterHits.SR.invertTrackPtRequirement        = false;
  trackPt_lostOuterHits.SR.invertNumOfLostOuterRequirement = false;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ABCD method with track Pt and DeDx
  ABCD DeDx_lostOuterHits("_DeDx_lostOuterHits",ofile,isSignal);
 
  DeDx_lostOuterHits.CR1.DeDxRequirement           = true;
  DeDx_lostOuterHits.CR1.NumOfLostOuterRequirement = true;
  DeDx_lostOuterHits.CR1.invertDeDxRequirement           = false;
  DeDx_lostOuterHits.CR1.invertNumOfLostOuterRequirement = true;

  DeDx_lostOuterHits.CR2.DeDxRequirement           = true;
  DeDx_lostOuterHits.CR2.NumOfLostOuterRequirement = true;
  DeDx_lostOuterHits.CR2.invertDeDxRequirement           = true;
  DeDx_lostOuterHits.CR2.invertNumOfLostOuterRequirement = true;

  DeDx_lostOuterHits.CR3.DeDxRequirement           = true;
  DeDx_lostOuterHits.CR3.NumOfLostOuterRequirement = true;
  DeDx_lostOuterHits.CR3.invertDeDxRequirement           = true;
  DeDx_lostOuterHits.CR3.invertNumOfLostOuterRequirement = false;

  DeDx_lostOuterHits.SR.DeDxRequirement           = true;
  DeDx_lostOuterHits.SR.NumOfLostOuterRequirement = true;
  DeDx_lostOuterHits.SR.invertDeDxRequirement           = false;
  DeDx_lostOuterHits.SR.invertNumOfLostOuterRequirement = false;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // ABCD method with track Pt and DeDx
  ABCD CaloIso_lostOuterHits("_CaloIso_lostOuterHits",ofile,isSignal);
 
  CaloIso_lostOuterHits.CR1.CaloIsolationRequirement  = true;
  CaloIso_lostOuterHits.CR1.NumOfLostOuterRequirement = true;
  CaloIso_lostOuterHits.CR1.invertCaloIsolationRequirement  = false;
  CaloIso_lostOuterHits.CR1.invertNumOfLostOuterRequirement = true;

  CaloIso_lostOuterHits.CR2.CaloIsolationRequirement  = true;
  CaloIso_lostOuterHits.CR2.NumOfLostOuterRequirement = true;
  CaloIso_lostOuterHits.CR2.invertCaloIsolationRequirement  = true;
  CaloIso_lostOuterHits.CR2.invertNumOfLostOuterRequirement = true;

  CaloIso_lostOuterHits.CR3.CaloIsolationRequirement  = true;
  CaloIso_lostOuterHits.CR3.NumOfLostOuterRequirement = true;
  CaloIso_lostOuterHits.CR3.invertCaloIsolationRequirement  = true;
  CaloIso_lostOuterHits.CR3.invertNumOfLostOuterRequirement = false;

  CaloIso_lostOuterHits.SR.CaloIsolationRequirement  = true;
  CaloIso_lostOuterHits.SR.NumOfLostOuterRequirement = true;
  CaloIso_lostOuterHits.SR.invertCaloIsolationRequirement  = false;
  CaloIso_lostOuterHits.SR.invertNumOfLostOuterRequirement = false;
  */
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //-------------------------------------------------------------------------------------
  // Read file for bad ecal cells and csc chambers 
  //-------------------------------------------------------------------------------------
  ifstream inputFile("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/BadCSCChambers.txt");
  int lines;
  int i=0;
  double n,m;
  while(inputFile>>n>>m){
    etaCSC.push_back(n);
    phiCSC.push_back(m);
    i++;
  }
  lines=i;
  inputFile.close();

  inputFile.open("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/DeadEcalChannelsNew.txt");
  i=0;
  while(inputFile>>n>>m){
    etaEcal.push_back(n);
    phiEcal.push_back(m);
    i++;
  }
  lines=i;
  
  //-------------------------------------------------------------------------------------
  // Read discriminator templates
  //-------------------------------------------------------------------------------------
  if(isData){
    template_pixel = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Discrim_template_pixel_data_2012.root");
    template_strip = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Data7TeV_Deco_SiStripDeDxMip_3D_Rcd.root");
  }
  else{
    template_pixel = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Discrim_template_pixel_mc_2012.root");
    template_strip = loadDeDxTemplate("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/Discrim_Templates_MC_2012.root");
  }

  //-------------------------------------------------------------------------------------
  // Declare additional branch addresses for hit information
  //-------------------------------------------------------------------------------------
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsDeDx",&HitsDeDx);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsPathlength",&HitsPathlength);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsShapetest",&HitsShapetest);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsSubdetId",&HitsSubdetid);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsEta",&HitsEta);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsPhi",&HitsPhi);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsTransverse",&HitsTransverse);
  
  //-------------------------------------------------------------------------------------
  // Loop over events
  //-------------------------------------------------------------------------------------

  clock_t start, stop;
  double time = 0.0;

  assert((start = clock())!=-1);
  // ***********************************************************************************************************************
  // ----------  Stuff for PU reweighing -----------------------------------------------------------------------------------
  //TFile file_PUdata("/nfs/dust/cms/user/rathjd/VBF-LS-tau/PU/DataPUFile_22Jan2013ReReco_Run2012.root", "read");
  //TFile file_PUmc("/afs/desy.de/user/t/tlenz/HSCPworkdir/PUhistos/TrueNumInteractions_0.root", "read");
  TFile *file_PUdata = new TFile("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/pileup/MyDataPileupHistogram.root", "read");
  TFile *file_PUmc   = new TFile("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/pileup/PileUpHistoCycle.MC.WJets.root", "read");

  if(PUunc){
    delete file_PUdata;
    if(up)   file_PUdata = new TFile("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/pileup/MyDataPileupHistogram_729.root", "read");
    if(down) file_PUdata = new TFile("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/pileup/MyDataPileupHistogram_659.root", "read");
  }


  //TH1F *PUweights    =   (TH1F*)file_PUdata.Get("analyzeHiMassTau/NVertices_0");
  cout<<"filename = "<<file_PUdata->GetName()<<endl;
  TH1F *PUweights    =   (TH1F*)file_PUdata->Get("pileup");
  PUweights          ->  Scale(1/PUweights->Integral());
  //TH1F *PUmc         =   (TH1F*)file_PUmc.Get("hTrueNumInteractions_0");
  TH1F *PUmc         =   (TH1F*)file_PUmc->Get("N_pileup_hist");
  PUmc               ->  Scale(1/PUmc->Integral());
  
  PUweights->Divide(PUmc);
 
  weight            = 1.;
  weightReweighting = 1.;
  // -----------------------------------------------------------------------------------------------------------------------
  // For Trigger Efficiency Uncertainty
  TH1D* hTrigger = 0;
  if(TrigEffunc){
    TFile* fileTrigger = new TFile("/afs/desy.de/user/t/tlenz/HighDeDx-DisappTrks-Analyzer/data/TriggerMetSF.root", "read");
    fileTrigger -> GetObject("effDiff",hTrigger);
  }
  // -----------------------------------------------------------------------------------------------------------------------
  // ***********************************************************************************************************************

  cout<<endl<<"Number Of Events = "<<nevents<<endl<<endl;
  
  randGenerator = new TRandom3(0);
  
  for(int entry=0; entry <nevents; ++entry)
    {

      // Read event into memory
      stream.read(entry);
      stream._chain->GetEntry(entry);

      // NB: call to clear object selection map (indexmap)
      initialize();
	  
      // Uncomment the following line if you wish to copy variables into
      // structs. See the header file analyzer.h to find out what structs
      // are available. Alternatively, you can call individual fill functions.
      fillObjects();
      
      /******************************************************************************************************************************
       ******************************************************************************************************************************
       ******************************************************************************************************************************
       *****************************************************************************************************************************/
      weight=1.;
      weightReweighting=1;
      findChiInGenParticleCollection();
      findChi0InGenParticleCollection();
      findChiInSimTrackCollection();
      findChiDecayVertex();
      
      
      // Chargino event reweighting:
      if(isSignal && TargetLifetime !=0){
	bool reweighted = true;
	CurrentLifetime = getCurrentLifetime(stream.filename());
	//cout<<"CurrentLifetime = "<<CurrentLifetime<<endl;

	// 1.) Get both charginos from genParticle collection and their proper lifetime

	if(CurrentLifetime != TargetLifetime){

	  if (ChiTrack.size() > 2){
	    cout << "Too many charginos!: " << ChiTrack.size() << endl;
	    continue;
	  }
	  if (ChiTrack.size() < 2 && !isSignalC1N1){
	    cout << "Too few charginos!: " << ChiTrack.size() << endl;
	    continue;
	  }
	  if (ChiTrack.size() < 1 && isSignalC1N1){
	    cout << "Too few charginos!: " << ChiTrack.size() << endl;
	    continue;
	  }

	  for(unsigned int i=0; i<ChiTrack.size();i++){

 	    if(ChiTrack[i].genenergy == 0){
	      cout<<"no matched genParticle found = skip event!"<<endl;
	      reweighted=false;
	      break;
	    }

	    double beta  = ChiTrack[i].genp/ChiTrack[i].genenergy;
	    double gamma = 1./std::sqrt(1.-pow(beta,2));

	    double rho      = TMath::Sqrt( TMath::Power(ChiTrack[i].SimVertexposition_x - Vertex[0].x,2) + TMath::Power(ChiTrack[i].SimVertexposition_y - Vertex[0].y,2));
	    double z        = TMath::Abs(ChiTrack[i].SimVertexposition_z - Vertex[0].z);
	    double distance = TMath::Sqrt( TMath::Power(rho,2) + TMath::Power(z,2) );
	    double ProperLifetime = distance/(beta*gamma);

	    double wtTarget       = (1. / TargetLifetime)  * TMath::Exp(-(ProperLifetime) /  TargetLifetime  );  
	    double wtCurrent      = (1. / CurrentLifetime) * TMath::Exp(-(ProperLifetime) /  CurrentLifetime );
	    double wt = wtTarget / wtCurrent;


	    if(ProperLifetime<0){
	      cout<<"Warning:  Found event with ctau<0."<<endl;
	      wt=0.;
	    }

 	    if(ChiTrack[i].SimVertexposition_z==0 && ChiTrack[i].SimVertexposition_y==0 && ChiTrack[i].SimVertexposition_x==0) wt=0;

	    weightReweighting *= wt;
	  }
	  if(!reweighted){
	    cout<<"reweighting not possible - skip event"<<endl;
	    continue;
	  }
	}
      }
      /******************************************************************************************************************************/
      if(isSignal){
	if (ChiTrack.size() != 2 && !isSignalC1N1){
	  cout << "Too many/few charginos!: " << ChiTrack.size() << endl;
	  return 1;
	}
	if (ChiTrack.size() != 1 && isSignalC1N1){
	  cout << "Too many/few charginos!: " << ChiTrack.size() << endl;
	  return 1;
	}

	TLorentzVector chi1 = lorentzVector(ChiTrack[0].genpt, ChiTrack[0].geneta, ChiTrack[0].genphi, ChiTrack[0].genenergy);
	TLorentzVector chi2;
	if(isSignalC1N1) chi2 = lorentzVector(chi0GenParticle[0].pt, chi0GenParticle[0].eta, chi0GenParticle[0].phi, chi0GenParticle[0].energy);
	else             chi2 = lorentzVector(ChiTrack[1].genpt, ChiTrack[1].geneta, ChiTrack[1].genphi, ChiTrack[1].genenergy);
	double ptSystem = (chi1+chi2).Pt();
   
	if(!ISRunc){
	  if(ptSystem>=0 && ptSystem<=120.)        weight *= 1.00;
	  else if(ptSystem>120 && ptSystem<=150.)  weight *= 0.95;
	  else if(ptSystem>150 && ptSystem<=250.)  weight *= 0.90;
	  else if(ptSystem>250)                    weight *= 0.85;
	}
	
	if(ISRunc){
	  if(down){
	    if(ptSystem>=0 && ptSystem<=120.)        weight *= 1.00;
	    else if(ptSystem>120 && ptSystem<=150.)  weight *= 0.90;
	    else if(ptSystem>150 && ptSystem<=250.)  weight *= 0.80;
	    else if(ptSystem>250)                    weight *= 0.60;
	  }
	  else if(up)                                weight *= 1.00;
	}
      }
      /******************************************************************************************************************************/
      if(isSignal && JECunc){

	for(unsigned int i=0; i<evt::Jet.size(); i++){

	  // correct jetPt
	  if(up)           Jet[i].pt = Jet[i].pt*(1 + Jet[i].jecUnc);
	  else if(down)    Jet[i].pt = Jet[i].pt*(1 - Jet[i].jecUnc);
	}
      }
      /******************************************************************************************************************************/
      if(isSignal && TrigEffunc){

	// 1.) Calculate MET
	double met = sqrt(pow(MET_pt,2) + pow(MET_pz,2));
	int bin = hTrigger->FindBin(met);
	weight *= (hTrigger->GetBinContent(bin) + 1) ;

      }
      //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(isSignal){

	double jerCentral[7]     = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
	double jerDown[7]        = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};
	double jerUp[7]          = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};

	double *jer = jerCentral;

	if(JERunc){

	  if(up)        jer = jerUp;
	  else if(down) jer = jerDown;

	}

	// 1.) Do matching to genJets
	for(unsigned int i=0; i<evt::Jet.size(); i++){

	  int c = 0;
	  if(abs(Jet[i].eta) >= 0.0 && abs(Jet[i].eta) < 0.5)          c = 0;
	  else if(abs(Jet[i].eta) >= 0.5 && abs(Jet[i].eta) < 1.1)     c = 1;
	  else if(abs(Jet[i].eta) >= 1.1 && abs(Jet[i].eta) < 1.7)     c = 2;
	  else if(abs(Jet[i].eta) >= 1.7 && abs(Jet[i].eta) < 2.3)     c = 3;
	  else if(abs(Jet[i].eta) >= 2.3 && abs(Jet[i].eta) < 2.8)     c = 4;
	  else if(abs(Jet[i].eta) >= 2.8 && abs(Jet[i].eta) < 3.2)     c = 5;
	  else if(abs(Jet[i].eta) >= 3.2 && abs(Jet[i].eta) < 5.0)     c = 6;
	  else c = 1;

	  for(unsigned int j=0; j<evt::GenJet.size(); j++){

	    double dPhi = std::abs(TVector2::Phi_mpi_pi(Jet[i].phi - GenJet[j].phi));
	    double dEta = std::abs(Jet[i].eta - GenJet[j].eta);
	    double dR   = std::sqrt(dPhi*dPhi + dEta*dEta);

	    if(dR<0.15){
	      Jet[i].pt = max(0.0, GenJet[j].pt + jer[c]*(Jet[i].pt - GenJet[j].pt) );
	    }
	  }
	}
      }
      /******************************************************************************************************************************
       ******************************************************************************************************************************
       ******************************************************************************************************************************
       *****************************************************************************************************************************/

      if(!edmEventHelper_isRealData){
	weight*=PUweights->GetBinContent(PUweights->FindBin(PileupSummaryInfo_getTrueNumInteractions[0]));
      }
      //weight =1.;
      ofile.count("NoCuts", weight);


      //------------- Calculate track dependent variables from hits and save them in trk coll ----
      for(unsigned int i=0; i<evt::Track.size();i++){
	
	double ASmiOnTheFly            = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], edmEventHelper_isRealData, template_strip, template_pixel,1,0); 
	double ASmiNPOnTheFly          = dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], edmEventHelper_isRealData, template_strip, template_pixel,0); 
	double ASmiOnTheFly_3          = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,3);
	double ASmiNPOnTheFly_3        = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,0,3);
	double ASmiOnTheFly_7          = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,7);
	double ASmiNPOnTheFly_7        = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,0,7);
	double ASmiOnTheFly_woLastHit  = -1;//dEdxOnTheFly(&(*HitsDeDx)[i], &(*HitsShapetest)[i], &(*HitsPathlength)[i], &(*HitsSubdetid)[i], &(*HitsTransverse)[i], 1, template_strip, template_pixel,1,-1);

	Track[i].ASmi           = ASmiOnTheFly;
       	Track[i].ASmiNP         = ASmiNPOnTheFly;
	Track[i].ASmi_3         = ASmiOnTheFly_3;
       	Track[i].ASmiNP_3       = ASmiNPOnTheFly_3;
	Track[i].ASmi_7         = ASmiOnTheFly_7;
       	Track[i].ASmiNP_7       = ASmiNPOnTheFly_7;
	Track[i].ASmi_woLastHit = ASmiOnTheFly_woLastHit;

      }
      
      
      //-------------------------------------------------------------- Cuts ---------------------------------------------------------
      nVertices_0->Fill(nVertex);
      hTrueNumInteractions_0->Fill(PileupSummaryInfo[0].getTrueNumInteractions);
      hPU_NumInteractions_0->Fill(PileupSummaryInfo[0].getPU_NumInteractions);

      //noSelection.Selection();
      //triggerRequirements.Selection();
      //preselection.Selection();
      //fullSelection.Selection();
      //chiTracksnoSelection.Selection();
      
      if(PDFunc && isSignal){
	// 2.) Run over the preselection and store pdf weights
	for (unsigned int i = 0; i < nbins_pdfweight; i++){ 
	  double PDFweight = (pdfWeights_cteq66[i]/pdfWeights_cteq66[0]) / pdf_weight_sum[i] * nevents;
	  weight *= PDFweight;
	  PDFselection[i].Selection();
	}
      } 
       
      chiTrackstriggerRequirements.Selection();
      //chiTrackstriggerRequirementsTrigger.Selection();  //Not Needed
      chiTracksQCDsupression.Selection();
      chiTracksQCDsupressionTrigger.Selection();
      chiTrackspreselection.Selection();
      chiTrackspreselectionTrigger.Selection();
      chiTrackspreselectionNoQCDCuts.Selection();
      chiTrackspreselectionTriggerNoQCDCuts.Selection();
      chiTrackspreselectionNoTrigger.Selection();
      chiTrackspreselectionNoQCDCutsNoTrigger.Selection();
      chiTracksfullSelection.Selection();
      chiTracksfullSelectionTrigger.Selection();
     

      //chiTrackspreselectionWjets.Selection();
      //chiTracksSMControlCalo.Selection();

      
    }//end of loop over events
 
  stop = clock();
  time = (double) (stop-start)/CLOCKS_PER_SEC;
  cout<<endl<<endl<<"time = "<<time/60.<<endl;

  stream.close();
  ofile.close();
  
  cout<<endl;
  cout<<"nChiInSimVertex = "<<nChiInSimVertex<<endl;
  cout<<"nChiInSimTrack  = "<<nChiInSimTrack<<endl<<endl;

  return 0;
}
