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
  bool isSignal = false;
  vector<string> filenames = getFilenames(cmdline.filelist);
  if(filenames[0].find("pMSSM") != std::string::npos) isSignal = true;
  itreestream stream(filenames, "Events");
  if ( !stream.good() ) error("unable to open ntuple file(s)");

  // Get number of events to be read

  int nevents = stream.size();
  //cout << "Number of events: " << nevents << endl;

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
  // Declaration of Variables
  //-------------------------------------------------------------------------------------
  cout<<endl<<endl<<"------------ Is it signal : "<<isSignal<<" --------------"<<endl<<endl;
  class Event noSelection("noSelection",ofile);
  class Event triggerRequirements("triggerRequirements",ofile);
  triggerRequirements.preselection        = true;
  triggerRequirements.triggerRequirements = true;
  class Event preselection("preselection",ofile);
  preselection.preselection               = true;
  preselection.triggerRequirements        = true;
  preselection.trackPreselection          = true;
  preselection.qcdSupression              = true;
  class Event fullSelection("fullSelection",ofile);
  fullSelection.preselection              = true;
  fullSelection.triggerRequirements       = true;
  fullSelection.trackPreselection         = true;
  fullSelection.qcdSupression             = true;
  fullSelection.trackCandidateCut         = true;
  class Event chiTracks("chiTracks",ofile);
  if(isSignal) chiTracks.onlyChi               = true;
  class Event chiTrackstriggerRequirements("chiTrackstriggerRequirements",ofile);
  if(isSignal) chiTrackstriggerRequirements.onlyChi               = true;
  chiTrackstriggerRequirements.preselection        = true;
  chiTrackstriggerRequirements.triggerRequirements = true;
  class Event chiTrackspreselection("chiTrackspreselection",ofile);
  if(isSignal) chiTrackspreselection.onlyChi             = true;
  chiTrackspreselection.preselection        = true;
  chiTrackspreselection.triggerRequirements = true;
  chiTrackspreselection.trackPreselection   = true;
  chiTrackspreselection.qcdSupression       = true;
  class Event chiTracksfullSelection("chiTracksfullSelection",ofile);
  if(isSignal) chiTracksfullSelection.onlyChi             = true;
  chiTracksfullSelection.preselection        = true;
  chiTracksfullSelection.triggerRequirements = true;
  chiTracksfullSelection.trackPreselection   = true;
  chiTracksfullSelection.qcdSupression       = true;
  chiTracksfullSelection.trackCandidateCut   = true;

  // Control region 1 small track pt and small dEdx
  class Event CR1("CR1",ofile);
  CR1.preselection              = true;
  CR1.triggerRequirements       = true;
  CR1.trackPreselection         = true;
  CR1.qcdSupression             = true;
  CR1.trackCandidateCut         = true;
  CR1.invertTrackPtRequirement         = false;
  CR1.invertNumOfValidHitsRequirement  = false;
  CR1.invertCaloIsolationRequirement   = false;
  CR1.invertNumOfLostOuterRequirement  = false;
  CR1.invertDeDxRequirement            = true;
  class Event chiTracksCR1("chiTracksCR1",ofile);
  if(isSignal) chiTracksCR1.onlyChi = true;
  chiTracksCR1.preselection         = true;
  chiTracksCR1.triggerRequirements  = true;
  chiTracksCR1.trackPreselection    = true;
  chiTracksCR1.qcdSupression        = true;
  chiTracksCR1.trackCandidateCut    = true;
  chiTracksCR1.invertTrackPtRequirement         = false;
  chiTracksCR1.invertNumOfValidHitsRequirement  = false;
  chiTracksCR1.invertCaloIsolationRequirement   = false;
  chiTracksCR1.invertNumOfLostOuterRequirement  = false;
  chiTracksCR1.invertDeDxRequirement            = true;
  

  //-------------------------------------------------------------------------------------
  // Read file for bad ecal cells and csc chambers 
  //-------------------------------------------------------------------------------------
  ifstream inputFile("data/BadCSCChambers.txt");
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

  inputFile.open("data/DeadEcalChannelsNew.txt");
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
  TH3 * template_pixel = loadDeDxTemplate("data/Discrim_template_pixel_data_2012.root");
  TH3 * template_strip = loadDeDxTemplate("data/Data7TeV_Deco_SiStripDeDxMip_3D_Rcd.root");

  //-------------------------------------------------------------------------------------
  // Declare additional branch addresses for hit information
  //-------------------------------------------------------------------------------------
  
  std::vector<std::vector<double> > *HitsDeDx = 0;
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsDeDx",&HitsDeDx);
  double electronenergy[1000] = {0};
  stream._chain->SetBranchAddress("patElectron_patElectronsLoosePFlow.energy",electronenergy);
  std::vector<std::vector<double> > *HitsPathlength = 0;
  std::vector<std::vector<int> >    *HitsShapetest  = 0;
  std::vector<std::vector<int> >    *HitsSubdetid   = 0;
  std::vector<std::vector<double> > *HitsEta        = 0;
  std::vector<std::vector<double> > *HitsPhi        = 0;
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsPathlength",&HitsPathlength);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsShapetest",&HitsShapetest);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsSubdetId",&HitsSubdetid);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsEta",&HitsEta);
  stream._chain->SetBranchAddress("recoTrackHelper_TrackRefitter_HitsPhi",&HitsPhi);
  

  //-------------------------------------------------------------------------------------
  // Loop over events
  //-------------------------------------------------------------------------------------

  clock_t start, stop;
  double time = 0.0;

  assert((start = clock())!=-1);
  // ***********************************************************************************************************************
  // ----------  Stuff for PU reweighing -----------------------------------------------------------------------------------
  TFile file_PUdata("/nfs/dust/cms/user/rathjd/VBF-LS-tau/PU/DataPUFile_22Jan2013ReReco_Run2012.root", "read");
  TFile file_PUmc("/afs/desy.de/user/t/tlenz/HSCPworkdir/PUhistos/TrueNumInteractions_0.root", "read");
  //TFile file_PUmc("/nfs/dust/cms/user/rathjd/VBF-LS-tau/PU/S10MC_PUFile.root", "read");
 
  //TH1F *PUweights = (TH1F*)file_PU.Get("ratio");
  TH1F *PUweights = (TH1F*)file_PUdata.Get("analyzeHiMassTau/NVertices_0");
  PUweights->Scale(1/PUweights->Integral());
  //cout<<"PUweights = "<<PUweights->Integral()<<endl;
  //  TH1F *PUmc = (TH1F*)file_PUmc.Get("analyzeHiMassTau/NVertices_0");
  TH1F *PUmc = (TH1F*)file_PUmc.Get("hTrueNumInteractions_0");
  PUmc->Scale(1/PUmc->Integral());
  //cout<<"PUmc = "<<PUmc->Integral()<<endl;
  
  PUweights->Divide(PUmc);
  //cout<<"PUweights = "<<PUweights->Integral()<<endl;
 
  weight = 1.;
  // -----------------------------------------------------------------------------------------------------------------------
  // ***********************************************************************************************************************

  cout<<endl<<"Number Of Events = "<<nevents<<endl<<endl;
 
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
      
      if(!edmEventHelper_isRealData){
	weight=PUweights->GetBinContent(PUweights->FindBin(PileupSummaryInfo_getTrueNumInteractions[0]));
	if(weight<0.0001) cout<<"weight = "<<weight<<endl;
      }
      //weight =1.;
      ofile.count("NoCuts", weight);
      
      //-------------------------------------------------------------- Cuts ---------------------------------------------------------
      findChiInGenParticleCollection();
      nVertices_0->Fill(nVertex);
      hTrueNumInteractions_0->Fill(PileupSummaryInfo[0].getTrueNumInteractions);
      hPU_NumInteractions_0->Fill(PileupSummaryInfo[0].getPU_NumInteractions);
      
  
      noSelection.Selection();
      triggerRequirements.Selection();
      preselection.Selection();
      fullSelection.Selection();
      chiTracks.Selection();
      chiTrackstriggerRequirements.Selection();
      chiTrackspreselection.Selection();
      chiTracksfullSelection.Selection();
      CR1.Selection();
      chiTracksCR1.Selection();
            
    }//end of loop over events
 
  stop = clock();
  time = (double) (stop-start)/CLOCKS_PER_SEC;
  cout<<endl<<endl<<"time = "<<time/60.<<endl;

  stream.close();
  ofile.close();
  
  cout<<endl;
  return 0;
}
