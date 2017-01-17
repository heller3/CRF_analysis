#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TLatex.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TGaxis.h"

#include <string>

#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;
TString inputdir = "/afs/cern.ch/work/r/rheller/CMSSW_8_0_20/src/HCALPFG/HcalTupleMaker/output/";
bool make_plots=true;
bool collisions=false;
const int nieta = 20;
const int niphi = 6;
const int nsections=20;
TString runnum;
int peakTS[nieta][niphi];

void plot_pulse(TH1F* h1_vec[],TH1F* h1_ped[], TString type);
void plot_vs_section(TH1F * h1_vec[][nsections], TString type);
void plot_distribution(TH1F* h1_ieta[], int ieta, TString type);
void plot_distribution(TH1F* h1_ieta[], TString type, TString name);


//find energy and time using peak time slice, n-1, n+1 and n+2
float energy(vector<float> fC_subt){
  int max_idx = distance(fC_subt.begin(),max_element(fC_subt.begin(),fC_subt.end()));
  //cout<<"Max index is "<<max_idx<<endl;
  int min = max_idx-1;
  int max = max_idx+2;
  if(min<0) min=0;
  if(max>9) max=9; 
  float energy=0.;
  for(int i=min;i<=max; i++){
    //cout<<"i is "<<i<<", fC is "<<fC_subt.at(i)<<endl;
    energy+=fC_subt.at(i);
  }
  // cout<<"total is "<<energy<<endl;
  return energy;
 
}


//Find energy using a constant window starting from TS2
float energy(vector<float> fC_subt, int window){
 
  int min = 2;
  int max = min+window;
  float energy=0.;
  for(int i=min;i<max; i++){
    energy+=fC_subt.at(i);
  }
  return energy;
}

//Find energy using a constant window starting from arbitrary TS
float energy(vector<float> fC_subt, int start, int length){
 
  int min = start;
  int max = min+length;
  float energy=0.;
  if(min<0) min=0;
  if(max>9) max=9; 
  for(int i=min;i<max; i++){
    energy+=fC_subt.at(i);
  }
  return energy;
}

//find energy and time given peak time slice to center
float time(vector<float> fC_subt, int start, int length){
  int min = start;
  int max =min+length;
  if(min<0) min=0;
  if(max>9) max=9; 
  float time=0.;
  float energy=0.;
  for(int i=min;i<=max; i++){
    time+= (i+1)*fC_subt.at(i); // offset index by 1 to avoid problems with i==0, subtract by one later
    energy+=fC_subt.at(i);
  }
  time = time/energy;
  time -= 1.;
  return time;
}


//find energy and time using peak time slice, n-1, n+1 and n+2
float time(vector<float> fC_subt){
  int max_idx = distance(fC_subt.begin(),max_element(fC_subt.begin(),fC_subt.end()));
  int min = max_idx-1;
  int max = max_idx+2;
  if(min<0) min=0;
  if(max>9) max=9; 
  float time=0.;
  float energy=0.;
  for(int i=min;i<=max; i++){
    time+= (i+1)*fC_subt.at(i); // offset index by 1 to avoid problems with i==0, subtract by one later
    energy+=fC_subt.at(i);
  }
  time = time/energy;
  time -= 1.;
  return time;
}

//tdc_time() is not debugged yet
float tdc_time(vector<int> tdc){
  
  int first=-1;
  int val=0;
  for(unsigned int i=0;i<tdc.size();i++){
    int thistdc = tdc.at(i);

    cout<<"this tdc, "<<i<<" "<<thistdc<<endl;
    if(thistdc == 63) continue;
    
    // Flag nonsense: direct transition from 63 to 62, or nonsense value
    if(first==-1){
      if(thistdc ==62) {first=11; break;}
      if(thistdc>=50 && thistdc!=62 && thistdc !=63) {first=12; break;}
      else if (thistdc<=49){
	first=i;
	val=thistdc;
      }
    }

    //Flag events with odd tdc results after peak
    if(first!=-1){
      if(thistdc!=62 && thistdc!=63) {first=13; break;}
      //This is the problem
    } 
  }
   
  float time;
  if(first==-1) time=14;
  else time= first+(val/50.);
  cout<<"first is "<<first<<endl;
  cout<<"time is "<<time<<endl;
  return time;
    
}





void make_hists(TString runnumber, bool calib=false, bool col=false){
  cout<<"Processing run "<<runnum<<endl;
  runnum=runnumber;
  collisions=col;

  if(collisions && calib){ cout<< "ERROR calib and collisions turned on "<<endl; return;}
  
  if(collisions)  TGaxis::SetMaxDigits(6);
  int run_int = runnumber.Atoi();
  gStyle->SetGridColor(16);

  //get list of all pedestal tables, use pedestal run with closest run number 
  TString peddir="pedestals/";
  TSystemDirectory pdir(peddir, peddir);
  TList *files = pdir.GetListOfFiles();
  TString pedname;
  if (files) {
    TSystemFile *file;
    TString fname;
    TIter next(files);
    int min_dist = 1000000;
    while ((file=(TSystemFile*)next())) {
      fname = file->GetName();
      fname.ReplaceAll("ped_","");
      fname.ReplaceAll(".txt","");
      int pedrun = fname.Atoi();
      if(fabs(pedrun-run_int) < min_dist){
	min_dist = fabs(pedrun-run_int);
	pedname = Form("%sped_%i.txt",peddir.Data(),pedrun);
      }
    }
  }


  cout<<"Using pedestal table "<< pedname<<endl;
  std::ifstream inputFile(pedname);
  std::string line;
  float peds[nieta][niphi][4];

   while(getline(inputFile, line)) {
      if (!line.length()) continue;
      std::istringstream iss(line);
      int eta,phi,dep;
      string det;
      float cap1,cap2,cap3,cap4;
      iss>>eta>>phi>>dep>>det>>cap1>>cap2>>cap3>>cap4;
      
      peds[eta][phi][0]=cap1;
      peds[eta][phi][1]=cap2;
      peds[eta][phi][2]=cap3;
      peds[eta][phi][3]=cap4;
      
   } 

   
  TString rootfile = inputdir+"outputFile_USC_"+runnum+".root";
  if(collisions) rootfile="/afs/cern.ch/work/r/rheller/test/CMSSW_8_0_20/src/HCALPFG/HcalTupleMaker/USC/HcalNZS_"+runnum+".root";
		   // "/afs/cern.ch/work/r/rheller/CMSSW_8_0_20/src/HCALPFG/HcalTupleMaker/ZeroBias/Run"+runnum+".root";
  //inputdir+"outputFile_collisions_"+runnum+".root";
  if(calib) rootfile = inputdir+"outputFile_calib_"+runnum+".root";
  cout<<"Opening "+rootfile<<endl;
  TFile *f = TFile::Open(rootfile, "READ");
  TDirectory* dir = f->GetDirectory("hcalTupleTree");
  dir->cd();
  TTree *tree = (TTree*)dir->Get("tree");
  int nevents = tree->GetEntries();
  
  cout<< Form("Number of events: %i",nevents)<<endl;


  // event,ls and run  
  UInt_t   event = 0;
  tree->SetBranchAddress("event", &event);
  UInt_t   ls = 0;
  tree->SetBranchAddress("ls", &ls);
  UInt_t   run = 0;
  tree->SetBranchAddress("run", &run);
  UInt_t   bx = 0;
  tree->SetBranchAddress("bx", &bx);

  // QIE11 

  vector<int>   *QIE11DigiSubdet = 0;
  tree->SetBranchAddress("QIE11DigiSubdet", &QIE11DigiSubdet);
  vector<int>   *QIE11DigiIEta = 0;
  tree->SetBranchAddress("QIE11DigiIEta", &QIE11DigiIEta);
  vector<int>   *QIE11DigiIPhi = 0;
  tree->SetBranchAddress("QIE11DigiIPhi", &QIE11DigiIPhi);
  vector<vector<int> >   *QIE11DigiCapID = 0;
  tree->SetBranchAddress("QIE11DigiCapID", &QIE11DigiCapID);
  vector<vector<float> >   *QIE11DigiFC = 0; 
  tree->SetBranchAddress("QIE11DigiFC", &QIE11DigiFC);
  // vector<vector<int> >   *QIE11DigiTDC = 0; 
  //tree->SetBranchAddress("QIE11DigiTDC", &QIE11DigiTDC);
  Int_t   laserType = 0;
  if(calib){
    tree->SetBranchAddress("laserType", &laserType);
  }

  TString rfilename= "hists/hists_"+runnum+".root";
  if(collisions)rfilename= "hists/hists_collisions_"+runnum+".root";
  TFile* file = new TFile(rfilename,"RECREATE");
  
  //Pulse shapes
  TH1F *h1_fC[nieta][niphi];  //Hists to store average fC in each time slice, before pedestal subtraction
  TH1F *h1_ped[nieta][niphi];  //Hists to store average pedestals
  TH1F *h1_fC_ped_subt[nieta][niphi]; //Hists to store average fC in each time slice, after pedestal subtraction


  //energy distributions (charge distributions)
  TH1F *h1_energy[nieta][niphi]; // Uses pedestal subtracted values and a window of 4 time slices around peak TS (for laser runs)
  TH1F *h1_energy_no_ps[nieta][niphi]; // No pedestal subtraction. Uses a window of 4 time slices around peak TS (for laser runs)

  TH1F *h1_energy_const[nieta][niphi]; //Filled with constant window, sum of TS 2-9. No pedestal subtraction (for plotting pedestal drift)
  TH1F *h1_energy_const_zoom[nieta][niphi]; //Same as above, with different binning
  TH1F *h1_energy_const4[nieta][niphi]; //Same as h1_energy_const, but uses sum of TS 2-5
  

  TH1F *h1_energy_laserconst[nieta][niphi]; //Window of 4, using constant window adjusted automatically run by run (and channel by channel) to match the laser arrival time (rather than event by event)

  TH1F *h1_energy_sections[nieta][niphi][nsections]; //this is filled with the same values as h1_energy_const, but separated into different nsections different periods throughout the run. This is to look for time dependence within one run.

  TH1F *h1_timing[nieta][niphi]; //Filled with charge-weighted average (uses peak finding)
  TH1F *h1_timing_laserconst[nieta][niphi]; //Filled with charge-weighted average, using run-by-run constant window
  
  //The following two are TDC info not used right now
  TH1F *h1_tdc_timing[nieta][niphi]; 
  TH1F *h1_delta_timing[nieta][niphi]; 




  //Range for charge distributions
  int enbins=1000;
  int eminx=0;
  int emaxx=10000;
  
  if(collisions){ emaxx*=20; enbins=500;}

  for(int ieta=0; ieta<nieta; ieta++) 
    { 
      for(int iphi=0; iphi<niphi; iphi++) 
        {
	  for(int isec=0; isec<nsections; isec++) {
	    h1_energy_sections[ieta][iphi][isec]= new TH1F(Form("h1_energy_ieta%i_iphi%i_sec%i",ieta,iphi,isec),Form("h1_energy_ieta%i_iphi%i_sec%i",ieta,iphi,isec),enbins,eminx,emaxx);
	    h1_energy_sections[ieta][iphi][isec]->Sumw2();
	    h1_energy_sections[ieta][iphi][isec]->StatOverflows(kTRUE);
	    // h1_energy_sections[ieta][iphi][isec]->SetTitle(Form("Run %s, Fiber %i, Channel %i, section %i; Total charge [fC]; Events",runnum.Data(),ieta+2,iphi,isec));
	  }
	  
	 
	  h1_fC[ieta][iphi]= new TH1F(Form("h1_fC_ieta%i_iphi%i",ieta,iphi),Form("h1_fC_ieta%i_iphi%i",ieta,iphi), 10, -0.5, 9.5);
	  h1_fC[ieta][iphi]->Sumw2();	
	  h1_fC[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; TS; Average charge [fC]",runnum.Data(),ieta+2,iphi));

	  h1_ped[ieta][iphi]= new TH1F(Form("h1_ped_ieta%i_iphi%i",ieta,iphi),Form("h1_ped_ieta%i_iphi%i",ieta,iphi), 10, -0.5, 9.5);
	  h1_ped[ieta][iphi]->Sumw2();

	  h1_fC_ped_subt[ieta][iphi]= new TH1F(Form("h1_fC_ieta%i_iphi%i_ped_subt",ieta,iphi),Form("h1_fC_ieta%i_iphi%i_ped_subt",ieta,iphi), 10, -0.5, 9.5);
	  h1_fC_ped_subt[ieta][iphi]->Sumw2();
	  h1_fC_ped_subt[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; TS; Average charge, ped subt [fC]",runnum.Data(),ieta+2,iphi));

	  h1_timing[ieta][iphi]= new TH1F(Form("h1_timing_ieta%i_iphi%i",ieta,iphi),Form("h1_timing_ieta%i_iphi%i",ieta,iphi), 50, -0.5, 9.5);
	  h1_timing[ieta][iphi]->Sumw2();
	  h1_timing[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Charge-averaged timing [TS]; Events",runnum.Data(),ieta+2,iphi));

	  h1_timing_laserconst[ieta][iphi]= new TH1F(Form("h1_timing_laserconst_ieta%i_iphi%i",ieta,iphi),Form("h1_timing_laserconst_ieta%i_iphi%i",ieta,iphi), 50, -0.5, 9.5);
	  h1_timing_laserconst[ieta][iphi]->Sumw2();
	  h1_timing_laserconst[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Charge-averaged timing [TS]; Events",runnum.Data(),ieta+2,iphi));
	    
	    
	  h1_tdc_timing[ieta][iphi]= new TH1F(Form("h1_tdc_timing_ieta%i_iphi%i",ieta,iphi),Form("h1_tdc_timing_ieta%i_iphi%i",ieta,iphi), 50, -0.5, 9.5);
	  h1_tdc_timing[ieta][iphi]->Sumw2();
	  h1_tdc_timing[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; TDC timing [TS]; Events",runnum.Data(),ieta+2,iphi));

	  h1_delta_timing[ieta][iphi]= new TH1F(Form("h1_delta_timing_ieta%i_iphi%i",ieta,iphi),Form("h1_delta_timing_ieta%i_iphi%i",ieta,iphi), 95, -9.5, 9.5);
	  h1_delta_timing[ieta][iphi]->Sumw2();
	  h1_delta_timing[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; TDC time - ADC time [TS]; Events",runnum.Data(),ieta+2,iphi));


	  h1_energy[ieta][iphi]= new TH1F(Form("h1_energy_ieta%i_iphi%i",ieta,iphi),Form("h1_energy_ieta%i_iphi%i",ieta,iphi),enbins,eminx,emaxx);
	  h1_energy[ieta][iphi]->Sumw2();
	  h1_energy[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Total charge [fC]; Events",runnum.Data(),ieta+2,iphi));
	  h1_energy[ieta][iphi]->StatOverflows(kTRUE);

	  h1_energy_laserconst[ieta][iphi]= new TH1F(Form("h1_energy_laserconst_ieta%i_iphi%i",ieta,iphi),Form("h1_energy_laserconst_ieta%i_iphi%i",ieta,iphi),enbins,eminx,emaxx);
	  h1_energy_laserconst[ieta][iphi]->Sumw2();
	  h1_energy_laserconst[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Total charge, constant window [fC]; Events",runnum.Data(),ieta+2,iphi));
	  h1_energy_laserconst[ieta][iphi]->StatOverflows(kTRUE);

	  h1_energy_const[ieta][iphi]= new TH1F(Form("h1_energy_const_ieta%i_iphi%i",ieta,iphi),Form("h1_energy_const_ieta%i_iphi%i",ieta,iphi),enbins,eminx,emaxx);
	  h1_energy_const[ieta][iphi]->Sumw2();
	  h1_energy_const[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Total charge, constant window [fC]; Events",runnum.Data(),ieta+2,iphi));
	  h1_energy_const[ieta][iphi]->StatOverflows(kTRUE);

	  h1_energy_const_zoom[ieta][iphi]= new TH1F(Form("h1_energy_const_zoom_ieta%i_iphi%i",ieta,iphi),Form("h1_energy_const_zoom_ieta%i_iphi%i",ieta,iphi),300,0,3000);
	  h1_energy_const_zoom[ieta][iphi]->Sumw2();
	  h1_energy_const_zoom[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Total charge, constant window [fC]; Events",runnum.Data(),ieta+2,iphi));
	  h1_energy_const_zoom[ieta][iphi]->StatOverflows(kTRUE);


	  h1_energy_const4[ieta][iphi]= new TH1F(Form("h1_energy_const4_ieta%i_iphi%i",ieta,iphi),Form("h1_energy_const4_ieta%i_iphi%i",ieta,iphi),enbins,eminx,emaxx);
	  h1_energy_const4[ieta][iphi]->Sumw2();
	  h1_energy_const4[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Total charge, constant window, 4 TS [fC]; Events",runnum.Data(),ieta+2,iphi));
	  h1_energy_const4[ieta][iphi]->StatOverflows(kTRUE);

	  h1_energy_no_ps[ieta][iphi]= new TH1F(Form("h1_energy_no_ps_ieta%i_iphi%i",ieta,iphi),Form("h1_energy_no_ps_ieta%i_iphi%i",ieta,iphi),enbins,eminx,emaxx);
	  h1_energy_no_ps[ieta][iphi]->Sumw2();
	  h1_energy_no_ps[ieta][iphi]->SetTitle(Form("Run %s, Fiber %i, Channel %i; Total charge, peak finding [fC]; Events",runnum.Data(),ieta+2,iphi));
	  h1_energy_no_ps[ieta][iphi]->StatOverflows(kTRUE);


        }
    }


  //Start event loop
  for(int ievent = 0; ievent<nevents; ievent++) 
    {
      tree->GetEntry(ievent); 
      if(calib && laserType==12) continue; //laserType == 12 means the laser fired into CRF during this event. This line is here take only CRF pedetal events in the abort gap.

      //For each event, loop over channels
      for(unsigned int ichan = 0; ichan<QIE11DigiIEta->size(); ichan++){
	int ieta =  QIE11DigiIEta->at(ichan);
	int iphi =  QIE11DigiIPhi->at(ichan);
	vector<float> QIE11DigiFC_ped_subt;

	//For each channel, loop over all time slices
        for(unsigned int its = 0; its<QIE11DigiFC->at(ichan).size(); its++){
	  float fC = QIE11DigiFC->at(ichan).at(its); //raw charge at this TS
	  int capid = QIE11DigiCapID->at(ichan).at(its); 
	  float fC_ped_subt = fC - peds[ieta][iphi][capid];
	  QIE11DigiFC_ped_subt.push_back(fC_ped_subt);
	 
	  h1_fC[ieta][iphi]->Fill(its,fC); //fill this timeslice to form the average pulse shapes
	  h1_ped[ieta][iphi]->Fill(its,peds[ieta][iphi][capid]);
	  h1_fC_ped_subt[ieta][iphi]->Fill(its,fC_ped_subt);

	  //cout<<Form("TS %i, fC %f",its,QIE11DigiFC->at(ichan).at(its))<<endl;
	}

	
	//Fill peak-finding charge distributions
	h1_energy[ieta][iphi]->Fill(energy(QIE11DigiFC_ped_subt));
	h1_energy_no_ps[ieta][iphi]->Fill(energy(QIE11DigiFC->at(ichan)));


	//Fill constant-window charge distributions (no pedestal subtraction). The second argument of energy() determines how many TS to use (starting from TS=2)
	h1_energy_const[ieta][iphi]->Fill(energy(QIE11DigiFC->at(ichan),8));
	h1_energy_const_zoom[ieta][iphi]->Fill(energy(QIE11DigiFC->at(ichan),8));
	h1_energy_const4[ieta][iphi]->Fill(energy(QIE11DigiFC->at(ichan),4));

	//find which event section this belongs in
	int section=0;
	for(int isec=0; isec<=nsections; isec++) {
	  if(isec*nevents/nsections > ievent) {section = isec-1; break;}
	}
	//	h1_energy_sections[ieta][iphi][section]->Fill(energy(QIE11DigiFC_ped_subt));
	h1_energy_sections[ieta][iphi][section]->Fill(energy(QIE11DigiFC->at(ichan),8)); //Currently using constant window, non-pedestal subtracted charge

	//fill timing info
	
	h1_timing[ieta][iphi]->Fill(time(QIE11DigiFC_ped_subt));
	//h1_tdc_timing[ieta][iphi]->Fill(tdc_time(QIE11DigiTDC->at(ichan)));
	//h1_delta_timing[ieta][iphi]->Fill(tdc_time(QIE11DigiTDC->at(ichan)) - time(QIE11DigiFC_ped_subt));

      }
    }
  //End event loop

  //Find peak timeslice for this run
  
  for(int ieta=0; ieta<nieta; ieta++) 
    { 
      for(int iphi=0; iphi<niphi; iphi++){ 
	peakTS[ieta][iphi] = h1_fC_ped_subt[ieta][iphi]->GetMaximumBin()-1; //Histogram bin numbers are offset by 1 from vector numbering...
	if(ieta==14&&iphi==5) cout<<"max is "<<peakTS[ieta][iphi]<<endl;
      }
    }

  //Redo event loop to fill the adjusted constant window

 for(int ievent = 0; ievent<nevents; ievent++) 
    {
      tree->GetEntry(ievent); 
      if(calib && laserType==12) continue; //laserType == 12 means the laser fired into CRF during this event. This line is here take only CRF pedetal events in the abort gap.

      //For each event, loop over channels
      for(unsigned int ichan = 0; ichan<QIE11DigiIEta->size(); ichan++){
	int ieta =  QIE11DigiIEta->at(ichan);
	int iphi =  QIE11DigiIPhi->at(ichan);
	vector<float> QIE11DigiFC_ped_subt;

	//For each channel, loop over all time slices
        for(unsigned int its = 0; its<QIE11DigiFC->at(ichan).size(); its++){
	  float fC = QIE11DigiFC->at(ichan).at(its); //raw charge at this TS
	  int capid = QIE11DigiCapID->at(ichan).at(its); 
	  float fC_ped_subt = fC - peds[ieta][iphi][capid];
	  QIE11DigiFC_ped_subt.push_back(fC_ped_subt);
	
	}

	//Fill constant-window charge distributions 
	h1_energy_laserconst[ieta][iphi]->Fill(energy(QIE11DigiFC_ped_subt,peakTS[ieta][iphi]-1,4));

	//fill timing info
	
	h1_timing_laserconst[ieta][iphi]->Fill(time(QIE11DigiFC_ped_subt,peakTS[ieta][iphi]-1,4));

      }
    }
  //End event loop


  
 


  if(calib) nevents = tree->GetEntries("laserType!=12");

  //Save histograms to root file, and normalize pulse shapes
  for(int ieta=0; ieta<nieta; ieta++) 
    { 
      if(ieta>=8&&ieta<=11) continue;
      for(int iphi=0; iphi<niphi; iphi++) 
        {	  
	  
	  h1_timing[ieta][iphi]->Write();
	  h1_timing_laserconst[ieta][iphi]->Write();
	  h1_tdc_timing[ieta][iphi]->Write();
	  h1_delta_timing[ieta][iphi]->Write();
	  h1_energy[ieta][iphi]->Write();
	  h1_energy_laserconst[ieta][iphi]->Write();
	  h1_energy_no_ps[ieta][iphi]->Write();
	  h1_energy_const[ieta][iphi]->Write();
	  h1_energy_const_zoom[ieta][iphi]->Write();
	  h1_energy_const4[ieta][iphi]->Write();

	  for(int isec=0;isec<nsections;isec++){
	    h1_energy_sections[ieta][iphi][isec]->Write();
	  }

	  //Global runs can have a different number of events for different channels, so find the correct normalization for each channel
	  if(collisions || calib) nevents = h1_energy_const[ieta][iphi]->GetEntries();

	  //normalize total pulse shapes to make average pulse shapes
	  h1_fC[ieta][iphi]->Scale(1./nevents);
	  h1_ped[ieta][iphi]->Scale(1./nevents);
	  h1_fC_ped_subt[ieta][iphi]->Scale(1./nevents);
	  h1_fC[ieta][iphi]->Write();	
	  h1_fC_ped_subt[ieta][iphi]->Write();
	

	
        }
    }

  if(make_plots){


    //This is a clunky way to pick which channels to show on a single frame
    // plot_distribution and plot_pulse take an array of TH1F

    //Remember that iEta = fiber num - 2, so iEta 14 is channel 16

    TH1F* petpen[6];
    petpen[0]=h1_fC[19][3];
    petpen[1]=h1_fC[19][5];
    petpen[2]=h1_fC[12][0];
    petpen[3]=h1_fC[19][4];
    petpen[4]=h1_fC[5][5];
    petpen[5]=h1_fC[2][3];

    TH1F* ped_petpen[6];
    ped_petpen[0]=h1_ped[19][3];
    ped_petpen[1]=h1_ped[19][5];
    ped_petpen[2]=h1_ped[12][0];
    ped_petpen[3]=h1_ped[19][4];
    ped_petpen[4]=h1_ped[5][5];
    ped_petpen[5]=h1_ped[2][3];


    TH1F* petpen2[6];
    petpen2[0]=h1_fC[12][1];
    petpen2[1]=h1_fC[19][1];
    petpen2[2]=h1_fC[19][2];
    petpen2[3]=h1_fC[18][5];
    petpen2[4]=h1_fC[3][0];
    petpen2[5]=h1_fC[1][0];

    TH1F* ped_petpen2[6];
    ped_petpen2[0]=h1_ped[12][1];
    ped_petpen2[1]=h1_ped[19][1];
    ped_petpen2[2]=h1_ped[19][2];
    ped_petpen2[3]=h1_ped[18][5];
    ped_petpen2[4]=h1_ped[3][0];
    ped_petpen2[5]=h1_ped[1][0];


    TH1F* other[6];
    other[0]=h1_fC[17][1]; //scintillator X
    other[1]=h1_fC[15][5]; //LS
    other[2]=h1_fC[14][2]; //LS
    other[3]=h1_fC[16][5]; //EJ200S
    other[4]=h1_fC[13][5]; //EJ200S
    other[5]=h1_fC[7][4];  //EJ200S

    TH1F* ped_other[6];
    ped_other[0]=h1_ped[17][1]; //scintillator X
    ped_other[1]=h1_ped[15][5]; //LS
    ped_other[2]=h1_ped[14][2]; //LS
    ped_other[3]=h1_ped[16][5]; //EJ200S
    ped_other[4]=h1_ped[13][5]; //EJ200S
    ped_other[5]=h1_ped[7][4];  //EJ200S




    //SCSN81-S channels in RM2, close to beam
    TH1F* rm2_close[6];
    rm2_close[0]=h1_fC[14][5];
    rm2_close[1]=h1_fC[13][1];
    rm2_close[2]=h1_fC[18][4];
    rm2_close[3]=h1_fC[18][5];
    rm2_close[4]=h1_fC[18][1];
    rm2_close[5]=h1_fC[0][4];
    
    TH1F* pedrm2_close[6];
    pedrm2_close[0]=h1_ped[14][5];
    pedrm2_close[1]=h1_ped[13][1];
    pedrm2_close[2]=h1_ped[18][4];
    pedrm2_close[3]=h1_ped[18][5];
    pedrm2_close[4]=h1_ped[18][1];
    pedrm2_close[5]=h1_ped[0][4];

    //SCSN81-S channels in RM2, far from beam
    TH1F * rm2_far[6]; 
    rm2_far[0]=h1_fC[14][5];
    rm2_far[1]=h1_fC[12][3];
    rm2_far[2]=h1_fC[14][1];
    rm2_far[3]=h1_fC[12][4];
    rm2_far[4]=h1_fC[13][4];
    rm2_far[5]=h1_fC[0][4];
    
    TH1F * pedrm2_far[6]; 
    pedrm2_far[0]=h1_ped[14][5];
    pedrm2_far[1]=h1_ped[12][3];
    pedrm2_far[2]=h1_ped[14][1];
    pedrm2_far[3]=h1_ped[12][4];
    pedrm2_far[4]=h1_ped[13][4];
    pedrm2_far[5]=h1_ped[0][4];


    //SCSN81-S channels in RM1, far from beam
    TH1F * rm1_far[6]; 
    rm1_far[0]=h1_fC[0][4];
    rm1_far[1]=h1_fC[5][2];
    rm1_far[2]=h1_fC[6][5];
    rm1_far[3]=h1_fC[1][2];
    rm1_far[4]=h1_fC[13][4];
    rm1_far[5]=h1_fC[14][5];

    TH1F * pedrm1_far[6]; 
    pedrm1_far[0]=h1_ped[0][4];
    pedrm1_far[1]=h1_ped[5][2];
    pedrm1_far[2]=h1_ped[6][5];
    pedrm1_far[3]=h1_ped[1][2];
    pedrm1_far[4]=h1_ped[13][4];
    pedrm1_far[5]=h1_ped[14][5];


    //SCSN81-S channels in RM1, close to beam
    TH1F * rm1_close[6]; 
    rm1_close[0]=h1_fC[0][4];
    rm1_close[1]=h1_fC[0][1];
    rm1_close[2]=h1_fC[3][0];
    rm1_close[3]=h1_fC[5][5];
    rm1_close[4]=h1_fC[1][0];
    rm1_close[5]=h1_fC[14][5];

    TH1F * pedrm1_close[6]; 
    pedrm1_close[0]=h1_ped[0][4];
    pedrm1_close[1]=h1_ped[0][1];
    pedrm1_close[2]=h1_ped[3][0];
    pedrm1_close[3]=h1_ped[5][5];
    pedrm1_close[4]=h1_ped[1][0];
    pedrm1_close[5]=h1_ped[14][5];

    

    //Empty channels, reference tiles, and 2 regular tiles
    TH1F * empty_ref_beam[6];
    empty_ref_beam[0]=h1_energy_const[2][5];
    empty_ref_beam[1]=h1_energy_const[0][4];
    empty_ref_beam[2]=h1_energy_const[3][0];
    empty_ref_beam[3]=h1_energy_const[4][5];
    empty_ref_beam[4]=h1_energy_const[14][5];
    empty_ref_beam[5]=h1_energy_const[18][4];
    

    TH1F * empty_ref_beam_zoom[6];
    empty_ref_beam_zoom[0]=h1_energy_const_zoom[2][5];
    empty_ref_beam_zoom[1]=h1_energy_const_zoom[0][4];
    empty_ref_beam_zoom[2]=h1_energy_const_zoom[3][0];
    empty_ref_beam_zoom[3]=h1_energy_const_zoom[4][5];
    empty_ref_beam_zoom[4]=h1_energy_const_zoom[14][5];
    empty_ref_beam_zoom[5]=h1_energy_const_zoom[18][4];
    

    TH1F * fC_empty_ref_beam[6];
    fC_empty_ref_beam[0]=h1_fC[2][5]; 
    fC_empty_ref_beam[1]=h1_fC[0][4]; 
    fC_empty_ref_beam[2]=h1_fC[3][0]; 
    fC_empty_ref_beam[3]=h1_fC[4][5]; 
    fC_empty_ref_beam[4]=h1_fC[14][5];
    fC_empty_ref_beam[5]=h1_fC[18][4];

    TH1F * ped_empty_ref_beam[6];
    ped_empty_ref_beam[0]=h1_ped[2][5];
    ped_empty_ref_beam[1]=h1_ped[0][4];
    ped_empty_ref_beam[2]=h1_ped[3][0];
    ped_empty_ref_beam[3]=h1_ped[4][5];
    ped_empty_ref_beam[4]=h1_ped[14][5];
    ped_empty_ref_beam[5]=h1_ped[18][4];
    
   
    
    //EJ260 channels
    TH1F * energy_ej260[6];
    energy_ej260[0]=h1_energy[14][5];
    energy_ej260[1]=h1_energy[7][2];
    energy_ej260[2]=h1_energy[1][1];
    energy_ej260[3]=h1_energy[6][2];
    energy_ej260[4]=h1_energy[3][5];

    TH1F * fC_ej260[6];
    fC_ej260[0]=h1_fC[14][5];
    fC_ej260[1]=h1_fC[7][2];
    fC_ej260[2]=h1_fC[1][1];
    fC_ej260[3]=h1_fC[6][2];
    fC_ej260[4]=h1_fC[3][5];

    TH1F * ped_ej260[6];
    ped_ej260[0]=h1_ped[14][5];
    ped_ej260[1]=h1_ped[7][2];
    ped_ej260[2]=h1_ped[1][1];
    ped_ej260[3]=h1_ped[6][2];
    ped_ej260[4]=h1_ped[3][5];




    //Plot all of these groups
    plot_pulse(petpen,ped_petpen,"petpen");
    plot_pulse(petpen2,ped_petpen2,"petpen2");
    plot_pulse(other,ped_other,"other");
    plot_pulse(fC_empty_ref_beam,ped_empty_ref_beam,"empty_ref_beam");
    plot_pulse(fC_ej260,ped_ej260,"EJ260S");
    plot_pulse(rm2_close,pedrm2_close,"rm2_close");
    plot_pulse(rm2_far,pedrm2_far,"rm2_far");
    plot_pulse(rm1_close,pedrm1_close,"rm1_close");
    plot_pulse(rm1_far,pedrm1_far,"rm1_far");

    plot_distribution(empty_ref_beam, "Energy_const","empty_ref_beam");
    plot_distribution(empty_ref_beam_zoom, "Energy_const","empty_ref_beam_zoom");
    plot_distribution(energy_ej260,"Energy","EJ260");


    //Besides the specific groups, also plot an entire fiber together (6 frames per page)
    for(int ieta=0; ieta<nieta; ieta++) 
      { 
	if(ieta>=8&&ieta<=11) continue;
	//	plot_pulse(h1_fC_ped_subt[ieta],Form("fiber%i",ieta+2));
	plot_pulse(h1_fC[ieta],h1_ped[ieta],Form("_no_ped_sub_fiber%i",ieta+2));
	plot_distribution(h1_energy[ieta],ieta, "Energy");
	plot_distribution(h1_energy_laserconst[ieta],ieta, "Energy_laser");
	//plot_distribution(h1_energy_no_ps[ieta],ieta, "Energy_no_ped_subt");
	plot_distribution(h1_energy_const[ieta],ieta, "Energy_const");
	//plot_distribution(h1_energy_const4[ieta],ieta, "Energy_const4");
       	plot_distribution(h1_timing[ieta],ieta, "Timing");
 	plot_distribution(h1_timing_laserconst[ieta],ieta, "Timing_laser");
      }
   


    //Even clunkier to specify channels for hists divided into run sections, since the arrays are 3D
    TH1F* rm1_close_sections[6][nsections];
    TH1F* rm1_far_sections[6][nsections];
    TH1F* rm2_close_sections[6][nsections];
    TH1F* rm2_far_sections[6][nsections];
    TH1F* empty_ref_beam_sections[6][nsections];

    for(int isec=0;isec<nsections;isec++){
      empty_ref_beam_sections[0][isec]=h1_energy_sections[2][5][isec];
      empty_ref_beam_sections[1][isec]=h1_energy_sections[0][4][isec];
      empty_ref_beam_sections[2][isec]=h1_energy_sections[3][0][isec];
      empty_ref_beam_sections[3][isec]=h1_energy_sections[4][5][isec];
      empty_ref_beam_sections[4][isec]=h1_energy_sections[14][5][isec];
      empty_ref_beam_sections[5][isec]=h1_energy_sections[18][4][isec];


      rm1_close_sections[0][isec]=h1_energy_sections[0][4][isec];
      rm1_close_sections[1][isec]=h1_energy_sections[0][1][isec];
      rm1_close_sections[2][isec]=h1_energy_sections[3][0][isec];
      rm1_close_sections[3][isec]=h1_energy_sections[5][5][isec];
      rm1_close_sections[4][isec]=h1_energy_sections[1][0][isec];
      rm1_close_sections[5][isec]=h1_energy_sections[14][5][isec];
     
      rm1_far_sections[0][isec]=h1_energy_sections[0][4][isec];
      rm1_far_sections[1][isec]=h1_energy_sections[5][2][isec];
      rm1_far_sections[2][isec]=h1_energy_sections[6][5][isec];
      rm1_far_sections[3][isec]=h1_energy_sections[1][2][isec];
      rm1_far_sections[4][isec]=h1_energy_sections[13][4][isec];
      rm1_far_sections[5][isec]=h1_energy_sections[14][5][isec];
    
      rm2_close_sections[0][isec]=h1_energy_sections[14][5][isec];
      rm2_close_sections[1][isec]=h1_energy_sections[13][1][isec];
      rm2_close_sections[2][isec]=h1_energy_sections[18][4][isec];
      rm2_close_sections[3][isec]=h1_energy_sections[18][5][isec];
      rm2_close_sections[4][isec]=h1_energy_sections[18][1][isec];
      rm2_close_sections[5][isec]=h1_energy_sections[0][4][isec];
  
      rm2_far_sections[0][isec]=h1_energy_sections[14][5][isec];
      rm2_far_sections[1][isec]=h1_energy_sections[12][3][isec];
      rm2_far_sections[2][isec]=h1_energy_sections[14][1][isec];
      rm2_far_sections[3][isec]=h1_energy_sections[12][4][isec];
      rm2_far_sections[4][isec]=h1_energy_sections[13][4][isec];
      rm2_far_sections[5][isec]=h1_energy_sections[14][5][isec];

    
    }
   



   
    plot_vs_section(h1_energy_sections[1],"fiber3");
    plot_vs_section(h1_energy_sections[18],"fiber20");
    plot_vs_section(empty_ref_beam_sections,"empty_ref_beam");

    plot_vs_section(rm2_close_sections,"rm2_close");
    plot_vs_section(rm2_far_sections,"rm2_far");
    plot_vs_section(rm1_close_sections,"rm1_close");
    plot_vs_section(rm1_far_sections,"rm1_far");



  }
  file->Close();
}


//Plot average pulse shape for 6 channels
void plot_pulse(TH1F* h1_vec[],TH1F* h1_ped[], TString type)
{
	TCanvas *c = new TCanvas("c", "c", 1200, 600); 

	c->Divide(3,2);
	int max = 6;
	if(type.Contains("260")) max=5;
	for(int iphi=0; iphi<max; iphi++) 
	  {
	    //  cout<<"pad number "<<iphi<<endl;
	    TPad *pad(NULL);
	    pad = static_cast<TPad *>(c->cd(iphi+1));
	    pad->SetLeftMargin(0.2);
	    pad->SetBottomMargin(0.2);
	    pad->SetGrid();
	    // TLatex *peak = new TLatex(0.7,0.84,Form("Peak TS: %i",));
	    //peak->SetNDC();
	    //peak->SetTextSize(textSize);
 

	    
	    h1_vec[iphi]->SetStats(0);
	    gStyle->SetTitleFontSize(0.1);
	    h1_vec[iphi]->GetXaxis()->SetTitleSize(0.07);
	    h1_vec[iphi]->GetXaxis()->SetLabelSize(0.06);
	    h1_vec[iphi]->GetYaxis()->SetTitleSize(0.07);
	    h1_vec[iphi]->GetYaxis()->SetLabelSize(0.06);
	    h1_vec[iphi]->GetYaxis()->SetTitleOffset(1.15);
	    h1_vec[iphi]->SetMinimum(0);
	    h1_vec[iphi]->Draw("EP"); 
	    h1_ped[iphi]->SetLineColor(kRed);
	    h1_ped[iphi]->SetLineStyle(3);
	    h1_ped[iphi]->Draw("hist same");

	  }
	
	if(!collisions) c->Print(Form("plots/Pulse_shape/Pulse_shape_run%s_%s.pdf",runnum.Data(),type.Data()));
	else  c->Print(Form("plots/Pulse_shape/Pulse_shape_run%s_%s_collisions.pdf",runnum.Data(),type.Data()));
	delete c;
}


void h1cosmetic(TH1F *hist){
  
  
  hist->GetXaxis()->SetTitleSize(0.07);
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.07);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetTitleOffset(1.15);
 

}



void plot_vs_section(TH1F * h1_vec[][nsections], TString type)
{
  TCanvas *c = new TCanvas("c", "c", 1200, 600); 
  c->Divide(3,2);
  
  for(int ifr=0; ifr<6; ifr++){
    TPad *pad(NULL);
    pad = static_cast<TPad *>(c->cd(ifr+1));
    pad->SetLeftMargin(0.2);
    pad->SetGrid();

    float x[nsections];
    float ex[nsections];
    float y[nsections];
    float ey[nsections];	  
    
    for(int isec=0;isec<nsections;isec++){
      y[isec] = h1_vec[ifr][isec]->GetMean();
      ey[isec] = h1_vec[ifr][isec]->GetMeanError();
      x[isec]=isec+1;
      ex[isec]=0.5;
    }

    TGraphErrors * g = new TGraphErrors(nsections,x,y,ex,ey);
   
    //Find channel from title string
    string rootsucks = h1_vec[ifr][0]->GetTitle();	   
    //Find out length of iEta by subtracting difference between phi and eta location
    int length = rootsucks.find("phi") - 4 - rootsucks.find("eta"); 
    TString fiber =  Form("%i",atoi(rootsucks.substr(rootsucks.find("eta")+3,length).c_str()) + 2 );	    
	    
    string phi = rootsucks.substr(rootsucks.find("phi")+3,1);
   
    g->SetMarkerSize(1);
    g->SetMarkerStyle(20);
    g->GetXaxis()->SetRangeUser(0.5,nsections+0.5);
    g->SetLineColor(kBlue+2);
    g->SetTitle(Form("Run %s, Fiber %s, Channel %s; Run section; Average charge [fC]",runnum.Data(),fiber.Data(),phi.c_str()));
    g->Draw("AEPZ");

  }
  if (!collisions) c->Print(Form("plots/gain_vs_section/gain_vs_section_run%s_%s.pdf",runnum.Data(),type.Data()));
  else  c->Print(Form("plots/gain_vs_section/gain_vs_section_run%s_%s_collisions.pdf",runnum.Data(),type.Data()));
  delete c;

}







void plot_distribution(TH1F* h1_ieta[], int ieta, TString type){
  TCanvas *c = new TCanvas("c", "c", 1200, 600); 
  c->Divide(3,2);
  
  for(int iphi=0; iphi<niphi; iphi++) 
    {
      TPad *pad(NULL);
      pad = static_cast<TPad *>(c->cd(iphi+1));
      pad->SetLeftMargin(0.2);
      pad->SetBottomMargin(0.2);
      pad->SetGrid();
      h1_ieta[iphi]->SetStats(0);
      
      h1cosmetic(h1_ieta[iphi]);
      gStyle->SetTitleFontSize(0.1);
      TString prec = "0";
      
      if(type.Contains("Timing")){
	//	xpos=0.25;
	prec="2";
      }

      TString n =Form("Entries: %."+prec+"f",h1_ieta[iphi]->GetEntries());
      TString m=Form("Mean: %."+prec+"f",h1_ieta[iphi]->GetMean());
      TString r= Form("RMS: %."+prec+"f",h1_ieta[iphi]->GetRMS());
      TString o= Form("Overflow: %.0f",h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1));
      TString u= Form("Underflow: %.0f",h1_ieta[iphi]->GetBinContent(0));

   
      if(type.Contains("Energy")){
	float xmax =  h1_ieta[iphi]->GetBinCenter(h1_ieta[iphi]->FindLastBinAbove(0)+1);
	if(!collisions){
	  if(xmax>900) h1_ieta[iphi]->Rebin();
	  if(xmax>2000) h1_ieta[iphi]->Rebin();
	  if(xmax>3000) h1_ieta[iphi]->Rebin();}
	h1_ieta[iphi]->GetXaxis()->SetRangeUser(0,xmax);
      }
      h1_ieta[iphi]->GetXaxis()->SetNdivisions(505,"X");
     

      if(h1_ieta[iphi]->GetBinContent(0) > 0) h1_ieta[iphi]->SetBinContent(1,h1_ieta[iphi]->GetBinContent(0)+h1_ieta[iphi]->GetBinContent(1));
      if(h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1) > 0) h1_ieta[iphi]->SetBinContent(h1_ieta[iphi]->GetNbinsX(),h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1)+h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()));


      h1_ieta[iphi]->Draw("hist"); 
	    
      float textSize = 0.055;
      float xpos=0.63;

      

      TLatex *ent = new TLatex(xpos,0.84,n);
      TLatex *mean = new TLatex(xpos,0.78,m);
      TLatex *rms = new TLatex(xpos,0.72,r);
      TLatex *of = new TLatex(xpos,0.66,o);
      TLatex *uf = new TLatex(xpos,0.6,u);



      ent->SetNDC();
      mean->SetNDC();
      of->SetNDC();
      rms->SetNDC();
      uf->SetNDC();
      ent->SetTextSize(textSize);
      mean->SetTextSize(textSize);
      of->SetTextSize(textSize); 
      rms->SetTextSize(textSize); 
      of->SetTextColor(kRed);
      uf->SetTextSize(textSize); 
      uf->SetTextColor(kRed);
      mean->Draw();
      rms->Draw();
      ent->Draw();
      if(h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1) > 0) of->Draw();
      if(h1_ieta[iphi]->GetBinContent(0) > 0) uf->Draw();
    }


  if(!collisions) c->Print(Form("plots/%s/%s_run%s_fiber%i.pdf",type.Data(),type.Data(),runnum.Data(),ieta+2));
  else c->Print(Form("plots/%s/%s_run%s_fiber%i_collisions.pdf",type.Data(),type.Data(),runnum.Data(),ieta+2));

  for(int iphi=0; iphi<niphi; iphi++) 
    {
      h1_ieta[iphi]->GetYaxis()->SetNdivisions(10); 
      TPad *pad(NULL);
      pad = static_cast<TPad *>(c->cd(iphi+1));
      pad->SetLogy();
      // if (collisions) pad->SetLogx();
    }
  if(!collisions) c->Print(Form("plots/%s/log_%s_run%s_fiber%i.pdf",type.Data(),type.Data(),runnum.Data(),ieta+2));
  else c->Print(Form("plots/%s/log_%s_run%s_fiber%i_collisions.pdf",type.Data(),type.Data(),runnum.Data(),ieta+2));
  delete c;
}

void plot_distribution(TH1F* h1_ieta[], TString type, TString name){
  TCanvas *c = new TCanvas("c", "c", 1200, 600); 
  
  int max = 4;
  if(name.Contains("260")){ max=5; c->Divide(3,2);}
  else if(name.Contains("empty_ref_beam")){ max= 6; c->Divide(3,2);}
  else c->Divide(2,2);
  for(int iphi=0; iphi<max; iphi++) 
    {
      TPad *pad(NULL);
      pad = static_cast<TPad *>(c->cd(iphi+1));
      pad->SetLeftMargin(0.2);
      pad->SetBottomMargin(0.2);
      pad->SetGrid();
      
      // gStyle->SetOptStat("Mo");
      //c->Update();
      h1_ieta[iphi]->SetStats(0);
      
      h1cosmetic(h1_ieta[iphi]);
      gStyle->SetTitleFontSize(0.1);

      TString prec = "0";
      
      
      if(type.Contains("Timing")){
	//	xpos=0.25;
	prec="2";
      }

      TString n =Form("Entries: %."+prec+"f",h1_ieta[iphi]->GetEntries());
      TString m=Form("Mean: %."+prec+"f",h1_ieta[iphi]->GetMean());
      TString r= Form("RMS: %."+prec+"f",h1_ieta[iphi]->GetRMS());
      TString o= Form("Overflow: %.0f",h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1));
      TString u= Form("Underflow: %.0f",h1_ieta[iphi]->GetBinContent(0));




      if(type.Contains("Energy")){
	float xmax =  h1_ieta[iphi]->GetBinCenter(h1_ieta[iphi]->FindLastBinAbove(0)+1);
       	if(!name.Contains("empty_ref_beam")){
	  if(xmax>900) h1_ieta[iphi]->Rebin();
	  if(xmax>1000) h1_ieta[iphi]->Rebin();
	  if(xmax>3000) h1_ieta[iphi]->Rebin();}
	
	if(!type.Contains("zoom")) h1_ieta[iphi]->GetXaxis()->SetRangeUser(0,xmax);
	h1_ieta[iphi]->GetXaxis()->SetNdivisions(505,"X");
	
      }

      if(h1_ieta[iphi]->GetBinContent(0) > 0) h1_ieta[iphi]->SetBinContent(1,h1_ieta[iphi]->GetBinContent(0)+h1_ieta[iphi]->GetBinContent(1));
      if(h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1) > 0) h1_ieta[iphi]->SetBinContent(h1_ieta[iphi]->GetNbinsX(),h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1)+h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()));
      h1_ieta[iphi]->Draw("hist"); 
	    
      float textSize = 0.055;
      float xpos=0.63;
     

      TLatex *ent = new TLatex(xpos,0.84,n);
      TLatex *mean = new TLatex(xpos,0.78,m);
      TLatex *rms = new TLatex(xpos,0.72,r);
      TLatex *of = new TLatex(xpos,0.66,o);
      TLatex *uf = new TLatex(xpos,0.6,u);


      ent->SetNDC();
      mean->SetNDC();
      of->SetNDC();
      uf->SetNDC();
      rms->SetNDC();
      ent->SetTextSize(textSize);
      mean->SetTextSize(textSize);
      of->SetTextSize(textSize); 
      rms->SetTextSize(textSize); 
      of->SetTextColor(kRed);
      uf->SetTextSize(textSize); 
      uf->SetTextColor(kRed);
      mean->Draw();
      rms->Draw();
      ent->Draw();
      if(h1_ieta[iphi]->GetBinContent(h1_ieta[iphi]->GetNbinsX()+1) > 0) of->Draw();
      if(h1_ieta[iphi]->GetBinContent(0) > 0) uf->Draw();
    }

  
  if(!collisions) c->Print(Form("plots/%s/%s_%s_run%s.pdf",type.Data(),type.Data(),name.Data(),runnum.Data()));
  else c->Print(Form("plots/%s/%s_%s_run%s_collisions.pdf",type.Data(),type.Data(),name.Data(),runnum.Data()));
  for(int iphi=0; iphi<max; iphi++) 
    {
      h1_ieta[iphi]->GetYaxis()->SetNdivisions(10); 
      TPad *pad(NULL);
      pad = static_cast<TPad *>(c->cd(iphi+1));
      pad->SetLogy();
      // if (collisions) pad->SetLogx();
    }
  if(!collisions) c->Print(Form("plots/%s/log_%s_%s_run%s.pdf",type.Data(),type.Data(),name.Data(),runnum.Data()));
  else c->Print(Form("plots/%s/log_%s_%s_run%s_collisions.pdf",type.Data(),type.Data(),name.Data(),runnum.Data()));
 
  /*for(int iphi=0; iphi<max; iphi++) 
    {
      //h1_ieta[iphi]->GetYaxis()->SetNdivisions(10); 
      TPad *pad(NULL);
      pad = static_cast<TPad *>(c->cd(iphi+1));
      pad->SetLogy(0);
      h1_ieta[iphi]->GetXaxis()->SetRangeUser(0,3000);
    }
  if(!collisions) c->Print(Form("plots/%s/zoom_%s_%s_run%s.pdf",type.Data(),type.Data(),name.Data(),runnum.Data()));
  else c->Print(Form("plots/%s/zoom_%s_%s_run%s_collisions.pdf",type.Data(),type.Data(),name.Data(),runnum.Data()));
  */


  delete c;
}
