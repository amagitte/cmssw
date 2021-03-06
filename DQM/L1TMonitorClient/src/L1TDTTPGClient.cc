#include "DQM/L1TMonitorClient/interface/L1TDTTPGClient.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DQMServices/Core/interface/QReport.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"
#include "TRandom.h"

#include <TF1.h>
#include <stdio.h>
#include <sstream>
#include <math.h>
#include <TProfile.h>
#include <TProfile2D.h>

using namespace edm;
using namespace std;

L1TDTTPGClient::L1TDTTPGClient(const edm::ParameterSet& ps)
{
  parameters_=ps;
  initialize();
}

L1TDTTPGClient::~L1TDTTPGClient(){
 LogInfo("TriggerDQM")<<"[TriggerDQM]: ending... ";
}

//--------------------------------------------------------
void L1TDTTPGClient::initialize(){ 

  counterLS_=0; 
  counterEvt_=0; 
  
  // get back-end interface
  dbe_ = Service<DQMStore>().operator->();
  
  // base folder for the contents of this job
  monitorName_ = parameters_.getUntrackedParameter<string>("monitorName","");
//  cout << "Monitor name = " << monitorName_ << endl;
  prescaleLS_ = parameters_.getUntrackedParameter<int>("prescaleLS", -1);
//  cout << "DQM lumi section prescale = " << prescaleLS_ << " lumi section(s)"<< endl;
  prescaleEvt_ = parameters_.getUntrackedParameter<int>("prescaleEvt", -1);
//  cout << "DQM event prescale = " << prescaleEvt_ << " events(s)"<< endl;
  output_dir_ = parameters_.getUntrackedParameter<string>("output_dir","");
//  cout << "DQM output dir = " << output_dir_ << endl;
  input_dir_ = parameters_.getUntrackedParameter<string>("input_dir","");
//  cout << "DQM input dir = " << input_dir_ << endl;
  
  LogInfo( "TriggerDQM");

      
}

//--------------------------------------------------------
void L1TDTTPGClient::beginJob(void){

  LogInfo("TriggerDQM")<<"[TriggerDQM]: Begin Job";

  // get backendinterface
  dbe_ = Service<DQMStore>().operator->();  

  dbe_->setCurrentFolder(output_dir_);

  // booking
  
  dttpgphmapcorrf = dbe_->book2D("DT_TPG_phi_map_corr_frac",
				 "Fraction of correlated best triggers per station",20,1,21,12,0,12);
  dttpgphmap2ndf = dbe_->book2D("DT_TPG_phi_map_2nd_frac",
				"Fraction of second tracks per station",20,1,21,12,0,12);
  dttpgphmapbxf[0] = dbe_->book2D("DT_TPG_phi_map_bx-1_frac",
				  "Fraction of triggers per station (BX=-1)",20,1,21,12,0,12);
  dttpgphmapbxf[1] = dbe_->book2D("DT_TPG_phi_map_bx0_frac",
				  "Fraction of triggers per station (BX=0)",20,1,21,12,0,12);
  dttpgphmapbxf[2] = dbe_->book2D("DT_TPG_phi_map_bx+1_frac",
				  "Fraction of triggers per station (BX=1)",20,1,21,12,0,12);
  setMapPhLabel(dttpgphmapcorrf);
  setMapPhLabel(dttpgphmap2ndf);
  setMapPhLabel(dttpgphmapbxf[0]);
  setMapPhLabel(dttpgphmapbxf[1]);
  setMapPhLabel(dttpgphmapbxf[2]);

  dttpgthmaphf = dbe_->book2D("DT_TPG_theta_map_corr_frac",
				 "Fraction of H quality best triggers per station",15,1,16,12,0,12);
  dttpgthmapbxf[0] = dbe_->book2D("DT_TPG_theta_map_bx-1_frac",
				  "Fraction of triggers per station (BX=-1)",15,1,16,12,0,12);
  dttpgthmapbxf[1] = dbe_->book2D("DT_TPG_theta_map_bx0_frac",
				  "Fraction of triggers per station (BX=0)",15,1,16,12,0,12);
  dttpgthmapbxf[2] = dbe_->book2D("DT_TPG_theta_map_bx+1_frac",
				  "Fraction of triggers per station (BX=1)",15,1,16,12,0,12);
  setMapThLabel(dttpgthmaphf);
  setMapThLabel(dttpgthmapbxf[0]);
  setMapThLabel(dttpgthmapbxf[1]);
  setMapThLabel(dttpgthmapbxf[2]);

  

}

//--------------------------------------------------------
void L1TDTTPGClient::beginRun(const Run& r, const EventSetup& context) {
}

//--------------------------------------------------------
void L1TDTTPGClient::beginLuminosityBlock(const LuminosityBlock& lumiSeg, const EventSetup& context) {
   // optionally reset histograms here
   // clientHisto->Reset();
}
//--------------------------------------------------------

void L1TDTTPGClient::endLuminosityBlock(const edm::LuminosityBlock& lumiSeg, 
                          const edm::EventSetup& c){
			  
}			  
//--------------------------------------------------------
void L1TDTTPGClient::analyze(const Event& e, const EventSetup& context){
//   cout << "L1TDTTPGClient::analyze" << endl;
   counterEvt_++;
   if (prescaleEvt_<1) return;
   if (prescaleEvt_>0 && counterEvt_%prescaleEvt_ != 0) return;

   string nName = "DT_TPG_phi_best_map_corr";
   string dName = "DT_TPG_phi_best_map";
   makeRatioHisto(dttpgphmapcorrf,nName,dName);
   dName = "DT_TPG_phi_map";
   nName = "DT_TPG_phi_map_2nd";
   makeRatioHisto(dttpgphmap2ndf,nName,dName);
   nName = "DT_TPG_phi_map_bx-1";
   makeRatioHisto(dttpgphmapbxf[0],nName,dName);
   nName = "DT_TPG_phi_map_bx0";
   makeRatioHisto(dttpgphmapbxf[1],nName,dName);
   nName = "DT_TPG_phi_map_bx+1";
   makeRatioHisto(dttpgphmapbxf[2],nName,dName);

   nName = "DT_TPG_theta_best_map_h";
   dName = "DT_TPG_theta_best_map";
   makeRatioHisto(dttpgthmaphf,nName,dName);
   dName = "DT_TPG_theta_map";
   nName = "DT_TPG_theta_map_bx-1";
   makeRatioHisto(dttpgthmapbxf[0],nName,dName);
   nName = "DT_TPG_theta_map_bx0";
   makeRatioHisto(dttpgthmapbxf[1],nName,dName);
   nName = "DT_TPG_theta_map_bx+1";
   makeRatioHisto(dttpgthmapbxf[2],nName,dName);


}

//--------------------------------------------------------
void L1TDTTPGClient::endRun(const Run& r, const EventSetup& context){
}

//--------------------------------------------------------
void L1TDTTPGClient::endJob(){
}

void L1TDTTPGClient::makeRatioHisto(MonitorElement *ratioME, string &nName, string &dName)
{

   TH2F *numerator;
   TH2F *denominator;

   denominator = this->get2DHisto(input_dir_+"/"+dName,dbe_);
   numerator   = this->get2DHisto(input_dir_+"/"+nName,dbe_);

   if (numerator && denominator) {

     TH2F * ratio = ratioME->getTH2F();
     if (ratio) {
       ratio->Divide(numerator,denominator);
     }
     else {
       LogInfo("TriggerDQM") << "[TriggerDQM]: ratio histo named \"" << ratioME->getName() << "\" not found!" << endl;
     }
   }
   else {
     if (!numerator)
       LogInfo("TriggerDQM") << "[TriggerDQM]: numerator histo \"" << nName << "\" not found!" << endl;
     if (!denominator)
       LogInfo("TriggerDQM") << "[TriggerDQM]: denominator histo \"" << dName << "\" not found!" << endl;
   }

}

TH1F * L1TDTTPGClient::get1DHisto(string meName, DQMStore * dbi)
{

  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogInfo("TriggerDQM") << "ME NOT FOUND.";
    return NULL;
  }

  return me_->getTH1F();
}

TH2F * L1TDTTPGClient::get2DHisto(string meName, DQMStore * dbi)
{


  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogInfo("TriggerDQM") << "ME NOT FOUND.";
    return NULL;
  }

  return me_->getTH2F();
}



TProfile2D * L1TDTTPGClient::get2DProfile(string meName, DQMStore * dbi)
{


  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogInfo("TriggerDQM") << "ME NOT FOUND.";
    return NULL;
  }

  return me_->getTProfile2D();
}


TProfile * L1TDTTPGClient::get1DProfile(string meName, DQMStore * dbi)
{


  MonitorElement * me_ = dbi->get(meName);

  if (!me_) { 
    LogInfo("TriggerDQM") << "ME NOT FOUND.";
    return NULL;
  }

  return me_->getTProfile();
}


void L1TDTTPGClient::setMapPhLabel(MonitorElement *me)
{

  me->setAxisTitle("DTTF Sector",2);
      for(int i=0;i<5;i++){
	ostringstream wheel;
	wheel << i-2;
	me->setBinLabel(1+i*4,"Wheel "+ wheel.str(),1);
      }
  
}


void L1TDTTPGClient::setMapThLabel(MonitorElement *me)
{

  me->setAxisTitle("DTTF Sector",2);
      for(int i=0;i<5;i++){
	ostringstream wheel;
	wheel << i-2;
	me->setBinLabel(1+i*3,"Wheel "+ wheel.str(),1);
      }
  
}
