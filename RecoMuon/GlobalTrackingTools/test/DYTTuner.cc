#include <memory>
#include <algorithm>
#include <sstream>
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CondFormats/RecoMuonObjects/interface/DYTThrObject.h"
#include "CondFormats/DataRecord/interface/DYTThrObjectRcd.h" 
#include "CondCore/DBOutputService/interface/PoolDBOutputService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/DYTInfo.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#define MAX_THR 1e7

class DYTTuner : public edm::EDAnalyzer {
public:
  explicit DYTTuner(const edm::ParameterSet&);
  ~DYTTuner();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void   beginJob() ;
  virtual void   analyze(const edm::Event&, const edm::EventSetup&);
  virtual void   endJob() ;
  virtual double doIntegral(std::vector<double>&, DetId&);
  virtual double doIntegral(std::vector<double>&, int);
  virtual void   writePlots();
  virtual void   beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void   endRun(edm::Run const&, edm::EventSetup const&);
  virtual void   beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void   endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  typedef edm::ValueMap<reco::DYTInfo> DYTestimators;
  edm::Service<cond::service::PoolDBOutputService> poolDbService;
  double MinEnVal, MaxEnVal, MaxEstVal;
  std::vector<double> IntegralCut;
  unsigned int MinNumValues, MaxValPlots, NBinsPlots;
  bool saveROOTfile;
  DYTThrObject* thresholds;
  std::map<DetId, std::vector<double> > mapId;
  double etaMu;
  std::map<double, double> mapEta1Station;
  std::map<double, double> mapEta2Station;
  std::map<double, double> mapEta3Station;
  std::map<double, double> mapEta4Station;
  std::map<int, std::vector<double> > mapRegionEta1Station;
  std::map<int, std::vector<double> > mapRegionEta2Station;
  std::map<int, std::vector<double> > mapRegionEta3Station;
  std::map<int, std::vector<double> > mapRegionEta4Station;
  std::map<int, double> mapRegionEta1StationThr;
  std::map<int, double> mapRegionEta2StationThr;
  std::map<int, double> mapRegionEta3StationThr;
  std::map<int, double> mapRegionEta4StationThr;
  std::map<DetId, TH1F*> EstPlots;
  edm::EDGetTokenT<DYTestimators> dytInfoToken;
  edm::EDGetTokenT<reco::MuonCollection> muonsToken;
};


DYTTuner::DYTTuner(const edm::ParameterSet& iConfig)
{
  MinEnVal     = iConfig.getParameter<double>("MinEnergyVal");
  MaxEnVal     = iConfig.getParameter<double>("MaxEnergyVal");
  IntegralCut  = iConfig.getParameter< std::vector<double> >("IntegralCut");
  MaxEstVal    = iConfig.getParameter<double>("MaxEstVal");
  MinNumValues = iConfig.getParameter<unsigned int>("MinNumValues");
  saveROOTfile = iConfig.getParameter<bool>("writePlots");  
  MaxValPlots  = iConfig.getParameter<unsigned int>("MaxValPlots");
  NBinsPlots   = iConfig.getParameter<unsigned int>("NBinsPlots");

  edm::ConsumesCollector iC  = consumesCollector();
  dytInfoToken=iC.consumes<DYTestimators>(edm::InputTag("tevMuons", "dytInfo"));
  muonsToken=iC.consumes<reco::MuonCollection>(edm::InputTag("muons"));

  if (MaxEstVal == -1) MaxEstVal = MAX_THR;

  if (!poolDbService->isNewTagRequest("DYTThrObjectRcd")) 
    throw cms::Exception("NotAvailable") << "The output file already contains a valid \"DYTThrObjectRcd\" record.\nPlease provide a different file name or tag.";
}


DYTTuner::~DYTTuner()
{
}


void DYTTuner::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace reco;

  Handle<DYTestimators> dytInfoH;
  iEvent.getByToken(dytInfoToken, dytInfoH);
  const DYTestimators &dytInfoC = *dytInfoH;
  Handle<MuonCollection> muons;
  iEvent.getByToken(muonsToken, muons);
  for(size_t i = 0; i != muons->size(); ++i) {
    
    // Energy range cut
    if (muons->at(i).pt() < MinEnVal || muons->at(i).pt() > MaxEnVal) continue;

    etaMu = muons->at(i).eta();
    const TrackRef& tkRef = muons->at(i).globalTrack();
    if(dytInfoC.contains(tkRef.id())){
      DYTInfo dytInfo = dytInfoC[muons->at(i).globalTrack()];
      vector<double> estimators = dytInfo.DYTEstimators();
      vector<DetId> ids         = dytInfo.IdChambers();  
      for (int j = 0; j < 4; j++) {
	if (ids[j].null()) continue;
	DetId chamberId = ids[j];
	double estValue = estimators[j];
	if (estValue >= 0 && estValue <= MaxEstVal)
	  mapId[chamberId].push_back(estValue);
          if( chamberId.subdetId() == MuonSubdetId::DT ){
             int station = DTChamberId(chamberId).station();
             if ( station == 1 ){ mapEta1Station[etaMu] = estValue; }
             else if ( station == 2 ){ mapEta2Station[etaMu] = estValue; }
             else if ( station == 3 ){ mapEta3Station[etaMu] = estValue; }
             else if ( station == 4 ){ mapEta4Station[etaMu] = estValue; }
          } else if( chamberId.subdetId() ==  MuonSubdetId::CSC ){
             int station = CSCDetId(chamberId).station();    
             if ( station == 1 ){ mapEta1Station[etaMu] = estValue; }
             else if ( station == 2 ){ mapEta2Station[etaMu] = estValue; }
             else if ( station == 3 ){ mapEta3Station[etaMu] = estValue; }
             else if ( station == 4 ){ mapEta4Station[etaMu] = estValue; }
          }

      }
    } else {continue;}
  }
}


void DYTTuner::beginJob()
{
}


void DYTTuner::endJob() 
{
  if (saveROOTfile) writePlots();
  thresholds = new DYTThrObject();

  // Full barrel/endcap computation
  DetId id;
  std::map<DetId, std::vector<double> >::iterator it;
  std::map<int, std::vector<double> > estBarrel, estEndcap;
  for ( it = mapId.begin() ; it != mapId.end(); it++ ) {
    id = (*it).first;
    std::vector<double> estValCh = (*it).second;
    if ((*it).first.subdetId() == MuonSubdetId::DT) 
      for (unsigned int b = 0; b < estValCh.size(); b++) {
	int station = DTChamberId(id).station();
	estBarrel[station].push_back(estValCh[b]);
      }
    if ((*it).first.subdetId() == MuonSubdetId::CSC) 
      for (unsigned int e = 0; e < estValCh.size(); e++) {
	int station = CSCDetId(id).station();
	estEndcap[station].push_back(estValCh[e]);
      }
  }
 
  double barrelCut[4], endcapCut[4];
  for (unsigned st = 1; st <= 4; st++) {
    barrelCut[st-1] = doIntegral(estBarrel[st], id);
    endcapCut[st-1] = doIntegral(estEndcap[st], id);
  } 

  //Full Eta Region Computation from mapEta1Station -> mapRegionEta1Station 
  std::map<double, double>::iterator st1, st2, st3, st4;
  for ( st1 = mapEta1Station.begin(); st1 != mapEta1Station.end(); st1++ ){
      double eta = st1->first;
      double est = st1->second;
      if ( fabs(eta) >= 0 && fabs(eta) < 0.4 ) { mapRegionEta1Station[1].push_back(est); }
      if ( fabs(eta) >= 0.4 && fabs(eta) < 0.8 ) { mapRegionEta1Station[2].push_back(est); }
      if ( fabs(eta) >= 0.8 && fabs(eta) < 1.0 ) { mapRegionEta1Station[3].push_back(est); }
      if ( fabs(eta) >= 1.0 && fabs(eta) < 1.2 ) { mapRegionEta1Station[4].push_back(est); }
      if ( fabs(eta) >= 1.2 && fabs(eta) < 1.4 ) { mapRegionEta1Station[5].push_back(est); }
      if ( fabs(eta) >= 1.4 && fabs(eta) < 1.6 ) { mapRegionEta1Station[6].push_back(est); }
      if ( fabs(eta) >= 1.6 && fabs(eta) < 1.8 ) { mapRegionEta1Station[7].push_back(est); }
      if ( fabs(eta) >= 1.8 && fabs(eta) < 2.0 ) { mapRegionEta1Station[8].push_back(est); }
      if ( fabs(eta) >= 2.0 && fabs(eta) < 2.2 ) { mapRegionEta1Station[9].push_back(est); }
      if ( fabs(eta) >= 2.2 && fabs(eta) < 2.4 ) { mapRegionEta1Station[10].push_back(est); }   
  }

  for ( st2 = mapEta2Station.begin(); st2 != mapEta2Station.end(); st2++ ){
      double eta = st2->first;
      double est = st2->second;
      if ( fabs(eta) >= 0 && fabs(eta) < 0.4 ) { mapRegionEta2Station[1].push_back(est); }
      if ( fabs(eta) >= 0.4 && fabs(eta) < 0.8 ) { mapRegionEta2Station[2].push_back(est); }
      if ( fabs(eta) >= 0.8 && fabs(eta) < 1.0 ) { mapRegionEta2Station[3].push_back(est); }
      if ( fabs(eta) >= 1.0 && fabs(eta) < 1.2 ) { mapRegionEta2Station[4].push_back(est); }
      if ( fabs(eta) >= 1.2 && fabs(eta) < 1.4 ) { mapRegionEta2Station[5].push_back(est); }
      if ( fabs(eta) >= 1.4 && fabs(eta) < 1.6 ) { mapRegionEta2Station[6].push_back(est); }
      if ( fabs(eta) >= 1.6 && fabs(eta) < 1.8 ) { mapRegionEta2Station[7].push_back(est); }
      if ( fabs(eta) >= 1.8 && fabs(eta) < 2.0 ) { mapRegionEta2Station[8].push_back(est); }
      if ( fabs(eta) >= 2.0 && fabs(eta) < 2.2 ) { mapRegionEta2Station[9].push_back(est); }
      if ( fabs(eta) >= 2.2 && fabs(eta) < 2.4 ) { mapRegionEta2Station[10].push_back(est); }
  }

  for ( st3 = mapEta3Station.begin(); st3 != mapEta3Station.end(); st3++ ){
      double eta = st3->first;
      double est = st3->second;
      if ( fabs(eta) >= 0 && fabs(eta) < 0.4 ) { mapRegionEta3Station[1].push_back(est); }
      if ( fabs(eta) >= 0.4 && fabs(eta) < 0.8 ) { mapRegionEta3Station[2].push_back(est); }
      if ( fabs(eta) >= 0.8 && fabs(eta) < 1.0 ) { mapRegionEta3Station[3].push_back(est); }
      if ( fabs(eta) >= 1.0 && fabs(eta) < 1.2 ) { mapRegionEta3Station[4].push_back(est); }
      if ( fabs(eta) >= 1.2 && fabs(eta) < 1.4 ) { mapRegionEta3Station[5].push_back(est); }
      if ( fabs(eta) >= 1.4 && fabs(eta) < 1.6 ) { mapRegionEta3Station[6].push_back(est); }
      if ( fabs(eta) >= 1.6 && fabs(eta) < 1.8 ) { mapRegionEta3Station[7].push_back(est); }
      if ( fabs(eta) >= 1.8 && fabs(eta) < 2.0 ) { mapRegionEta3Station[8].push_back(est); }
      if ( fabs(eta) >= 2.0 && fabs(eta) < 2.2 ) { mapRegionEta3Station[9].push_back(est); }
      if ( fabs(eta) >= 2.2 && fabs(eta) < 2.4 ) { mapRegionEta3Station[10].push_back(est); }
  }

  for ( st4 = mapEta4Station.begin(); st4 != mapEta4Station.end(); st4++ ){
      double eta = st4->first;
      double est = st4->second;
      if ( fabs(eta) >= 0 && fabs(eta) < 0.4 ) { mapRegionEta4Station[1].push_back(est); }
      if ( fabs(eta) >= 0.4 && fabs(eta) < 0.8 ) { mapRegionEta4Station[2].push_back(est); }
      if ( fabs(eta) >= 0.8 && fabs(eta) < 1.0 ) { mapRegionEta4Station[3].push_back(est); }
      if ( fabs(eta) >= 1.0 && fabs(eta) < 1.2 ) { mapRegionEta4Station[4].push_back(est); }
      if ( fabs(eta) >= 1.2 && fabs(eta) < 1.4 ) { mapRegionEta4Station[5].push_back(est); }
      if ( fabs(eta) >= 1.4 && fabs(eta) < 1.6 ) { mapRegionEta4Station[6].push_back(est); }
      if ( fabs(eta) >= 1.6 && fabs(eta) < 1.8 ) { mapRegionEta4Station[7].push_back(est); }
      if ( fabs(eta) >= 1.8 && fabs(eta) < 2.0 ) { mapRegionEta4Station[8].push_back(est); }
      if ( fabs(eta) >= 2.0 && fabs(eta) < 2.2 ) { mapRegionEta4Station[9].push_back(est); }
      if ( fabs(eta) >= 2.2 && fabs(eta) < 2.4 ) { mapRegionEta4Station[10].push_back(est); }


  }

  // Chamber by chamber computation
  for ( it = mapId.begin() ; it != mapId.end(); it++ ) {
    DetId id = (*it).first;
    std::vector<double> estValCh = (*it).second;
    DYTThrObject::DytThrStruct obj;
    obj.id = id;
    if (estValCh.size() < MinNumValues) {
      if (id.subdetId() == MuonSubdetId::DT) {
	int station = DTChamberId(id).station();
	obj.thr = barrelCut[station-1];
      }
      if (id.subdetId() == MuonSubdetId::CSC) {
	int station = CSCDetId(id).station();
	obj.thr = endcapCut[station-1];
      }
      thresholds->thrsVec.push_back(obj);
      continue;
    }
    obj.thr = doIntegral(estValCh, id);
    thresholds->thrsVec.push_back(obj);
  }

  // Writing to DB
  edm::Service<cond::service::PoolDBOutputService> poolDbService;
  if( poolDbService.isAvailable() )
  {
    poolDbService->writeOne( thresholds, poolDbService->beginOfTime(), "DYTThrObjectRcd"  ); 
  }
  else throw cms::Exception("NotAvailable") << "PoolDBOutputService is not available.";

  //Integral in Eta region
  std::map<int, std::vector<double> >::iterator rg1, rg2, rg3, rg4;

  for ( rg1 = mapRegionEta1Station.begin(); rg1 != mapRegionEta1Station.end(); rg1++ ){
     int etarg1 = (*rg1).first;
     std::vector<double> vthr1 = (*rg1).second;
     if( vthr1.size() > MinNumValues ) mapRegionEta1StationThr[etarg1] = doIntegral(vthr1, etarg1);     
  }

  for ( rg2 = mapRegionEta2Station.begin(); rg2 != mapRegionEta2Station.end(); rg2++ ){
     int etarg2 = (*rg2).first;
     std::vector<double> vthr2 = (*rg2).second;
     if( vthr2.size() > MinNumValues ) mapRegionEta2StationThr[etarg2] = doIntegral(vthr2, etarg2);
  }

  for ( rg3 = mapRegionEta3Station.begin(); rg3 != mapRegionEta3Station.end(); rg3++ ){
     int etarg3 = (*rg3).first;
     std::vector<double> vthr3 = (*rg3).second;
     if( vthr3.size() > MinNumValues ) mapRegionEta3StationThr[etarg3] = doIntegral(vthr3, etarg3);
  }

  for ( rg4 = mapRegionEta4Station.begin(); rg4 != mapRegionEta4Station.end(); rg4++ ){
     int etarg4 = (*rg4).first;
     std::vector<double> vthr4 = (*rg4).second;
     if( vthr4.size() > MinNumValues ) mapRegionEta4StationThr[etarg4] = doIntegral(vthr4, etarg4); 
  }

  std::map<int, double>::iterator thr1, thr2, thr3, thr4;
  std::cout << "The format is the following: Station, Eta Region, Thr" <<std::endl;
  for ( thr1 = mapRegionEta1StationThr.begin(); thr1 != mapRegionEta1StationThr.end(); thr1++ ){
     std::cout << "1 " << (*thr1).first << " " << (*thr1).second << std::endl;
  }

  for ( thr2 = mapRegionEta2StationThr.begin(); thr2 != mapRegionEta2StationThr.end(); thr2++ ){
     std::cout << "2 " << (*thr2).first << " " << (*thr2).second << std::endl;
  }

  for ( thr3 = mapRegionEta3StationThr.begin(); thr3 != mapRegionEta3StationThr.end(); thr3++ ){
     std::cout << "3 " << (*thr3).first << " " << (*thr3).second << std::endl;
  }

  for ( thr4 = mapRegionEta4StationThr.begin(); thr4 != mapRegionEta4StationThr.end(); thr4++ ){
     std::cout << "4 " << (*thr4).first << " " << (*thr4).second << std::endl;
  }


}


double DYTTuner::doIntegral(std::vector<double>& estValues, DetId& id) {
  double cutOnE = -1;
  int nPosVal = 0;
  double IntegralCutUsed = 0;

  if( id.subdetId() == MuonSubdetId::DT ){
    IntegralCutUsed = IntegralCut[0];
  } else if( id.subdetId() == MuonSubdetId::CSC ){
    IntegralCutUsed = IntegralCut[1];
  }
  
  sort( estValues.begin(), estValues.end() );
  for (unsigned int j = 0; j < estValues.size(); j++) 
    if (estValues[j] > 0 && estValues[j] < MaxEstVal) nPosVal++;
  double limit = nPosVal * IntegralCutUsed;
  int nVal = 0; 
  for (unsigned int j = 0; j < estValues.size(); j++) {
    if (estValues[j] < 0) continue;
    nVal++;
    if (nVal >= limit) {
      cutOnE = estValues[j-1];
      break;
    }
  }
  //std::cout << "Det Id: " << id.rawId() << " - Threshold:: " << cutOnE << std::endl;
  return cutOnE;
}

double DYTTuner::doIntegral(std::vector<double>& estValues, int etaregion){
  double cutOnE = -1;
  int nPosVal = 0;
  double IntegralCutUsed = 0;
  
  if( etaregion > 0 && etaregion < 3 ){
     IntegralCutUsed = IntegralCut[0];
  } else if( etaregion > 2 && etaregion < 11 ){
     IntegralCutUsed = IntegralCut[1];
  }
 
  sort( estValues.begin(), estValues.end() );
  for (unsigned int j = 0; j < estValues.size(); j++)
    if (estValues[j] > 0 && estValues[j] < MaxEstVal) nPosVal++;
  double limit = nPosVal * IntegralCutUsed;
  int nVal = 0;
  for (unsigned int j = 0; j < estValues.size(); j++) {
    if (estValues[j] < 0) continue;
    nVal++;
    if (nVal >= limit) {
      cutOnE = estValues[j-1];
      break;
    }
  }
  //std::cout << "Det Id: " << id.rawId() << " - Threshold:: " << cutOnE << std::endl;
  return cutOnE;
}

void DYTTuner::writePlots() {
  edm::Service<TFileService> fs;
  std::map<DetId, std::vector<double> >::iterator it;
  for ( it = mapId.begin() ; it != mapId.end(); it++ ) {
    DetId id = (*it).first;
    std::vector<double> estValCh = (*it).second;
    int sector, station, wheel, ring;
    std::ostringstream sector_osstr, station_osstr, wheel_osstr, ring_osstr;
    std::string sector_str, station_str, wheel_str, ring_str;
    std::string plotName;
    if( id.subdetId() == MuonSubdetId::DT ){
      station = DTChamberId(id).station();
      wheel   = DTChamberId(id).wheel();
      sector  = DTChamberId(id).sector();
      station_osstr << station;
      wheel_osstr << wheel;
      sector_osstr << sector;
      station_str = station_osstr.str();
      wheel_str   = wheel_osstr.str();
      sector_str  = sector_osstr.str();
      plotName = "DT_Wheel_" + wheel_str + "_Station_" + station_str + "_Sector_" + sector_str;
 
    } else if( id.subdetId() == MuonSubdetId::CSC ){
      station = CSCDetId(id).station();
      sector   = CSCDetId(id).chamber();
      ring  = CSCDetId(id).ring();
      station_osstr << station;
      sector_osstr << sector;
      ring_osstr << ring;
      station_str = station_osstr.str();
      sector_str   = sector_osstr.str();
      ring_str  = ring_osstr.str();
      plotName = "CSC_Ring_" + ring_str + "_Station_" + station_str + "_Sector_" + sector_str;

    }

    TH1F* tmpPlot = new TH1F(plotName.c_str(), plotName.c_str(), NBinsPlots, 0., MaxValPlots);
    for (unsigned int i = 0; i < estValCh.size(); i++) 
      tmpPlot->Fill(estValCh[i]);
    EstPlots[id] = fs->make<TH1F>(*tmpPlot);
    delete tmpPlot;
  }
}


void DYTTuner::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

void DYTTuner::endRun(edm::Run const&, edm::EventSetup const&)
{
}


void DYTTuner::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


void DYTTuner::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}


void DYTTuner::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(DYTTuner);
