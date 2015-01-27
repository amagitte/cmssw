#include "FWCore/Framework/interface/EventSetup.h"
#include "CondFormats/RecoMuonObjects/interface/DYTThrObject.h"
#include "CondFormats/DataRecord/interface/DYTThrObjectRcd.h"
#include "CondFormats/RecoMuonObjects/interface/DYTParamObject.h"
#include "CondFormats/RecoMuonObjects/interface/DYTParamsObject.h"
#include "CondFormats/DataRecord/interface/DYTParamsObjectRcd.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include "DataFormats/TrackingRecHit/interface/AlignmentPositionError.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
#include "CondFormats/Alignment/interface/Alignments.h"
#include "CondFormats/Alignment/interface/AlignTransform.h"
#include "CondFormats/AlignmentRecord/interface/DTAlignmentRcd.h"
#include "CondFormats/AlignmentRecord/interface/CSCAlignmentRcd.h"
#include "CondFormats/Alignment/interface/AlignmentErrors.h"
#include "CondFormats/Alignment/interface/AlignTransformError.h"
#include "CondFormats/AlignmentRecord/interface/DTAlignmentErrorRcd.h"
#include "CondFormats/AlignmentRecord/interface/CSCAlignmentErrorRcd.h"

class ThrParameters {

 public:
  ThrParameters(const edm::EventSetup*);
  ~ThrParameters();

  //Thrs Method
  void setInitialThr(double thr0) {x0 = thr0;};
  const bool isValidThdDB() {return isValidThdDB_;};
  const std::map<DTChamberId, GlobalError>& GetDTApeMap() {return dtApeMap;};
  const std::map<CSCDetId, GlobalError>& GetCSCApeMap() {return cscApeMap;};
  const DYTThrObject* getInitialThresholds() {return dytThresholds;}

  //Parameters Method
  const bool isValidParamDB() {return isValidParamDB_;};
  const std::map<uint32_t, std::vector<double> >& GetParametersMap() {return ParametersMap;}

 private:
  double x0;
  bool isValidThdDB_;
  bool isValidParamDB_;
  const DYTThrObject* dytThresholds;
  const DYTParamsObject* dytParams;
  std::vector<DYTParamObject> dytParameters;
  std::string function;
  std::map<DTChamberId, GlobalError> dtApeMap;
  std::map<CSCDetId, GlobalError> cscApeMap;
  std::map<uint32_t,std::vector<double> > ParametersMap;
};
