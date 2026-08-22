#include <FillHistogram.h>
#include <GHistogramer.h>
#include <GBCS.h>
#include <ddasHit.h>
#include <TOFCorrector.h>
#include <globals.h>
#include <utils.h>

#include <TString.h>

#include <cstdio>
#include <memory>
#include <filesystem>
#include <algorithm>


namespace {
  std::unique_ptr<TOFCorrector> gTOF;
  bool gTOFLoaded = false;
}


void LoadTOF(int run) {

  std::string path = programPath();
  std::string correctionFile = Form("%s/../tof/tof%04d.tof", path.c_str(), run);

  if(std::filesystem::exists(correctionFile)) {
    printf("Loading TOF correction: %s\n", correctionFile.c_str());
    gTOF = std::make_unique<TOFCorrector>(correctionFile);
  } else {
    printf(YELLOW "No TOF correction for run %04d - PID will use RAW tof\n"
           RESET_COLOR, run);
  }

  gTOFLoaded = true;
}


double CorrectedTOF(const GBCS& bcs) {

  if(!gTOFLoaded) {
    printf(RED "CorrectedTOF() called before LoadTOF() - using RAW tof\n"
           RESET_COLOR);
    gTOFLoaded = true;   // complain once, not once per event
  }

  double rawtof = bcs.fI2SPin1.fCharge;
  if(!gTOF) return rawtof;

  return gTOF->Correct(rawtof, bcs.fPin1.fTimestamp / 1.e8);
}


// Fill Histograms

void FillHistograms(const GBCS& bcs, const std::vector<ddasHit>& hits){

  Channel(hits);
  PID(bcs);
  Gamma(bcs);
  Other(bcs, hits);
}


void Channel(const std::vector<ddasHit>& hits){

  for(const auto& hit : hits) {
    GHistogramer::Get().Fill("All_Channel/ecal", 10000,0,32000, hit.GetEcal(),
                                                 300,0,300,     hit.GetId());

    GHistogramer::Get().Fill("All_Channel/raw",  10000,0,32000, hit.GetCharge(),
                                                 300,0,300,     hit.GetId());
  }
}


void PID(const GBCS& bcs){
  double tof = CorrectedTOF(bcs);
  double dE  = bcs.fPin1.fEcal;

  GHistogramer::Get().Fill("PID/PID_Total", 3600,0,24000, tof,
                                            1800,0,12000, dE);

}


void Gamma(const GBCS& bcs){

}


void Other(const GBCS& bcs, const std::vector<ddasHit>& hits){

}
