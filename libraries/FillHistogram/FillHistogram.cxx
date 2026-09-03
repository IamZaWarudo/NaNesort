#include <FillHistogram.h>
#include <GHistogramer.h>
#include <GBCS.h>
#include <ddasHit.h>
#include <TOFCorrector.h>
#include <globals.h>
#include <utils.h>

#include <TString.h>
#include <TKey.h>
#include <TCutG.h>

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
  double runtime = bcs.fPin1.fTimestamp / 1.e8 ;

  GHistogramer::Get().Fill("TOF/tof_corrected",500,0,5000, runtime,
                                               1500,0,25000, tof);
/*
if(bcs.EventType() == 1){
  GHistogramer::Get().Fill("TOF/tof_Implant",500,0,5000, runtime,
                                               1500,0,25000, tof);
}
else
if(bcs.EventType() == 3) {
  GHistogramer::Get().Fill("TOF/tof_Lightion",500,0,5000, runtime,
                                               1500,0,25000, tof);
}
else
if(bcs.EventType() == 4) {
  GHistogramer::Get().Fill("TOF/tof_Veto",500,0,5000, runtime,
                                               1500,0,25000, tof);
}
*/
  GHistogramer::Get().Fill("PID/PID_Total", 3600,0,24000, tof,
                                            1800,0,12000, dE);
 
  if(bcs.EventType() == 1) {   // IMPLANT
  GHistogramer::Get().Fill("PID/PID_Implant", 3600,0,24000, tof,
                                                1800,0,12000, dE);
  } 
}
 

void Gamma(const GBCS& bcs){

  if(bcs.EventType() == 1) {   // IMPLANT
  GHistogramer::Get().Fill("Clover/Gamma_Implant", 3600,0,24000, bcs.fDSSD.fEnergy,
                                                  3600,0,24000, bcs.fClover.fEcal);
  } else
  if(bcs.EventType() == 2) {  // DECAY
   GHistogramer::Get().Fill("Clover/Gamma_Decay", 3600,0,24000, bcs.fDSSD.fEnergy,
                                                  3600,0,24000, bcs.fClover.fEcal);
  } 
}


void LoadPIDs() { //TODO make me find the file better
  TDirectory *current = gDirectory;
  
  pids = new TList;

  TFile *cuts = TFile::Open("pid.cuts");
  TIter iter(cuts->GetListOfKeys());
  while(TKey *key =(TKey*)iter.Next()) {
    pids->Add(key->ReadObj());    
  }
  current->cd();
};



// ---------------------------------------------------------------
// Correlations
// ---------------------------------------------------------------
void Correlation(const GBCS& implant, const GBCS& decay, double dt){

  if(!pids) 
    LoadPIDs();

  TIter iter(pids);
  while(TCutG *cut = (TCutG*)iter.Next()) {
  
  const double tof = CorrectedTOF(implant);
  const double dE  = implant.fPin1.fEcal;
  if(!cut->IsInside(tof,dE)) continue;


  GHistogramer::Get().Fill(Form("%s/dt",cut->GetName()), 4000, -1000, 3000, dt);


  GHistogramer::Get().Fill(Form("%s/implant_dt_DSSD",cut->GetName()), 4000,-1000,3000, dt,
                           3600,0,64000, implant.fDSSD.fEnergy);

 for(const auto& g : decay.fClover.fCloverHits) {
  double energy = g.GetEcal();
  if(energy <= 0) continue;

  GHistogramer::Get().Fill(Form("%s/decay_gamma",cut->GetName()),1000,0,4000, energy);
  GHistogramer::Get().Fill(Form("%s/decay_dt_gamma",cut->GetName()),4000, -1000, 3000,dt,
                           1000,0,4000, energy);
 }

  for(const auto& g : implant.fClover.fCloverHits) {
  double energy = g.GetEcal();
  if(energy <= 0) continue;

  GHistogramer::Get().Fill(Form("%s/implant_gamma",cut->GetName()),1000,0,4000, energy);
  GHistogramer::Get().Fill(Form("%s/implant_tof_gamma",cut->GetName()), 1000,0,4000, energy,
                                                            2000,0,24000, CorrectedTOF(implant));

  }
  
  GHistogramer::Get().Fill(Form("%s/position",cut->GetName()), 40,0,40, decay.fDSSD.Xpos,
                           40,0,40, decay.fDSSD.Ypos);
  }
}


void CorrelationAllPairs(const GBCS& implant, const GBCS& decay, double dt){

  if(!pids) 
    LoadPIDs();

  TIter iter(pids);
  while(TCutG *cut = (TCutG*)iter.Next()) {

  GHistogramer::Get().Fill(Form("%s/dt_allPairs",cut->GetName()), 4000, -1000, 3000, dt);

   for(const auto& g : decay.fClover.fCloverHits) {
    double energy = g.GetEcal();
    if(energy <= 0) continue;
  
    GHistogramer::Get().Fill(Form("%s/decay_dt_gamma_allPairs",cut->GetName()),4000, -1000, 3000,dt,
                                                                                1000,0,4000, energy);
   }

 }
}

void Other(const GBCS& bcs, const std::vector<ddasHit>& hits){


//-----------------DSSD Position----------------------//

  if(bcs.EventType() == 1) {   // IMPLANT
  GHistogramer::Get().Fill("Position/DSSD_Implant", 40,0,40, bcs.fDSSD.Xpos,
                                                    40,0,40, bcs.fDSSD.Ypos); } 

  if(bcs.EventType() == 2) {   // DECAY
  GHistogramer::Get().Fill("Position/DSSD_Decay", 40,0,40, bcs.fDSSD.Xpos,
                                                  40,0,40, bcs.fDSSD.Ypos); }


}
