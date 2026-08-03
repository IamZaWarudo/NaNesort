#include <ddasHit.h>
#include <GPixel.h>

GPixel::GPixel() { Reset(); }
GPixel::~GPixel() {}


void GPixel::Reset() {
  ImplantPixel = {-1.0, -1.0};
  DecayPixel   = {-1.0, -1.0};

  ImplantFrontStrip = -1;  
  ImplantBackStrip  = -1;
  DecayFrontStrip   = -1;  
  DecayBackStrip    = -1;

  ImplantTime   = -1.0;  
  DecayTime     = -1.0;
  ImplantEnergy =  0.0;  
  DecayEnergy   =  0.0;

  ImplantFrontE = 0.0;  
  ImplantBackE  = 0.0;
  DecayFrontE   = 0.0;  
  DecayBackE    = 0.0;

  ImplantFrontMult = 0;  
  ImplantBackMult  = 0;
  DecayFrontMult   = 0;  
  DecayBackMult    = 0;
}


// Implant - front low 40...79 & back low 120...159
void GPixel::BuildImplantPixel(const std::vector<ddasHit>& hits) {
  double frontSumE = 0.0, frontSumEStrip = 0.0;
  double backSumE  = 0.0, backSumEStrip  = 0.0;

  double frontMaxE = -1, frontMaxT = -1;
  double backMaxE  = -1, backMaxT  = -1;
  int    frontMaxStrip = -1, backMaxStrip = -1;
  int    frontMult = 0, backMult = 0;

  for(const auto& hit : hits) {
    const int    id     = hit.GetId();
    const double energy = hit.GetEcal();
    if(energy <= 0.0) continue;

    if(id >= 40 && id <= 79) {                 // front low gain
      const int strip = id - 40;
      frontSumE      += energy;
      frontSumEStrip += energy * strip;
      ++frontMult;
      if(energy > frontMaxE) {
        frontMaxE = energy; frontMaxT = hit.GetTime(); frontMaxStrip = strip;
      }
    }
    else if(id >= 120 && id <= 159) {          // back low gain
      const int strip = id - 120;
      backSumE      += energy;
      backSumEStrip += energy * strip;
      ++backMult;
      if(energy > backMaxE) {
        backMaxE = energy; backMaxT = hit.GetTime(); backMaxStrip = strip;
      }
    }
  }

  const double dt = frontMaxT - backMaxT;
  if(dt < 2 || dt > 7) return;                 // LG: 2-7 ticks

  ImplantPixel[0] = frontSumE > 0.0 ? frontSumEStrip / frontSumE : -1.0;
  ImplantPixel[1] = backSumE  > 0.0 ? backSumEStrip  / backSumE  : -1.0;

  ImplantTime       = frontMaxT;
  ImplantEnergy     = 0.5 * (frontSumE + backSumE);
  ImplantFrontStrip = frontMaxStrip;
  ImplantBackStrip  = backMaxStrip;
  ImplantFrontE     = frontSumE;
  ImplantBackE      = backSumE;
  ImplantFrontMult  = frontMult;
  ImplantBackMult   = backMult;
}

// Decay - front high 0...39 & back high 80...119
void GPixel::BuildDecayPixel(const std::vector<ddasHit>& hits) {
  double frontSumE = 0.0, frontSumEStrip = 0.0;
  double backSumE  = 0.0, backSumEStrip  = 0.0;

  double frontMaxE = -1, frontMaxT = -1;
  double backMaxE  = -1, backMaxT  = -1;
  int    frontMaxStrip = -1, backMaxStrip = -1;
  int    frontMult = 0, backMult = 0;

  for(const auto& hit : hits) {
    const int    id     = hit.GetId();
    const double energy = hit.GetEcal();
    if(energy <= 0.0) continue;

    if(id >= 0 && id <= 39) {                  // front high gain
      const int strip = id;
      frontSumE      += energy;
      frontSumEStrip += energy * strip;
      ++frontMult;
      if(energy > frontMaxE) {
        frontMaxE = energy; frontMaxT = hit.GetTime(); frontMaxStrip = strip;
      }
    }
    else if(id >= 80 && id <= 119) {           // back high gain
      const int strip = id - 80;
      backSumE      += energy;
      backSumEStrip += energy * strip;
      ++backMult;
      if(energy > backMaxE) {
        backMaxE = energy; backMaxT = hit.GetTime(); backMaxStrip = strip;
      }
    }
  }

  const double dt = frontMaxT - backMaxT;
  if(dt < 4 || dt > 13) return;                // HG: 4-13 ticks

  DecayPixel[0] = frontSumE > 0.0 ? frontSumEStrip / frontSumE : -1.0;
  DecayPixel[1] = backSumE  > 0.0 ? backSumEStrip  / backSumE  : -1.0;

  DecayTime       = frontMaxT;
  DecayEnergy     = 0.5 * (frontSumE + backSumE);
  DecayFrontStrip = frontMaxStrip;
  DecayBackStrip  = backMaxStrip;
  DecayFrontE     = frontSumE;
  DecayBackE      = backSumE;
  DecayFrontMult  = frontMult;
  DecayBackMult   = backMult;
}
