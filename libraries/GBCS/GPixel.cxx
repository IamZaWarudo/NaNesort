#include <ddasHit.h>
#include <GPixel.h>
#include <GBCS.h>
#include <cmath>

GPixel::GPixel() { Reset();}
GPixel::~GPixel() {}

void GPixel::Reset() {
  ImPixel = {-1.0, -1.0};
  DePixel = {-1.0, -1.0};

  ImTime = -1.0;
  ImEnergy = 0.0;
  ImfMaxStrip = -1;
  ImbMaxStrip = -1;
  ImfMult = 0;
  ImbMult = 0;

  DeTime = -1.0;
  DeEnergy = 0.0;
  DefMaxStrip = -1;
  DebMaxStrip = -1;
  DefMult = 0;
  DebMult = 0;

}

void GPixel::BuildImplantPixel(const std::vector<ddasHit>& hits) {
  double feSum = 0;
  double feSumStrip = 0;
  double beSum = 0;
  double beSumStrip = 0;
  double fMaxE = -1.0;
  double bMaxE = -1.0;
  double fMaxTime = -1.0;
  double bMaxTime = -1.0;
  int fMult = 0;
  int bMult = 0;
  int fMaxStrip = -1;
  int bMaxStrip = -1;


  for(const auto& hit : hits){
    const int id = hit.GetId();
    const double e = hit.GetEcal();

    if(id >= 40 && id <= 79){
      const int strip = id - 40;
      feSum += e;
      feSumStrip += e * strip;
      ++fMult;
      if(e > fMaxE) {fMaxE = e; fMaxTime = hit.GetTime(); fMaxStrip = strip;}
    } else 
    if(id >= 120 && id <= 159) {
      const int strip = id - 120;
      beSum += e;
      beSumStrip += e * strip;
      ++bMult;
       if(e > bMaxE) {bMaxE = e; bMaxTime = hit.GetTime(); bMaxStrip = strip;}     
       }
  }

  const double dt = fMaxTime - bMaxTime;
  if(dt < 2 || dt > 7) return;

  ImPixel[0] = feSum > 0.0 ? feSumStrip/feSum : -1;
  ImPixel[1] = beSum > 0.0 ? beSumStrip/beSum : -1;

  ImTime = fMaxTime;
  ImEnergy = 0.5 * (feSum + beSum);
  ImfMaxStrip = fMaxStrip;
  ImbMaxStrip = bMaxStrip;
  ImfMult = fMult;
  ImbMult = bMult;
}


void GPixel::BuildDecayPixel(const std::vector<ddasHit>& hits) {
  double feSum = 0;
  double feSumStrip = 0;
  double beSum = 0;
  double beSumStrip = 0;
  double fMaxE = -1.0;
  double bMaxE = -1.0;
  double fMaxTime = -1.0;
  double bMaxTime = -1.0;
  int fMult = 0;
  int bMult = 0;
  int fMaxStrip = -1;
  int bMaxStrip = -1;

  for(const auto& hit : hits){
    const int id = hit.GetId();
    const double e = hit.GetEcal();

    if(id >= 0 && id <= 39){
      const int strip = id;
      feSum += e;
      feSumStrip += e * strip;
      ++fMult;
      if(e > fMaxE) {fMaxE = e; fMaxTime = hit.GetTime(); fMaxStrip = strip;}
    }else 
    if(id >= 80 && id <= 119) {
      const int strip = id - 80;
      beSum += e;
      beSumStrip += e * strip;
      ++bMult;
       if(e > bMaxE) {bMaxE = e; bMaxTime = hit.GetTime(); bMaxStrip = strip;}     
       }

     }

  const double dt = fMaxTime - bMaxTime;
  if(dt < 4 || dt > 12) return;

  DePixel[0] = feSum > 0.0 ? feSumStrip/feSum : -1;
  DePixel[1] = beSum > 0.0 ? beSumStrip/beSum : -1;

  DeTime = fMaxTime;
  DeEnergy = 0.5 * (feSum + beSum);
  DefMaxStrip = fMaxStrip;
  DebMaxStrip = bMaxStrip;
  DefMult = fMult;
  DebMult = bMult;
}

