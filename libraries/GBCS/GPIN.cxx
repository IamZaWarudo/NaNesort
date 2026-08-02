
#include <GPIN.h>

GPIN::GPIN() { Reset(); }
GPIN::~GPIN() {}

void GPIN::Reset() {
  fCharge = 0;
  fEcal   = 0;
  fTime   = 0;
}

void GPIN::Get(const ddasHit& hit){
  
  fTime   = hit.GetTime();
  fCharge = hit.GetCharge();
  fEcal   = hit.GetEcal();
}
