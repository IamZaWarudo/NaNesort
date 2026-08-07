
#include <GDSSD.h>

GDSSD::GDSSD() { Reset();  }
GDSSD::~GDSSD() {}

void GDSSD::Reset() {
  fCharge = 0;
  fEcal   = 0;
  fTime   = 0;
  Strip   = 0;
}

void GDSSD::Get(int strip, const ddasHit& hit) {
  
  Strip   = strip;
  fTime   = hit.GetTime();
  fCharge = hit.GetCharge();
  fEcal   = hit.GetEcal();

}
