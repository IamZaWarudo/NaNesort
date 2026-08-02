
#include <GSSSD.h>

GSSSD::GSSSD() { Reset(); }
GSSSD::~GSSSD() {}

void GSSSD::Reset() {
  fCharge = 0;
  fEcal   = 0;
  fTime   = 0;
}

void GSSSD::Get(int strip, const ddasHit& hit) {
  
  fTime   = hit.GetTime();
  fCharge = hit.GetCharge();
  fEcal   = hit.GetEcal();

}
