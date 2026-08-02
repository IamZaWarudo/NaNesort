

#include <GClover.h>

GClover::GClover() { Reset();  }
GClover::~GClover() {}

void GClover::Reset() {
  fCharge = 0;
  fEcal   = 0;
  fTime   = 0;
}


void GClover::Get(int Crystal, const ddasHit& hit) {
  
  fTime   = hit.GetTime();
  fCharge = hit.GetCharge();
  fEcal   = hit.GetEcal();

}
