#include <ddasHit.h>
#include <GPIN.h>

GPIN::GPIN() { Reset();  }
GPIN::~GPIN() {}

void GPIN::Reset() {
  fEcal      = 0;
  Timestamp  = 0;
  fCharge    = 0;
}

void GPIN::Unpack(const ddasHit& hit) {   
 
  fCharge = hit.GetCharge();
  Timestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();
}

bool GPIN::HasHit(){ 
  if(fEcal >0)
    return true;
  return false;
}

