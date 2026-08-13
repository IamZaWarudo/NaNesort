#include <ddasHit.h>
#include <GPIN.h>

GPIN::GPIN() { Reset();  }
GPIN::~GPIN() {}

void GPIN::Reset() {
  fEcal      = 0;
  Timestamp  = 0;
  fCharge    = 0;
  Hit     = false;
}

void GPIN::Unpack(const ddasHit& hit) {   
 
  fCharge = hit.GetCharge();
  Timestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();
  Hit = true;   
}

bool GPIN::HasHit(){ 
return Hit;
}

