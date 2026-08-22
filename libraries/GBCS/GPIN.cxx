#include <ddasHit.h>
#include <GPIN.h>

GPIN::GPIN() { Reset();  }
GPIN::~GPIN() {}

void GPIN::Reset() {
  fEcal       = 0;
  fTimestamp  = 0;
  fCharge     = 0;
}

void GPIN::Unpack(const ddasHit& hit) {   
 
  fCharge = hit.GetCharge();
  fTimestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();
}

bool GPIN::HasHit() const{ 
  if(fEcal >0)
    return true;
  return false;
}

