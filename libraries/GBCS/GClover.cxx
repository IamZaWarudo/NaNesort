#include <GClover.h>

GClover::GClover()  { Reset(); }
GClover::~GClover() {}


void GClover::Reset() {
  fId        = -1;
  fEcal      = 0;
  fTimestamp = 0;
  fCharge    = 0;
}

void GClover::Unpack(const ddasHit& hit, int id) {

  fId = id;
  fCharge = hit.GetCharge();
  fTimestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();
}

bool GClover::HasHit() const{
  if(fEcal >0)
    return true;
  return false;
}


