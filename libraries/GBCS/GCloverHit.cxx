#include <GCloverHit.h>

GCloverHit::GCloverHit()  { Reset(); }
GCloverHit::~GCloverHit() {}


void GCloverHit::Reset() {
  fId        = -1;
  fEcal      = 0;
  fTimestamp = 0;
  fCharge    = 0;
}

void GCloverHit::Unpack(const ddasHit& hit, int id) {

  fId = id;
  fCharge = hit.GetCharge();
  fTimestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();
}



