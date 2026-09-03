#include <GClover.h>

GClover::GClover()  { Reset(); }
GClover::~GClover() {}


void GClover::Reset() {
  fId        = -1;
  fEcal      = 0;
  fTimestamp = 0;
  fCharge    = 0;
  fCloverHits.clear();
}

void GClover::Unpack(const ddasHit& hit) {

  fId = hit.GetId() - 224;
  fCharge = hit.GetCharge();
  fTimestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();  // overwritten every loop
}

bool GClover::HasHit() const{
  if(fEcal >0)
    return true;
  return false;
}

void GClover::AddBack(){


 for(const auto& c : fCloverHits){
  
   // Do Something

  int AddBackClover = c.GetId() / 4 ;
  if(AddBackClover < 0 || AddBackClover >= 16) continue;


 }


}
