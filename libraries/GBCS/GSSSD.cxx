#include <ddasHit.h>
#include <GSSSD.h>

GSSSD::GSSSD() { Reset();  }
GSSSD::~GSSSD() {}

void GSSSD::Reset() {
  fSum        = 0;
  fTimestamp  = 0;
  fStrips.clear();
  fEcal.clear();
}

void GSSSD::Unpack(const ddasHit& hit) {   
  const int Id = hit.GetId();
  if(fTimestamp==0) 
    fTimestamp = hit.GetTimestamp();
  else if(hit.GetTimestamp()<fTimestamp)
    fTimestamp = hit.GetTimestamp();
  
  fSum += hit.GetEcal();
  fEcal.push_back(hit.GetEcal());
  if(Id >= 160 && Id <= 175){
    fStrips.push_back(Id -160);
  }else if(Id >= 192 && Id <=207){
    fStrips.push_back(Id -192);
  }  
}

bool GSSSD::HasHit(){ 
  if(fSum>0) 
    return true;
  return false;
}

