#include <ddasHit.h>
#include <GSSSD.h>

GSSSD::GSSSD() { Reset();  }
GSSSD::~GSSSD() {}

void GSSSD::Reset() {
  fSum        = 0;
  fTimestamp  = 0;
  fStripTime.clear();
  fStrips.clear();
  fEcal.clear();
}

void GSSSD::Unpack(const ddasHit& hit,int Id) {   
  if(fTimestamp==0) 
    fTimestamp = hit.GetTimestamp();
  else if(hit.GetTimestamp()<fTimestamp)
    fTimestamp = hit.GetTimestamp();
  fSum += hit.GetEcal();
  
  fStripTime.push_back(hit.GetTimestamp());
  fEcal.push_back(hit.GetEcal());
  fStrips.push_back(Id);

  //if(Id >= 160 && Id <= 175){
  //  fStrips.push_back(Id -160);
  //}else if(Id >= 192 && Id <=207){
  //  fStrips.push_back(Id -192);
  //}  
}

bool GSSSD::HasHit(){ 
  if(fSum>0) 
    return true;
  return false;
}

