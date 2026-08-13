#include <ddasHit.h>
#include <GSSSD.h>

GSSSD::GSSSD() { Reset();  }
GSSSD::~GSSSD() {}

void GSSSD::Reset() {
  fEcal      = 0;
  Timestamp  = 0;
  strip      = 0;
  Hit     = false;
}

void GSSSD::Unpack(const ddasHit& hit) {   
  int Id = hit.GetId();
  Timestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();
    
  if(Id >= 160 && Id <= 175){
    strip = Id -160;
  }else if(Id >= 192 && Id <=207){
    strip = Id -192;
    }  
  Hit = true;
}

bool GSSSD::HasHit(){ 
return Hit && fEcal > 500;
}

