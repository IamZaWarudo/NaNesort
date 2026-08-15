#include <ddasHit.h>
#include <GDSSD.h>
#include <Unpacker.h>
#include <cmath>

GDSSD::GDSSD() { Reset();  }
GDSSD::~GDSSD() {}

void GDSSD::Reset() {
  fId        = -1;
  fEcal      = 0;
  frontTimestamp  = 0;
  backTimestamp   = 0;
  xpos       = 0;
  ypos       = 0;
  Xpos       = 0;
  Ypos       = 0;
  GoodX      = false;
  GoodY      = false;

  FrontHit   = false;
  BackHit    = false;

  fTimestamp  = 0;

  fxEsum     = 0; 
  fxE        = 0; 
  fxMaxE     = 0; 
  fxMaxStrip = 0;
  fyEsum     = 0; 
  fyE        = 0; 
  fyMaxE     = 0; 
  fyMaxStrip = 0;
  
  fxMult     = 0;
  fyMult     = 0;
}

void GDSSD::UnpackFront(const ddasHit& hit) {   
// if (hit.GetEcal() >= 18000 ) return;

  fId = hit.GetId();
  if (frontTimestamp == 0 || hit.GetTimestamp() < frontTimestamp)
  frontTimestamp = hit.GetTimestamp();

  fEcal = hit.GetEcal();

   if(fId >= 0 && fId <= 39){
    xpos = fId;
   }else if(fId >= 40 && fId <= 79){
    xpos = fId - 40;
   }  

  fxEsum += xpos * fEcal;
  fxE    += fEcal;
  fxMult++;

  if(fEcal >= fxMaxE) {
    fxMaxE = fEcal;
    fxMaxStrip = xpos;}

  FrontHit = true;  
}

void GDSSD::UnpackBack(const ddasHit& hit) {   
// if (hit.GetEcal() >= 18000 ) return;

  fId = hit.GetId();
  if (backTimestamp == 0 || hit.GetTimestamp() < backTimestamp)
  backTimestamp = hit.GetTimestamp();

  fEcal = hit.GetEcal();

   if(fId >= 80 && fId <= 119){
    ypos = fId - 80;
   }else if(fId >= 120 && fId <= 159){
    ypos = fId - 120;
   }  
  
  fyEsum += ypos * fEcal;
  fyE    += fEcal;
  fyMult++;

  if(fEcal >= fyMaxE) {
    fyMaxE = fEcal;
    fyMaxStrip = ypos;}


  BackHit = true;
}


bool GDSSD::HasHit(){
 return FrontHit && BackHit;
} 


bool GDSSD::HasGoodPosition() {
  double dt = frontTimestamp - backTimestamp;

  fTimestamp = frontTimestamp;

  if (dt < 2 || dt > 14) {
    GoodX = false;
    GoodY = false;
    return false;
  }


  if (fxE > 0){
    Xpos = fxEsum / fxE;
    GoodX = fabs(Xpos - fxMaxStrip) <= 3;
  }
  if (fyE > 0){
    Ypos = fyEsum / fyE;
    GoodY = fabs(Ypos - fyMaxStrip) <= 3;
  }


return GoodX && GoodY;
}



