#include <ddasHit.h>
#include <GDSSD.h>
#include <Unpacker.h>


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

  fxEsum     = 0; 
  fxE        = 0; 
  fxMaxE     = 0; 
  fxMaxStrip = 0;
  fyEsum     = 0; 
  fyE        = 0; 
  fyMaxE     = 0; 
  fyMaxStrip = 0;

}

void GDSSD::UnpackFront(const ddasHit& hit) {   
  fId = hit.GetId();
  frontTimestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();

   if(fId >= 0 && fId <= 39){
    xpos = fId;
   }else if(fId >= 40 && fId <= 79){
    xpos = fId - 40;
   }  

  fxEsum += xpos * fEcal;
  fxE    += fEcal;

  if(fEcal >= fxMaxE) {
    fxMaxE = fEcal;
    fxMaxStrip = xpos;}

}

void GDSSD::UnpackBack(const ddasHit& hit) {   
  fId = hit.GetId();
  backTimestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();

   if(fId >= 80 && fId <= 119){
    ypos = fId - 80;
   }else if(fId >= 120 && fId <= 159){
    ypos = fId - 120;
   }  
  
  fyEsum += ypos * fEcal;
  fyE    += fEcal;

  if(fEcal >= fyMaxE) {
    fyMaxE = fEcal;
    fyMaxStrip = ypos;}

}

bool GDSSD::HasGoodPosition() {
  double dt = frontTimestamp - backTimestamp;

  if (dt < 2 || dt > 14) {
    GoodX = false;
    GoodY = false;
    return false;
  }

  if (fxE > 0){
    Xpos = fxEsum / fxE;
    GoodX = fabs(Xpos - fxMaxStrip) <= 2;
  }
  if (fyE > 0){
    Ypos = fyEsum / fyE;
    GoodY = fabs(Ypos - fyMaxStrip) <= 2;
  }

//printf("[GDSSD::Finalize] Xpos=%d (sum=%.2f, weight=%.2f, maxStrip=%.2f) fGoodX=%d | "
//       "Ypos=%d (sum=%.2f, weight=%.2f, maxStrip=%.2f) fGoodY=%d\n",
//       Xpos, fxEsum, fxE, fxMaxStrip, GoodX,
//       Ypos, fyEsum, fyE, fyMaxStrip, GoodY);

return GoodX && GoodY;
}



