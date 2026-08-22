#include <ddasHit.h>
#include <GDSSD.h>
#include <Unpacker.h>
#include <cmath>

GDSSD::GDSSD() { Reset();  }
GDSSD::~GDSSD() {}

void GDSSD::Reset() {

  GoodPosition = false;
  GoodTimeGate = false;
  GoodposGate  = false;

  fId    = -1;
  xpos   = -1;
  ypos   = -1;
  
  fEcal   = 0;
  fEnergy = 0;

  fxWEsum = 0;
  fyWEsum = 0;
  
  fxEsum  = 0; 
  fyEsum  = 0;
  fxMaxE  = 0;
  fyMaxE  = 0; 
  
  fxMaxStrip = -1;
  fyMaxStrip = -1;
  
  FrontTimestamp = 0;
  BackTimestamp  = 0;
  fTimestamp     = 0;
  
  XStripDeviation = 0;
  YStripDeviation = 0;
  
}

void GDSSD::UnpackFront(const ddasHit& hit) {   

  fId = hit.GetId();
  fEcal = hit.GetEcal();

  if(fId >= 0 && fId < 39){
    xpos = fId;
  }else
  if(fId >= 40 && fId <= 79){
    xpos = fId - 40;
  }

  fxWEsum += xpos * fEcal;
  fxEsum  += fEcal;
   if(fxEsum <= 0) return;  // Energy gate

  Xpos = fxWEsum / fxEsum;  // final X co-ord

  if(fEcal >= fxMaxE){
    fxMaxE = fEcal;
    fxMaxStrip = xpos;
    FrontTimestamp = hit.GetTimestamp();  
  }

  XStripDeviation = std::abs(Xpos - fxMaxStrip);

}

void GDSSD::UnpackBack(const ddasHit& hit) {   

  fId = hit.GetId();
  fEcal = hit.GetEcal();

  if(fId >= 80 && fId < 119){
    ypos = fId - 80;
  }else
  if(fId >= 120 && fId <= 159){
    ypos = fId - 120;
  }

  fyWEsum += ypos * fEcal;
  fyEsum  += fEcal;
   if(fyEsum <= 0) return;  // Energy gate

  Ypos = fyWEsum / fyEsum;  // final Y co-ord

  if(fEcal >= fyMaxE){
    fyMaxE = fEcal;
    fyMaxStrip = ypos;
    BackTimestamp = hit.GetTimestamp(); 
  }

  YStripDeviation = std::abs(Ypos - fyMaxStrip);

}



bool GDSSD::HasPosition() {
  GoodPosition = false;

  fTimestamp = FrontTimestamp; // front Strip - Max Energy used for Time ( consistent across all )

  double dt = FrontTimestamp - BackTimestamp;
  if( dt > 2 && dt < 14) {
    GoodTimeGate = true;
  }

  double drift = 3;
  if(XStripDeviation < drift && YStripDeviation < drift){
    GoodposGate = true;
  }
 
  if(GoodTimeGate && GoodposGate){
    GoodPosition = true;
  }

  fEnergy = (fxEsum + fyEsum) * 0.5;  // DSSD Energy for comparison/ploting purposes 


 return GoodPosition; 
}



