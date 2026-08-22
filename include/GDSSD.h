#ifndef GDSSD_H
#define GDSSD_H

#include <ddasHit.h>

class GDSSD {
  public:
  GDSSD();
  ~GDSSD();
  
  void UnpackFront(const ddasHit& hit);
  void UnpackBack(const ddasHit& hit);
  
  void Reset();

  void Finalize();
  bool HasPosition() const { return GoodPosition; }

//  private:
 
  bool GoodPosition;
  bool GoodTimeGate;
  bool GoodposGate;

  double fId;
  double xpos;
  double ypos;
  
  double Xpos;
  double Ypos;

  double fEcal;
  double fEnergy;

  double fxWEsum;
  double fyWEsum;
   
  double fxEsum; 
  double fyEsum;
  double fxMaxE;
  double fyMaxE; 
   
  double fxMaxStrip;
  double fyMaxStrip;
   
  double FrontTimestamp;
  double BackTimestamp;
  double fTimestamp;
   
  double XStripDeviation;
  double YStripDeviation;
 

 ClassDef(GDSSD,1);
};


#endif
