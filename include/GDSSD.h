#ifndef GDSSD_H
#define GDSSD_H

#include <ddasHit.h>

class GDSSD {
  public:
  GDSSD();
  ~GDSSD();
  
  void Get(int strip, const ddasHit& hit);
  void Reset();


//  private:

  double fCharge;
  double fEcal;
  double fTime;

  ClassDef(GDSSD,1);
};


#endif
