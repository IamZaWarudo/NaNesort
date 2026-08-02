#ifndef GPIN_H
#define GPIN_H

#include <ddasHit.h>

class GPIN {
  public:
  GPIN();
  ~GPIN();
  
  void Get(const ddasHit& hit);
  void Reset();

  double fCharge;
  double fEcal;
  double fTime;

  ClassDef(GPIN,1);
};


#endif
