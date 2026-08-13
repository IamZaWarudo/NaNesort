#ifndef GPIN_H
#define GPIN_H

#include <ddasHit.h>

class GPIN {
  public:
  GPIN();
  ~GPIN();
  
  void Unpack(const ddasHit& hit);
  void Reset();

  bool HasHit();
  bool Hit = false;

  double fCharge;
  double fEcal;
  double Timestamp;

  ClassDef(GPIN,1);
};


#endif
