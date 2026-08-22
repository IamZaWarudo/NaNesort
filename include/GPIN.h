#ifndef GPIN_H
#define GPIN_H

#include <ddasHit.h>

class GPIN {
  public:
  GPIN();
  ~GPIN();
  
  void Unpack(const ddasHit& hit);
  void Reset();

  bool HasHit() const;

 //private
  double fCharge;
  double fEcal;
  double fTimestamp;

  ClassDef(GPIN,1);
};


#endif
