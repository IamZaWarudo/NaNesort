#ifndef GCLOVER_H
#define GCLOVER_H

#include <ddasHit.h>


class GClover {
  public:
  GClover();
  ~GClover();
  
  void Get(int Crystal, const ddasHit& hit);
  void Reset();
 
  private: 

  double fCharge;
  double fEcal;
  double fTime;

 ClassDef(GClover,1);
};


#endif
