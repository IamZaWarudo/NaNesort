#ifndef GSSSD_H
#define GSSSD_H

#include <ddasHit.h>


class GSSSD {
  public:
  GSSSD();
  ~GSSSD();
 
  void Get(int strip, const ddasHit& hit);
  void Reset();

//  private:

  int Strip;
  double fCharge;
  double fEcal;
  double fTime;
  int fId;

  ClassDef(GSSSD,1);
};


#endif
