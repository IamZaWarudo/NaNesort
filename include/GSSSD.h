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

  double fCharge;
  double fEcal;
  double fTime;


  ClassDef(GSSSD,1);
};


#endif
