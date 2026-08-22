#ifndef GSSSD_H
#define GSSSD_H

#include <ddasHit.h>


class GSSSD {
  public:
  GSSSD();
  ~GSSSD();
 
  void Unpack(const ddasHit& hit,int id);
  void Reset();

//  private:

  bool HasHit() const;
  //bool Hit = false;

  double fSum;
  double fTimestamp;
  //int strip;
  std::vector<double> fStripTime;
  std::vector<int>    fStrips;
  std::vector<double> fEcal;

  ClassDef(GSSSD,2);
};


#endif
