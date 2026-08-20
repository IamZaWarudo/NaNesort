#ifndef __GBCS_H__
#define __GBCS_H__

#include <vector>
#include <ddasHit.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>
#include <GCloverHit.h>


class GBCS {
 public:
  GBCS();
  ~GBCS() = default;

  void Reset();


  //GDSSD DSSDLow;

  double Timestamp() const { return fDSSD.fTimestamp; }
  double TOF() const { return fI2SPin1.fEcal>0 ? fI2SPin1.fEcal : -1; } 
  int    X()   const { return fDSSD.Xpos; }
  int    Y()   const { return fDSSD.Ypos; }

  void Print(int type=-1) const;
  int EventType() const;


//private:
  
  GDSSD fDSSD; // DSSD ->> is DSSD High Gain
  GDSSD fDSSDLow;

  GPIN fPin1;
  GPIN fPin2;
  GPIN fI2TAC;
  GPIN fI2SPin1;
  GPIN fI2SPin2;

  GSSSD fSSSDLow;
  GSSSD fSSSDHigh;

  GClover fCrystal;
  std::vector<GClover> fCloverHits;


  int fEventType;


};




#endif
