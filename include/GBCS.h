#ifndef __GBCS_H__
#define __GBCS_H__

#include <vector>
#include <ddasHit.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>


class GBCS {
 public:
  GBCS();
  ~GBCS() = default;

  void Reset();

  void Classification();
  int EventType() const { return fEventType; }
  int fEventType = 5;

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

  GClover fClover;
  GLaBr fLaBr; 


};




#endif
