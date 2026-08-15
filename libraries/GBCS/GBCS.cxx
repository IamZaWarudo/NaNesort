

#include <ddasHit.h>
#include <Unpacker.h>
#include <GBCS.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>

GBCS::GBCS() { }

void GBCS::Reset() {

 // DSSDLow.Reset();  // Not using this - unreliable
 DSSD.Reset();

 Pin1.Reset();
 Pin2.Reset();

 SSSDLow.Reset();
 SSSDHigh.Reset();
  
}


// fDSSDHigh filled in DSSD because we are going to use High only

void GBCS::Fill(const Unpacker& Event) {
  

 // DSSDLow = Event.fDSSDLow;
 DSSD = Event.fDSSDHigh;

 Pin1 = Event.fPin1;
 Pin2 = Event.fPin2;
 I2SPin1 = Event.fI2SPin1;
 I2SPin2 = Event.fI2SPin2;


 SSSDLow = Event.fSSSDLow;
 SSSDHigh = Event.fSSSDHigh;

}


/****************
Implant -----> 1
*****************
decay -------> 2
*****************
LightIon ----> 3
*****************
Undentified -> 4
*****************/

 
int GBCS::EventType() {

int condition = 0101010;


bool hasPin      = false;
bool hasDSSD     = false;
bool hasDSSDGood = false;
bool hasSSSDL    = false;
bool hasSSSDH    = false;

bool implant  = false;
bool veto     = false;
bool lightion = false;



double dtPin = Pin1.fTimestamp - Pin2.fTimestamp;
if(Pin1.HasHit() && Pin2.HasHit() && dtPin > 5 && dtPin < 11) {
  hasPin = true;
}


if(DSSD.HasHit()){
  hasDSSD = true;
}

if(hasDSSD == true && DSSD.HasGoodPosition()){
  hasDSSDGood = true;
}

if(SSSDLow.HasHit()) {
 hasSSSDL = true;
}


  return condition;
}
