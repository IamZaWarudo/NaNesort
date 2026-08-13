

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

 DSSDLow.Reset();
 DSSDHigh.Reset();

 Pin1.Reset();
 Pin2.Reset();

 SSSDLow.Reset();
 SSSDHigh.Reset();
  
}


void GBCS::Fill(const Unpacker& Event) {
  

 DSSDLow = Event.fDSSDLow;
 DSSDHigh = Event.fDSSDHigh;

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

int GBCS::EventType(){

 int condition = 123456789;

//signal flags
 bool DSSDLowHasGoodPosition = false;
 bool DSSDHighHasGoodPosition = false;
 
 bool Pin1HasHit = false;
 bool Pin2HasHit = false;

 bool SSSDLowHasHit = false;
 bool SSSDHighHasHit = false;

 bool HasToF = false;

 bool HasClover = false;
 bool HasLaBr = false;


//event identification
 bool IsImplant = false;        // 1
 bool IsDecay = false;          // 2
 bool IsLightIon = false;       // 3  
 bool IsUnidentified = false;   // 4


// add ifs here
 
  if(DSSDLow.HasGoodPosition()){
   DSSDLowHasGoodPosition = true;
  }
  
  if(DSSDHigh.HasGoodPosition()){
   DSSDHighHasGoodPosition = true;
  }

  if(Pin1.HasHit()){
    Pin1HasHit = true;
  }
  
  if(Pin2.HasHit()){
    Pin2HasHit = true;
  }

  if(I2SPin1.HasHit() && Pin1.HasHit()){
    HasToF = true;
  }

  if(SSSDLow.HasHit()){
    SSSDLowHasHit = true;
  }
  
  if(SSSDHigh.HasHit()){
    SSSDHighHasHit = true;
  }


if(HasToF && Pin1HasHit && Pin2HasHit && DSSDLowHasGoodPosition && !SSSDLowHasHit && !SSSDHighHasHit){
 IsImplant = true;
 condition = 1;
} else
if(!Pin1HasHit && !Pin2HasHit && !DSSDLowHasGoodPosition && DSSDHighHasGoodPosition && !SSSDLowHasHit && !SSSDHighHasHit){
 IsDecay = true;
 condition = 2;
} else
if(Pin1HasHit && Pin2HasHit && !DSSDLowHasGoodPosition && DSSDHighHasGoodPosition && !SSSDLowHasHit && SSSDHighHasHit){
 IsLightIon = true;
 condition = 3;
}

if(!IsImplant && !IsDecay && !IsLightIon) { 
 IsUnidentified = true;
 condition = 4;
}


// PID Diagnostic 

// Detectors : PIN1 PIN2 DSSDHigh SSSDLow SSSDHigh     (DSSDLow = unreliable)

// noise 
if(Pin1HasHit && Pin2HasHit && !DSSDHighHasGoodPosition && !SSSDLowHasHit && !SSSDHighHasHit) {
  condition = 5;}

// possible implants
if(Pin1HasHit && Pin2HasHit && DSSDHighHasGoodPosition && !SSSDLowHasHit && !SSSDHighHasHit) {
  condition = 6;}

// possible light ions
if(Pin1HasHit && Pin2HasHit && DSSDHighHasGoodPosition && !SSSDLowHasHit && SSSDHighHasHit) {
  condition = 7;}

  return condition;
}
