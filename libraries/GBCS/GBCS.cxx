#include <ddasHit.h>
#include <Unpacker.h>
#include <GBCS.h>
#include <GPIN.h>
#include <GDSSD.h>
#include <GSSSD.h>
#include <GLaBr.h>
#include <GClover.h>

#include <globals.h>

GBCS::GBCS() { }

void GBCS::Reset() {

 fEventType = 5;

 fDSSD.Reset(); //------> Used for Position
 fDSSDLow.Reset();

 fPin1.Reset(); //------> Used for TOF
 fPin2.Reset();
  
 fI2SPin1.Reset(); //---> Used for TOF
 fI2SPin2.Reset();

 fSSSDLow.Reset();
 fSSSDHigh.Reset();

 fClover.Reset();
}
 
/*

 +----------------+
 | 1 - IMPLANT    |
 | 2 - DECAY      |
 | 3 - LIGHTIONS  |
 | 4 - VETO       |
 | 5 - UNKNOWN    |
 +----------------+

*/

void GBCS::Classification() {

  fEventType = 5;  // Unknown



bool hasGoodPosition  = false;
bool hasBadPosition   = false;
bool hasDSSD          = false;
bool hasSSSDL         = false;
bool hasSSSDH         = false;

bool hasPin           = false;
bool hasNoPin         = false;
bool hasTOF           = false;

bool PinLightIon      = false;
bool SSSDHLightIon    = false;

bool IsImplant        = false; 
bool IsDecay          = false;
bool IsLightIon       = false;
bool IsVeto           = false;

bool TOFGate          = false;



if(fPin1.HasHit() == true && fPin2.HasHit() == true) {
  hasPin = true;
}

if(fPin1.HasHit() == false && fPin2.HasHit() == false && fI2SPin1.HasHit() == false) {
  hasNoPin = true;
}


if(fPin1.HasHit() == true && fI2SPin1.HasHit() == true) {
  hasTOF = true;
}


if(fDSSD.HasPosition() == true) { 
  hasGoodPosition = true; 
}else
if(fDSSD.HasPosition() == false) {
  hasBadPosition  = true;
}

if(fSSSDLow.HasHit() == true) {
  hasSSSDL = true;
}

if(fSSSDHigh.HasHit() == true) {
  hasSSSDH = true;
}

if(hasPin && fPin1.fCharge < 3500) {
  PinLightIon = true;
}

if(fSSSDHigh.fSum > 2600) {
  SSSDHLightIon = true;
}


if(PinLightIon || SSSDHLightIon) {
  IsLightIon = true;   //-------------> LightIons
  fEventType = 3;
}

if(fDSSD.HasPosition()) {
  hasDSSD = true;
}


if(hasPin && hasDSSD && hasSSSDL) {
  IsVeto = true;   //-----------------> Veto
  fEventType = 4;
}


if(fI2SPin1.fCharge > 8000 && fI2SPin1.fCharge < 18000) {
  TOFGate = true;
}

if(hasPin && hasTOF && hasGoodPosition && !hasSSSDL && !IsLightIon && TOFGate) {
  IsImplant = true; //----------------> Implant
  fEventType = 1;
}

if(hasNoPin && !hasTOF && hasGoodPosition && !hasSSSDL && !IsLightIon) {
  IsDecay = true;
  fEventType = 2;
}



}




/*
void GBCS::Print(int type) const {
  switch(fEventType) {
    case 1: //implant
      printf(DRED);
      break;
    case 2: //decay
      printf(DGREEN);
      break;
    case 3: //veto
      printf(DYELLOW);
      break;
    case 4: //unknown
    default:
      printf(CYAN);
      break;
  }
  
  if(type<0 || type==fEventType) {
    printf("BCS [%02i][%02i] @ %.1f\n",X(),Y(),Timestamp());
    printf("EvnetType = %i\n",fEventType);
    printf("\tPin1: %f\n",Pin1.fEcal);
    printf("\tTOF:  %f\n",TOF());
    printf("\tDSSD MULTI [%i][%i]\n",DSSD.fxMult,DSSD.fyMult);
    printf("\tDSSD charege %.1f\n",DSSD.fxE + DSSD.fyE);
    printf("\tSSSDLow MULTI [%lu]\n",SSSDLow.fStrips.size());
    printf("\tSSSDHigh MULTI [%lu]\n",SSSDHigh.fStrips.size());
  }
  printf(RESET_COLOR);
  return;
}   */
