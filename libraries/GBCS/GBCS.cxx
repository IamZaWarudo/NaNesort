

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

 // DSSDLow.Reset();  // Not using this - unreliable
 DSSD.Reset();

 Pin1.Reset();
 Pin2.Reset();
  
 I2SPin1.Reset();
 I2SPin2.Reset();

 SSSDLow.Reset();
 SSSDHigh.Reset();
}

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
}


// fDSSDHigh filled in DSSD because we are going to use High only

//GBCS::Unpack(std::vector<ddasHit> event) { Unpacker::Unpack(*this,event); }

void GBCS::Fill(const Unpacker& Event) {
  
  //GBCS::Fill(const std::vector<ddasHit> event) { // this comes from the eventbuilder. 
  


 // DSSDLow = Event.fDSSDLow;
 DSSD = Event.fDSSDHigh;

 Pin1 = Event.fPin1;
 Pin2 = Event.fPin2;
 
 I2SPin1 = Event.fI2SPin1;
 I2SPin2 = Event.fI2SPin2;
 
 //I2NPin1 = Event.fI2NPin1;
 //I2NPin2 = Event.fI2NPin2;

 SSSDLow = Event.fSSSDLow;
 SSSDHigh = Event.fSSSDHigh;

 fCloverHits = Event.fCloverHits; 



/****************
Implant -----> 1
*****************
Decay -------> 2
*****************
Veto ----> 3
*****************
Undentified -> 4
*****************/

 

int condition = 4;  // Unidentified


bool hasPin      = false;
bool hasDSSD     = false;
bool hasDSSDGood = false;
bool hasSSSDL    = false;
bool hasSSSDH    = false;

bool hasI2       = false;

bool implant  = false;
bool veto     = false;
bool decay    = false;
bool lightion = false;



double dtPin = Pin1.fTimestamp - Pin2.fTimestamp;
if(Pin1.HasHit() && Pin2.HasHit() && dtPin > 5 && dtPin < 11) {
  hasPin = true;
}

if(I2SPin1.HasHit() || I2SPin2.HasHit()) {
  hasI2 = true;
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

if(SSSDHigh.HasHit()){
 hasSSSDH = true;
}



if(!hasPin && !hasI2 && hasDSSDGood && !hasSSSDL){
  decay = true;
  condition = 2; // decay 
}

if((hasPin || hasI2) && hasDSSD && (hasSSSDL || hasSSSDH)){
  veto = true;
  condition = 3; //veto
}

if(hasPin && hasI2 && hasDSSDGood & !hasSSSDL && !veto) {
  implant = true;
  condition = 1; // implant
}

/*
if(hasDSSD && hasSSSDH){
  lightion = true;
} */

  static int betterImplantCounter = 0;
  //if(condition==2) printf(DYELLOW "found decay" RESET_COLOR "\n");
  if(condition==1) printf(DGREEN  "found implant" RESET_COLOR "\n");
  //if(betterImplantCounter++>25) exit(0);


  fEventType =  condition;
}
