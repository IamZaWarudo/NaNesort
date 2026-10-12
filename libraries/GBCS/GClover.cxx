#include <GClover.h>
#include <fstream>


GClover::GClover()  { Reset(); }
GClover::~GClover() {}


void GClover::Reset() {
  fId        = -1;
  fEcal      = 0;
  fTimestamp = 0;
  fCharge    = 0;
  fCloverHits.clear();

  for(int clover = 0; clover < 16; ++clover){
    fAddback[clover] = 0.0;
  for (int crystal = 0; crystal < 4; ++ crystal){
    CrystalEcal[clover][crystal] = 0.0;
    EcalCorrection[clover][crystal] = 0.0;
    ActiveCrystal[clover][crystal] = false;
   }
  }

}


void GClover::Unpack(const ddasHit& hit) {

  fId = hit.GetId() - 224;
  fCharge = hit.GetCharge();
  fTimestamp = hit.GetTimestamp();
  fEcal = hit.GetEcal();  // overwritten every loop
  
}

bool GClover::HasHit() const{
  if(fEcal >0)
    return true;
  return false;
}

bool GClover::LoadXTalk(const char* filename) {  
  std::ifstream input(filename);
  if(!input) return false;

// read the Co-effs
  double CoEffs[16][4][4] = {};
  for(int block = 0; block < 16; ++block){
  
  for(int row = 0; row < 4; ++row){
    int clover;
    if(!(input >> clover) || clover != block) return false;
  
  for(int column = 0; column < 4; ++column){
    if (!(input >> CoEffs[clover][row][column])) return false;
  }}}

// fill the Co-effs in XTalk
  for(int clover = 0; clover < 16; ++clover)
  for(int row = 0; row < 4; ++row)
  for(int column = 0; column < 4; ++column)
    XTalk[clover][row][column] = CoEffs[clover][row][column];

  return true;
}

void GClover::AddBack(){
  

  for(int clover = 0; clover < 16; ++clover){
    fAddback[clover] = 0.0;
  for (int crystal = 0; crystal < 4; ++ crystal){
    CrystalEcal[clover][crystal] = 0.0;
    EcalCorrection[clover][crystal] = 0.0;
    ActiveCrystal[clover][crystal] = false;
   }
  }


 for(const auto& hit : fCloverHits){
  int id = hit.GetId() - 224;
  if(id < 0 || id >=64) continue;

  int clover = id / 4;
  int crystal = id % 4;
// might need to add a time gate for consecutive hits
  CrystalEcal[clover][crystal] += hit.GetEcal();
  ActiveCrystal[clover][crystal] = true;
 }

 for(int clover = 0; clover < 16; ++clover){
   for(int i = 0; i < 4; ++i){
    if(!ActiveCrystal[clover][i]) continue;

    double GammaBaby = CrystalEcal[clover][i];

    for(int j = 0; j < 4; ++j){
      GammaBaby += XTalk[clover][i][j] * CrystalEcal[clover][j];
    }

    EcalCorrection[clover][i] = GammaBaby;
    fAddback[clover] += GammaBaby;
   }
 }

}



