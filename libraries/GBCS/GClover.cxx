#include <GClover.h>

GClover::GClover()  { Reset(); }
GClover::~GClover() {}


void GClover::Reset() {
  fId        = -1;
  fEcal      = 0;
  fTimestamp = 0;
  fCharge    = 0;
  fCloverHits.clear();
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

void GClover::LoadTalk(const char* filename) {  
  std::ifstream input(filename);
  int clover;

  for(int block = 0; block < 16; ++block) {
  for(int row = 0; row < 4; ++row) {
    input >> clover;

  for(int column = 0; column < 4; ++column)
    input >> Xtalk[clover][row][column];
  }
}

}

void GClover::AddBack(){

// Build individual clovers

// double Co-eff = XTalk(CloverID][i][j]

 for(const auto& c : fCloverHits){
  int ID = c.GetId() -224;
  if(ID < 0 || ID >= 16);  

  int CloverID = (c.GetId() - 224) / 4 ;
  int CrystalID = (c.GetId() - 224) % 4;




 }
}



