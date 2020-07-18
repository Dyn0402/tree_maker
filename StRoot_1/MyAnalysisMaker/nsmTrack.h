#ifndef ROOT_nsmTrack
#define ROOT_nsmTrack

#include "TObject.h"

class nsmTrack : public TObject
{
private:

  Float_t pt;
  Float_t phi;
  Float_t eta;
  Float_t dca;
  Float_t nsigma;
  Float_t beta;
  Int_t charge;

public:
 nsmTrack();
 nsmTrack(   Float_t mpt, 
	     Float_t mphi, 
	     Float_t meta, 
	     Float_t mdca,
         Float_t mnsigma,
         Float_t mbeta,
         Int_t mcharge
	     );


 Float_t Pt() { return pt; }
 Float_t Phi() { return phi; }
 Float_t Eta() { return eta; }
 Float_t Dca() { return dca; }
 Float_t Nsigma() {return nsigma; }
 Float_t Beta() { return beta; }
 Int_t Charge() {return charge;}
 ClassDef(nsmTrack,2) 
};

#endif
