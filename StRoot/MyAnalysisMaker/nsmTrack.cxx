#include "nsmTrack.h"
ClassImp(nsmTrack)

nsmTrack::nsmTrack()
{

   pt = 0.;
   p = 0.;
   phi = 0.;
   eta = 0.;
   dca = 0.;
   nsigmapr=0.;
   nsigmapi=0.;
   beta=0.;
   charge=0.;
    
}

nsmTrack::nsmTrack(
			 Float_t mpt, 
			 Float_t mp,
			 Float_t mphi, 
			 Float_t meta, 
			 Float_t mdca,
			 Float_t mnsigmapr,
			 Float_t mnsigmapi,
			 Float_t mbeta,
             Float_t mcharge
	               
			 ) {
 
  pt = mpt;
  p = mp;
  phi = mphi;
  eta = meta;
  dca = mdca;
  nsigmapr=mnsigmapr;
  nsigmapi=mnsigmapi;
  beta=mbeta;
  charge=mcharge;

  
}
