#include "nsmEvent.h"
ClassImp(nsmEvent)

nsmEvent::nsmEvent()
{
    vtx_x 		= 0.;
    vtx_y 		= 0.;
    vtx_z 		= 0.;
    Nprim 		= 0;
    run   		= 0;
    ref2  		= 0;
    btof 	= 0;
    event_plane = 0.;
}

nsmEvent::nsmEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
                   UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt_t mbtof, Float_t mevent_plane
                   )
{
    vtx_x = mvtx_x;
    vtx_y = mvtx_y;
    vtx_z = mvtx_z;
    
    Nprim = mNprim;
    run  = mrun;
    ref2 = mref2;
    
    btof = mbtof;

    event_plane = mevent_plane;
}

void nsmEvent::SetEventData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
                            UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt_t mbtof, FLoat_t mevent_plane
                            )
{
    vtx_x = mvtx_x;
    vtx_y = mvtx_y;
    vtx_z = mvtx_z;
    
    Nprim = mNprim;
    run  = mrun;
    ref2 = mref2;
    
    btof = mbtof;

    event_plane = mevent_plane;
}
