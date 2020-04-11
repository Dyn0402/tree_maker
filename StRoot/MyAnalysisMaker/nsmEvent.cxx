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
    ref3        = 0;
    btof 	= 0;
    event_plane_ref2 = 0.;
    event_plane_ref3 = 0.;
}

nsmEvent::nsmEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
                   UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt_t mref3, UInt_t mbtof,
				   Float_t mevent_plane_ref2, Float_t mevent_plane_ref3)
{
    vtx_x = mvtx_x;
    vtx_y = mvtx_y;
    vtx_z = mvtx_z;
    
    Nprim = mNprim;
    run  = mrun;

    ref2 = mref2;
    ref3 = mref3;
    
    btof = mbtof;

    event_plane_ref2 = mevent_plane_ref2;
    event_plane_ref3 = mevent_plane_ref3;
}

void nsmEvent::SetEventData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z,
                            UInt_t mNprim, UInt_t mrun, UInt_t mref2, UInt_t mref3, UInt_t mbtof,
							Float_t mevent_plane_ref2, Float_t mevent_plane_ref3
                            )
{
    vtx_x = mvtx_x;
    vtx_y = mvtx_y;
    vtx_z = mvtx_z;
    
    Nprim = mNprim;
    run  = mrun;

    ref2 = mref2;
    ref3 = mref3;
    
    btof = mbtof;

    event_plane_ref2 = mevent_plane_ref2;
    event_plane_ref3 = mevent_plane_ref3;
}
