#include "nsmEvent.h"
ClassImp(nsmEvent)

nsmEvent::nsmEvent()
{
    vtx_x			= 0.;
    vtx_y			= 0.;
    vtx_z			= 0.;
    dca_xy_avg		= 0.;
    dca_xy_sd	= 0.;
    Nprim 			= 0;
    run   			= 0;
    event_id    	= 0;
    ref2  			= 0;
    ref3        	= 0;
    btof 	    	= 0;
    qx          	= 0.;
    qy          	= 0.;
}

nsmEvent::nsmEvent(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z, Float_t mdca_xy_avg, Float_t mdca_xy_sd,
                   UInt_t mNprim, UInt_t mrun, UInt_t mevent_id, UInt_t mref2, UInt_t mref3, UInt_t mbtof,
				   Float_t mqx, Float_t mqy)
{
    vtx_x = mvtx_x;
    vtx_y = mvtx_y;
    vtx_z = mvtx_z;
    
    dca_xy_avg = mdca_xy_avg;
    dca_xy_sd = mdca_xy_sd;

    Nprim = mNprim;
    run  = mrun;
    event_id = mevent_id;

    ref2 = mref2;
    ref3 = mref3;
    
    btof = mbtof;

    qx = mqx;
    qy = mqy;
}

void nsmEvent::SetEventData(Float_t mvtx_x, Float_t mvtx_y, Float_t mvtx_z, Float_t mdca_xy_avg, Float_t mdca_xy_sd,
                            UInt_t mNprim, UInt_t mrun, UInt_t mevent_id, UInt_t mref2, UInt_t mref3, UInt_t mbtof,
							Float_t mqx, Float_t mqy
                            )
{
    vtx_x = mvtx_x;
    vtx_y = mvtx_y;
    vtx_z = mvtx_z;
    
    dca_xy_avg = mdca_xy_avg;
    dca_xy_sd = mdca_xy_sd;

    Nprim = mNprim;
    run  = mrun;
    event_id = mevent_id;

    ref2 = mref2;
    ref3 = mref3;
    
    btof = mbtof;

    qx = mqx;
    qy = mqy;
}
