//Dylan Edited 06/25/19

#include "MyAnalysisMaker.h"
#include <iostream>

#include "StMuDSTMaker/COMMON/StMuDstMaker.h"
#include "StMuDSTMaker/COMMON/StMuTrack.h"
#include "StMuDSTMaker/COMMON/StMuEvent.h"
#include "StBTofHeader.h"
#include "StMessMgr.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TMath.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TObjArray.h"
#include "TComplex.h"
#include "nsmEvent.h"
#include "nsmTrack.h"
#include "TClonesArray.h"
#include "TTree.h"



using namespace std;

ClassImp(MyAnalysisMaker)                       // Macro for CINT compatibility
MyAnalysisMaker::MyAnalysisMaker(StMuDstMaker* maker) : StMaker("MyAnalysisMaker")
{                                               // Initialize data members here.
    mMuDstMaker      = maker ;                    // Pass MuDst pointer to DstAnlysisMaker Class member functions
    histogram_output = NULL  ;                    // Zero the Pointer to histogram output file
    mEventsRead = 0     ;                    // Zero the Number of Events read by the maker
    mEventsProcessed = 0     ;                    // Zero the Number of Events processed by the maker
    OutputFileName = "" ;                         // Output File Name( will be set inside the "readMuDst".C )
    energy = 0 ;
}


MyAnalysisMaker::~MyAnalysisMaker()
{/* */}


Int_t MyAnalysisMaker::Init()
{
    //----------------------------------------
    runnumber        =      -999;  //just a no.
    
    //-----------------------------------------------------------------------------------------
    histogram_output = new TFile(OutputFileName,"RECREATE") ;
    
    tree = new TTree("tree","tree");//

    protonArr = new TClonesArray("nsmTrack", 1000);
    pionArr = new TClonesArray("nsmTrack", 1000);
    tree->Branch("Proton", &protonArr, 256000, 99);
    tree->Branch("Pion", &pionArr, 256000, 99);

    levent = new nsmEvent();//                                                                                                                          
    tree->Branch("Event", "nsmEvent", &levent, 256000, 99);
    
    event_cut_hist = new TH1I("Event Cut Hist", "Event Cut Hist", 8, -0.5, 7.5);
	event_cut_hist->GetXaxis()->SetBinLabel(1, "Original");
	event_cut_hist->GetXaxis()->SetBinLabel(2, "Is muEvent");
	event_cut_hist->GetXaxis()->SetBinLabel(3, "Good Trigger");
	event_cut_hist->GetXaxis()->SetBinLabel(4, "Good Run");
	event_cut_hist->GetXaxis()->SetBinLabel(5, "Good Vz");
	event_cut_hist->GetXaxis()->SetBinLabel(6, "Good Vr");
	event_cut_hist->GetXaxis()->SetBinLabel(7, "Vertex Non-Zero");
	event_cut_hist->GetXaxis()->SetBinLabel(8, "Good VPD Vz");

	track_cut_hist = new TH1I("Track Cut Hist", "Track Cut Hist", 13, -0.5, 12.5);
	track_cut_hist->GetXaxis()->SetBinLabel(1, "Original");
	track_cut_hist->GetXaxis()->SetBinLabel(2, "Charge");
	track_cut_hist->GetXaxis()->SetBinLabel(3, "p_low");
	track_cut_hist->GetXaxis()->SetBinLabel(4, "ratio_low");
	track_cut_hist->GetXaxis()->SetBinLabel(5, "ratio_high");
	track_cut_hist->GetXaxis()->SetBinLabel(6, "eta");
	track_cut_hist->GetXaxis()->SetBinLabel(7, "nHitsFit");
	track_cut_hist->GetXaxis()->SetBinLabel(8, "nHitsDedx");
	track_cut_hist->GetXaxis()->SetBinLabel(9, "dca");
	track_cut_hist->GetXaxis()->SetBinLabel(10, "pt_low");
	track_cut_hist->GetXaxis()->SetBinLabel(11, "pt_high");
	track_cut_hist->GetXaxis()->SetBinLabel(12, "nsigma");
	track_cut_hist->GetXaxis()->SetBinLabel(13, "m");

//	de_dx_pq_hist = new TH2F("dedx_pq_pid", "Dedx PID", 1000, -3, 3, 1000, 0, 0.5e-4);
//	beta_pq_hist = new TH2F("beta_pq_pid", "Beta PID", 1000, -3, 3, 1000, 0, 5);

	VertexZPos = -100.0;
	VpdVzPos   = -100.0;

	return kStOK ;
}

void MyAnalysisMaker::SetEnergy(int energy_in) {
	energy = energy_in;
}


Bool_t MyAnalysisMaker::IsBadEvent(StMuEvent *muEvent)
{
    if(!muEvent) {
    	return kTRUE;
    }
    
    event_cut_hist->Fill("Is muEvent", 1);

    // Only accept events with good trigger.

    //Vpd-Zdc-mb & Vpd-Zdc-mb-protected  AuAu2011
    if(energy == 200) {
		if(!(muEvent->triggerIdCollection().nominal().isTrigger(350001) || muEvent->triggerIdCollection().nominal().isTrigger(350011) || muEvent->triggerIdCollection().nominal().isTrigger(350003) || muEvent->triggerIdCollection().nominal().isTrigger(350013) || muEvent->triggerIdCollection().nominal().isTrigger(350023) || muEvent->triggerIdCollection().nominal().isTrigger(350033) || muEvent->triggerIdCollection().nominal().isTrigger(350043)))
			return kTRUE;
    }
    
    //mb AuAu62
    else if(energy == 62) {
		if(!(muEvent->triggerIdCollection().nominal().isTrigger(270001) || muEvent->triggerIdCollection().nominal().isTrigger(270011) || muEvent->triggerIdCollection().nominal().isTrigger(270021)))
			return kTRUE;
    }

    //mb AuAu39
    else if(energy == 39) {
		if(!(muEvent->triggerIdCollection().nominal().isTrigger(280001)))
		return kTRUE;
    }
   
    //mb AuAu27
    else if(energy == 27) {
		if(!(muEvent->triggerIdCollection().nominal().isTrigger(360001)))
		return kTRUE;
    }
   
    //mb AuAu19
    else if(energy == 19) {
		if(!(muEvent->triggerIdCollection().nominal().isTrigger(340001) || muEvent->triggerIdCollection().nominal().isTrigger(340011) || muEvent->triggerIdCollection().nominal().isTrigger(340021))){
			return kTRUE;
		}
    }

    //mb AuAu14
    else if(energy == 14) {
	   if(!(muEvent->triggerIdCollection().nominal().isTrigger(440001) || muEvent->triggerIdCollection().nominal().isTrigger(440004) || muEvent->triggerIdCollection().nominal().isTrigger(440005)))
	   return kTRUE;
    }
    
    //mb AuAu11
    else if(energy == 11) {
		 if(!(muEvent->triggerIdCollection().nominal().isTrigger(310004) || muEvent->triggerIdCollection().nominal().isTrigger(310014)))
		 return kTRUE;
    }
   
    //mb AuAu7
    else if(energy == 7) {
		if(!(muEvent->triggerIdCollection().nominal().isTrigger(290004) || muEvent->triggerIdCollection().nominal().isTrigger(290001)))
		return kTRUE;
    }
    
    event_cut_hist->Fill("Good Trigger", 1);


    runnumber = muEvent->runId();

	int BadRunList[1000] = {11114085, 11114089, 11114095, 11115005, 11115042, 11115064, 11115080, 11116027, 11116028, 11116060, 11116068, 11116070, 11117044, 11117075, 11117088, 11117090, 11117098, 11118007, 11118008, 11118015, 11119024, 11119062, 11120040, 11120042, 11121044, 11121054, 11122037, 11122049, 11122050, 11122073, 11123028, 11123051, 11123053, 11123054, 11123055, 11123056, 11123057, 11123058, 11123059, 11123060, 11123076, 11123095, 11123102, 11124041, 11124046, 11124047, 11124048, 11124050, 11124062, 11128056, 11129028, 11133006, 11136005, 11136012, 11136013, 11142117, 11147006, 11148001, 11148045, 11149015, 11149017, 11149018, 11149047, 11150025, 11150029, 11151050, 11151057, 11152016, 11152036, 11152078, 11153032, 11154026, 11156009, 11156036, 11156043, 11156044, 11156045, 11157050, 11158011, 11158012, 11158015, 11158016, 11158017, 11158019, 11158020, 11158021, 11158022, 11158024, 15046096, 15046106, 15046107, 15046108, 15046110, 15046111, 15047004, 15047015, 15047019, 15047021, 15047023, 15047027, 15047028, 15047029, 15047030, 15047039, 15047041, 15047044, 15047050, 15047056, 15047057, 15047061, 15047062, 15047064, 15047068, 15047069, 15047071, 15047074, 15047075, 15047082, 15047093, 15051131, 15051132, 15053052, 15053053, 15053054, 15053055, 12113091, 12114007, 12115014, 12115015, 12115016, 12115018, 12115019, 12115020, 12115022, 12115023, 12115026, 12115073, 12116012, 12116015, 12116016, 12120026, 12121017, 12121022, 12121034, 12122019, 12172050, 12172051, 12172055, 12173009, 12173030, 12173031, 12173032, 12173033, 12173034, 12174077, 12174085, 12174096, 12175062, 12175087, 12175113, 12175114, 12175115, 12176054, 12176067, 12176069, 12177092, 12177099, 12177101, 12177106, 12177107, 12177108, 12178003, 12178004, 12178005, 12178006, 12178051, 12178099, 12178123, 12179068, 11099124, 11100034, 11100045, 11101046, 11101104, 11102012, 11102098, 11103008, 11103009, 11103046, 11103058, 11103062, 11103065, 11105011, 11105018, 11105029, 11106026, 11106027, 11106028, 11106029, 11106030, 11106040, 11106041, 11107061, 11108101, 11109013, 11109088, 11109090, 11109092, 11109097, 11109098, 11109100, 11109101, 11109102, 11109103, 11109104, 11109105, 11110041, 11110073, 11084009, 11085015, 11087001, 11087002, 11087003, 11087004, 11087057, 11087058, 11087059, 11089026, 11089048, 11089055, 11092012, 11092031, 11092032, 11092033, 11092034, 11095022, 11095048, 11095050, 11095051, 11095061, 11095062, 11095063, 11095064, 11002145, 11003001, 11003002, 11003006, 11003012, 11003013, 11023046, 11023048, 11025034, 11025054, 11028004, 11028005, 11028006, 11028018, 11028019, 11028020, 11028021, 11028022, 11057035, 11064003, 11070051, 11070052, 11070053, 11070054, 11070055, 11070056, 11070057, 11071001, 11071002, 11071003, 11071004, 11071006, 11071007, 11071013, 11071014, 11071015, 11071016, 11071017, 11071018, 11071019, 11071020, 11072030, 11072031, 11072032, 11072044, 11072045, 11073001, 11073003, 0, 0};

	int iii = 0;
	while(iii != 0) //0 end of bad run list
	{
		if(runnumber==BadRunList[iii]) return kTRUE; // Filter out events in bad runs.
		++iii;
	}

	event_cut_hist->Fill("Good Run", 1);

	double vx = muEvent->primaryVertexPosition().x();
	double vy = muEvent->primaryVertexPosition().y();
	double vz = muEvent->primaryVertexPosition().z();
    
    if(energy == 7) {
    	if(fabs(vz)>50.0) { // change to 50 for 7 GeV
			return kTRUE;
		}
    } else if(fabs(vz)>30.0) { // change to 50 for 7 GeV
		return kTRUE; // Vertex within 30cm of detector center along beam pipe.
	}

    event_cut_hist->Fill("Good Vz", 1);

    if(energy == 14) {
		if(sqrt(pow(vx,2.)+pow((vy+0.89),2.))>1.) //for 14 GeV
			return kTRUE;
    } else if(sqrt(vx*vx+vy*vy)>2.0) {
    	return kTRUE; // Vertex within 2cm radially of detector center axis.
    }

    event_cut_hist->Fill("Good Vr", 1);

	if( (vx < 1.e-5 && vx > -1.e-5) &&
	   (vy < 1.e-5 && vy > -1.e-5) &&
	   (vz < 1.e-5 && vz > -1.e-5)  ) {
		return kTRUE; // Too close to zero?
	}

	event_cut_hist->Fill("Vertex Non-Zero", 1);

    return kFALSE;
}


Int_t MyAnalysisMaker::Make()
{
	++mEventsRead;
	event_cut_hist->Fill("Original", 1);
    StMuEvent* muEvent  =  mMuDstMaker->muDst()->event();

    if(IsBadEvent(muEvent))  {                                     //Nominal Event cuts and trigger cut
    	return           kStOK;
    }

    //----------------------------------------------------

    // Filter out events with disagreement between vpd and vertex reconstruction.
    if(energy >= 39) {
    	if(mMuDstMaker->muDst()->btofHeader()) {
    		VpdVzPos    =  mMuDstMaker->muDst()->btofHeader()->vpdVz();
    		VertexZPos  =  muEvent-> primaryVertexPosition().z();
    		if(fabs(VpdVzPos-VertexZPos) > 3) return kStOK; // for 39,62 GeV
    	} else { return kStOK; }
    }
    
    //---------------------------------------------------------
    
    event_cut_hist->Fill("Good VPD Vz", 1);

	int nHitsFit, nHitsDedx;
	float ratio, dca, eta, pt, nsigmapr, nsigmapi, phi, charge, Qx, Qy;
	double beta, p, m;

	float dca_xy_avg = 0.;
	int dca_xy_count = 0;

	int protonp = 0; int pionp = 0;
	int ref2 = 0, ref3 = 0;
	Qx = 0; Qy = 0;

	TObjArray* tracks = mMuDstMaker->muDst()->primaryTracks() ;    // Create a TObject array containing the primary tracks
	TObjArrayIter  GetTracks(tracks) ;                              // Create an iterator to step through the tracks
	StMuTrack*                 track ;                              // Pointer to a track

	while((track = (StMuTrack*)GetTracks.Next()))
	{
		track_cut_hist->Fill("Original", 1);
		// Track quality cuts----------------------
		charge = track->charge();
		if(fabs(charge)!=1) continue; // Eliminates neutral particles
		track_cut_hist->Fill("Charge", 1);

		p = track->p().mag();
		if (p < 0.15) continue;
		track_cut_hist->Fill("p_low", 1);

		nHitsFit =  track->nHitsFit();
		nHitsFit =  fabs(nHitsFit)+1;
		ratio    =  (float) nHitsFit / (float) track->nHitsPoss();
		if(ratio < 0.52) continue;
		track_cut_hist->Fill("ratio_low", 1);
		if(ratio > 1.05) continue;
		track_cut_hist->Fill("ratio_high", 1);

//		de_dx_pq_hist->Fill(charge*p, track->dEdx());

		dca_xy_avg += track->dcaD(); // Fix
		dca_xy_count++;

		eta = track->eta();

		cout << "pre-global track dcaD: " << track->dcaD() << endl;
		cout << "pre-global track dcaZ: " << track->dcaZ() << endl;

		dca = track->dcaGlobal().mag();
		cout << "dcaGlobal.mag: " << dca << endl;

		cout << "post-global track dcaD: " << track->dcaD() << endl;
		cout << "post-global track dcaZ: " << track->dcaZ() << endl << endl;

		pt = track->pt();
		phi = track->phi();
		nsigmapr = track->nSigmaProton();
		nHitsDedx = track->nHitsDedx();
		if(phi < 0) phi = phi + 2 * TMath::Pi();

		beta = -999;
		beta = track->btofPidTraits().beta();
		m = -999;
		if(beta > 1.e-5) {
			m = p*p*(1./(beta*beta) - 1.);
//			beta_pq_hist->Fill(charge*p, 1 / beta);
		}

		// ref2
		if(nHitsFit > 10 && dca < 3.0 && fabs(eta) > 0.5 && fabs(eta) < 1.0) ref2++;

		// ref3
		if(nHitsFit > 10 && dca < 3.0 && fabs(eta) < 1.0 && m < 0.4 && nsigmapr < -3.0) ref3++;

		// Q vector for event plane
		if(nHitsFit > 15 && dca < 2.0 && fabs(eta) < 1.0 && pt > 0.2 && pt < 2.) {
			Qx += cos(2*phi); Qy += sin(2*phi);
		}

		if(fabs(eta) > 1.0) continue;
		track_cut_hist->Fill("eta", 1);

		if(nHitsFit < 20) continue;
		track_cut_hist->Fill("nHitsFit", 1);
		if(nHitsDedx <= 5) continue;
		track_cut_hist->Fill("nHitsDedx", 1);

		if(dca < 0 || dca > 2.2) continue;
		track_cut_hist->Fill("dca", 1);

		if(pt < 0.3) continue;
		track_cut_hist->Fill("pt_low", 1);
		if(pt > 2.5) continue;
		track_cut_hist->Fill("pt_high", 1);

		nsigmapi = track->nSigmaPion();

		if(energy == 27) {
			if(fabs(nsigmapr) <= 1.2 && fabs(eta) <= 1.0) {
				track_cut_hist->Fill("nsigma", 1);
				if( (m > 0.75 && m < 1.05) || m == -999) {
					track_cut_hist->Fill("m", 1);
					new((*protonArr)[protonp++]) nsmTrack(pt,phi,eta,dca,nsigmapr,beta,charge);
				}
			} if(fabs(nsigmapi <= 1.0 && fabs(eta) <= 0.5)) {
				track_cut_hist->Fill("nsigma", 1);
				if( (m > -0.15 && m < 0.15) || m == -999) {
					track_cut_hist->Fill("m", 1);
					new((*pionArr)[pionp++]) nsmTrack(pt,phi,eta,dca,nsigmapi,beta,charge);
				}
			}
		} else {
			if(fabs(nsigmapr) <= 2.2) {
				track_cut_hist->Fill("nsigma", 1);
				if( (m > 0.75 && m < 1.05) || m == -999) {
					track_cut_hist->Fill("m", 1);
					new((*protonArr)[protonp++]) nsmTrack(pt,phi,eta,dca,nsigmapr,beta,charge);
				}
			} if(fabs(nsigmapi <= 2.0)) {
				track_cut_hist->Fill("nsigma", 1);
				if( (m > -0.15 && m < 0.15) || m == -999) {
					track_cut_hist->Fill("m", 1);
					new((*pionArr)[pionp++]) nsmTrack(pt,phi,eta,dca,nsigmapi,beta,charge);
				}
			}
		}

		// Cuts selecting relevant particles----------------------

    }//==================track loop ends=========================

	cout << "pre dca_xy_count: " << dca_xy_count << "  |  dca_xy_avg: " << dca_xy_avg << endl;
	if(dca_xy_count > 0) { dca_xy_avg /= dca_xy_count; }
	else { dca_xy_avg = -9999; }

    levent->SetEventData(muEvent->primaryVertexPosition().x(), muEvent->primaryVertexPosition().y(), muEvent->primaryVertexPosition().z(), dca_xy_avg, muEvent->refMult(), runnumber, muEvent->eventId(), ref2, ref3, muEvent->btofTrayMultiplicity(), Qx, Qy);
    
    //fill tree
    tree->Fill();
    protonArr->Clear();
    pionArr->Clear();
 
    mEventsProcessed++ ;
    return kStOK ;
    
}//--------------- Make (Event) loop ends--------------------------


Int_t MyAnalysisMaker::Finish()
{
    cout<<" Inside Finish and writing histograms..."<<endl;
    cout << endl;
    cout << endl;

    histogram_output -> Write();
    histogram_output ->Close();
    
    cout << "donzo" << endl;

    cout <<"\n ======> All done <======"<<endl;
    cout<<" Acutal #Events Read = " <<mEventsRead<<"\n###### Thank You ######\n"<< endl ;
    cout<<" Acutal #Events Processed = " <<mEventsProcessed<<"\n###### Thank You ######\n"<< endl ;
    
    return kStOk ;
}
