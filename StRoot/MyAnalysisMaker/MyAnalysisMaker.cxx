//Dylan Edited 06/25/19

#include "MyAnalysisMaker.h"
#include "BES1_QA_Parameters.h"
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
    run_num        =      -999;  //just a no.
    
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
    
    // Check if trigger is good
    vector<int> good_triggers = triggers[energy];
    bool good_trig = false;
    for(int trig_index = 0; trig_index < (int)good_triggers.size(); trig_index++) {
    	if(muEvent->triggerIdCollection().nominal().isTrigger(good_triggers[trig_index])) {
    		good_trig = true;
    		break;
    	}
    }
    if(!good_trig) { return kTRUE; }

    event_cut_hist->Fill("Good Trigger", 1);

    // Check if run number is good
    run_num = muEvent->runId();
//    vector<int> bad_runs_energy = bad_runs[energy];
//    for(int bad_run_index = 0; bad_run_index < (int)bad_runs_energy.size(); bad_run_index++) {
//    	if(run_num == bad_run_energy[bad_run_index]) { return kTRUE; }
//    }

	event_cut_hist->Fill("Good Run", 1);

	// Check if vertex is good
	double vx = muEvent->primaryVertexPosition().x();
	double vy = muEvent->primaryVertexPosition().y();
	double vz = muEvent->primaryVertexPosition().z();
    
    if(energy == 7) {
    	if(fabs(vz)>50.0) {
			return kTRUE;
		}
    } else if(fabs(vz)>30.0) {
		return kTRUE; // Vertex within 30cm of detector center along beam pipe.
	}

    event_cut_hist->Fill("Good Vz", 1);

    if(energy == 14) {
		if(sqrt(pow(vx,2.)+pow((vy+0.89),2.))>1.)
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
	float dca_xy_sd = 0.;
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
		if(fabs(track->dcaD()) <= 4) {
			dca_xy_avg += track->dcaD();  // Check
			dca_xy_sd += pow(track->dcaD(), 2);  // Calculate second raw moment first
			dca_xy_count++;
		}

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
	if(dca_xy_count > 0) { dca_xy_avg /= dca_xy_count; dca_xy_sd = dca_xy_sd / dca_xy_count - pow(dca_xy_avg, 2); }
	else { dca_xy_avg = -899; dca_xy_sd = -899; }

    levent->SetEventData(muEvent->primaryVertexPosition().x(), muEvent->primaryVertexPosition().y(), muEvent->primaryVertexPosition().z(), dca_xy_avg, dca_xy_sd, muEvent->refMult(), run_num, muEvent->eventId(), ref2, ref3, muEvent->btofTrayMultiplicity(), Qx, Qy);
    
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
