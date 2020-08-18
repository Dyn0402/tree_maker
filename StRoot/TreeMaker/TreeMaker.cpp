/*
 * TreeMaker.cpp
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */


#include "TreeMaker.h"

ClassImp(TreeMaker)


TreeMaker::TreeMaker(StMuDstMaker *maker) : StMaker("TreeMaker") {
	muDst_maker = maker;
	muDst = NULL;

	picoDst_maker = NULL;
	picoDst = NULL;

	out_file_name = "";
	out_file = NULL;
	tree = NULL;

	event_cut_hist = NULL;
	track_cut_hist = NULL;
	de_dx_pq_hist = NULL;
	beta_pq_hist = NULL;

	// Temp QA plots
//	flag_diff_hist = NULL;
//	nHitsFit_diff_hist = NULL;
//	nHitsPoss_diff_hist = NULL;
//	dca_diff_hist = NULL;
//	dca_prim_glob_hist = NULL;
//	nHitsFit_diff_post_hist = NULL;
//	nHitsPoss_diff_post_hist = NULL;
//	dca_diff_post_hist = NULL;
//	dca_prim_glob_post_hist = NULL;

	events_read = 0;
	events_processed = 0;
	energy = 0;
	bes_phase = 1;

	pars.set_energy_bes(energy, bes_phase);
}

TreeMaker::TreeMaker(StMuDstMaker *maker, string name, int energy_in, int bes_phase) : StMaker("TreeMaker") {
	muDst_maker = maker;
	muDst = NULL;

	picoDst_maker = NULL;
	picoDst = NULL;

	out_file_name = name;
	out_file = NULL;
	tree = NULL;

	event_cut_hist = NULL;
	track_cut_hist = NULL;
	de_dx_pq_hist = NULL;
	beta_pq_hist = NULL;


	// Temp QA plots
//	flag_diff_hist = NULL;
//	nHitsFit_diff_hist = NULL;
//	nHitsPoss_diff_hist = NULL;
//	dca_diff_hist = NULL;
//	dca_prim_glob_hist = NULL;
//	nHitsFit_diff_post_hist = NULL;
//	nHitsPoss_diff_post_hist = NULL;
//	dca_diff_post_hist = NULL;
//	dca_prim_glob_post_hist = NULL;

	events_read = 0;
	events_processed = 0;
	energy = energy_in;
	this->bes_phase = bes_phase;

	pars.set_energy_bes(energy, bes_phase);
}


// Structors
TreeMaker::TreeMaker(StPicoDstMaker *maker) : StMaker("TreeMaker") {
	picoDst_maker = maker;
	picoDst = NULL;

	muDst_maker = NULL;
	muDst = NULL;

	out_file_name = "";
	out_file = NULL;
	tree = NULL;

	event_cut_hist = NULL;
	track_cut_hist = NULL;
	de_dx_pq_hist = NULL;
	beta_pq_hist = NULL;

	// Temp QA plots
//	flag_diff_hist = NULL;
//	nHitsFit_diff_hist = NULL;
//	nHitsPoss_diff_hist = NULL;
//	dca_diff_hist = NULL;
//	dca_prim_glob_hist = NULL;
//	nHitsFit_diff_post_hist = NULL;
//	nHitsPoss_diff_post_hist = NULL;
//	dca_diff_post_hist = NULL;
//	dca_prim_glob_post_hist = NULL;

	events_read = 0;
	events_processed = 0;
	energy = 0;
	bes_phase = 1;

	pars.set_energy_bes(energy, bes_phase);
}

TreeMaker::TreeMaker(StPicoDstMaker *maker, string name, int energy_in, int bes_phase) : StMaker("TreeMaker") {
	picoDst_maker = maker;
	picoDst = NULL;

	muDst_maker = NULL;
	muDst = NULL;

	out_file_name = name;
	out_file = NULL;
	tree = NULL;

	event_cut_hist = NULL;
	track_cut_hist = NULL;
	de_dx_pq_hist = NULL;
	beta_pq_hist = NULL;


	// Temp QA plots
//	flag_diff_hist = NULL;
//	nHitsFit_diff_hist = NULL;
//	nHitsPoss_diff_hist = NULL;
//	dca_diff_hist = NULL;
//	dca_prim_glob_hist = NULL;
//	nHitsFit_diff_post_hist = NULL;
//	nHitsPoss_diff_post_hist = NULL;
//	dca_diff_post_hist = NULL;
//	dca_prim_glob_post_hist = NULL;

	events_read = 0;
	events_processed = 0;
	energy = energy_in;
	this->bes_phase = bes_phase;

	pars.set_energy_bes(energy, bes_phase);
}

TreeMaker::~TreeMaker() {
	// Nothing
}



// Setters

void TreeMaker::set_energy(int energy_in) {
	energy = energy_in;
}

void TreeMaker::set_out_file_name(string name) {
	out_file_name = name;
}



// St Doers

Int_t TreeMaker::Init() {
	out_file = new TFile(out_file_name.data(),"RECREATE") ;
	tree = new TTree("tree", "tree");

	tree->Branch("run_num", &event.run_num, "run_num/I");
	tree->Branch("event_id", &event.event_id, "event_id/I");
	tree->Branch("refmult", &event.refmult, "refmult/S");
	tree->Branch("refmult2", &event.refmult2, "refmult2/S");
	tree->Branch("refmult3", &event.refmult3, "refmult3/S");
	tree->Branch("btof", &event.btof, "btof/S");
	tree->Branch("vx", &event.vx, "vx/F");
	tree->Branch("vy", &event.vy, "vy/F");
	tree->Branch("vz", &event.vz, "vz/F");
	tree->Branch("qx", &event.qx, "qx/F");
	tree->Branch("qy", &event.qy, "qy/F");
	tree->Branch("dca_xy_avg", &event.dca_xy_avg, "dca_xy_avg/F");
	tree->Branch("dca_xy_err", &event.dca_xy_err, "dca_xy_err/F");

	tree->Branch("proton.pt", &protons.pt, pars.branch_buffer, pars.branch_split);
	tree->Branch("proton.phi", &protons.phi, pars.branch_buffer, pars.branch_split);
	tree->Branch("proton.eta", &protons.eta, pars.branch_buffer, pars.branch_split);
	tree->Branch("proton.dca", &protons.dca, pars.branch_buffer, pars.branch_split);
	tree->Branch("proton.nsigma", &protons.nsigma, pars.branch_buffer, pars.branch_split);
	tree->Branch("proton.beta", &protons.beta, pars.branch_buffer, pars.branch_split);
	tree->Branch("proton.charge", &protons.charge, pars.branch_buffer, pars.branch_split);

	tree->Branch("pion.pt", &pions.pt, pars.branch_buffer, pars.branch_split);
	tree->Branch("pion.phi", &pions.phi, pars.branch_buffer, pars.branch_split);
	tree->Branch("pion.eta", &pions.eta, pars.branch_buffer, pars.branch_split);
	tree->Branch("pion.dca", &pions.dca, pars.branch_buffer, pars.branch_split);
	tree->Branch("pion.nsigma", &pions.nsigma, pars.branch_buffer, pars.branch_split);
	tree->Branch("pion.beta", &pions.beta, pars.branch_buffer, pars.branch_split);
	tree->Branch("pion.charge", &pions.charge, pars.branch_buffer, pars.branch_split);

	event_cut_hist = new TH1D("Event Cut Hist", "Event Cut Hist", 9, -0.5, 8.5);
	event_cut_hist->GetXaxis()->SetBinLabel(1, "Expected");
	event_cut_hist->GetXaxis()->SetBinLabel(2, "Events Read");
	event_cut_hist->GetXaxis()->SetBinLabel(3, "Is dstEvent");
	event_cut_hist->GetXaxis()->SetBinLabel(4, "Good Trigger");
	event_cut_hist->GetXaxis()->SetBinLabel(5, "Good Run");
	event_cut_hist->GetXaxis()->SetBinLabel(6, "Good Vz");
	event_cut_hist->GetXaxis()->SetBinLabel(7, "Good Vr");
	event_cut_hist->GetXaxis()->SetBinLabel(8, "Vertex Non-Zero");
	event_cut_hist->GetXaxis()->SetBinLabel(9, "Good VPD Vz");

	track_cut_hist = new TH1D("Track Cut Hist", "Track Cut Hist", 16, -0.5, 15.5);
	track_cut_hist->GetXaxis()->SetBinLabel(1, "Tracks Read");
	track_cut_hist->GetXaxis()->SetBinLabel(2, "Is Track");
	track_cut_hist->GetXaxis()->SetBinLabel(3, "Primary Track/Flags");
	track_cut_hist->GetXaxis()->SetBinLabel(4, "Charge");
	track_cut_hist->GetXaxis()->SetBinLabel(5, "nHitsRatio Min");
	track_cut_hist->GetXaxis()->SetBinLabel(6, "nHitsRatio Max");
	track_cut_hist->GetXaxis()->SetBinLabel(7, "eta");
	track_cut_hist->GetXaxis()->SetBinLabel(8, "nHitsFit");
	track_cut_hist->GetXaxis()->SetBinLabel(9, "nHitsDedx");
	track_cut_hist->GetXaxis()->SetBinLabel(10, "dca");
	track_cut_hist->GetXaxis()->SetBinLabel(11, "pt_low");
	track_cut_hist->GetXaxis()->SetBinLabel(12, "pt_high");
	track_cut_hist->GetXaxis()->SetBinLabel(13, "nsigma_proton");
	track_cut_hist->GetXaxis()->SetBinLabel(14, "m_proton");
	track_cut_hist->GetXaxis()->SetBinLabel(15, "nsigma_pion");
	track_cut_hist->GetXaxis()->SetBinLabel(16, "m_pion");

	de_dx_pq_hist = new TH2F("dedx_pq_pid", "Dedx PID", 1000, -3, 3, 1000, 0, 0.5e-4);
	beta_pq_hist = new TH2F("beta_pq_pid", "Beta PID", 1000, -3, 3, 1000, 0, 5);

	// Temp QA plots
//	flag_diff_hist = new TH1D("flag_diff_hist", "Flag Primary - Flag Global", 801, -400.5, 400.5);
//	nHitsFit_diff_hist = new TH1D("nHitsFit_diff_hist", "nHitsFit Primary - nHitsFit Global", 16, -10.5, 5.5);
//	nHitsPoss_diff_hist = new TH1D("nHitsPoss_diff_hist", "nHitsPoss Primary - nHitsPoss Global", 16, -10.5, 5.5);
//	dca_diff_hist = new TH1D("dca_diff_hist", "Dca Primary - Dca Global", 200, -10.0, 2.0);
//	dca_prim_glob_hist = new TH2D("dca_prim_glob_hist", "Dca Primary vs Dca Global", 1000, 0, 5, 2000, 0, 10);
//	nHitsFit_diff_post_hist = new TH1D("nHitsFit_diff_post_hist", "nHitsFit Primary - nHitsFit Global", 16, -10.5, 5.5);
//	nHitsPoss_diff_post_hist = new TH1D("nHitsPoss_diff_post_hist", "nHitsPoss Primary - nHitsPoss Global", 16, -10.5, 5.5);
//	dca_diff_post_hist = new TH1D("dca_diff_post_hist", "Dca Primary - Dca Global", 200, -10.0, 2.0);
//	dca_prim_glob_post_hist = new TH2D("dca_prim_glob_post_hist", "Dca Primary vs Dca Global", 1000, 0, 5, 2000, 0, 10);

	return kStOK;
}


Int_t TreeMaker::Make() {
	events_read++;
	event_cut_hist->Fill("Events Read", 1);

	event.clear(); protons.clear(); pions.clear();  // Clear event/particle objects before processing new event

	if(muDst_maker) {
		muDst = muDst_maker->muDst();
		StMuEvent* mu_event = muDst->event();  // Get muEvent from maker

		if(is_bad_event(mu_event)) { return kStOk; }  // Check if event is good, save event vars to event

		track_loop(mu_event);  // Loop over tracks in mu_event, save track vars to protons/pions
	}

	else if(picoDst_maker) {
		picoDst = picoDst_maker->picoDst();
		StPicoEvent* pico_event = picoDst->event();  // Get picoEvent from maker

		if(is_bad_event(pico_event)) { return kStOk; }  // Check if event is good, save event vars to event

		track_loop(pico_event);  // Loop over tracks in pico_event, save track vars to protons/pions
	}

	tree->Fill();  // Fill tree with event/protons/pions

	events_processed++;

	return kStOk;
}


Int_t TreeMaker::Finish() {
	cout << endl;
	cout << "Finishing and writing histograms to file... " << endl;
	cout << endl;

	out_file->Write();
	out_file->Close();

	cout <<"\n ======> Finished <======"<<endl;
	cout<<" Acutal #Events Read = " << events_read << endl ;
	cout<<" Acutal #Events Processed = " << events_processed << endl ;

	return kStOk;
}



// Doers
bool TreeMaker::is_bad_event(StMuEvent *mu_event) {
	if(!mu_event) { return true; }
	event_cut_hist->Fill("Is dstEvent", 1);

	// Check for good trigger
	vector<int> good_triggers = pars.triggers;
	bool good_trig = false;
	for(int trig_index = 0; trig_index < (int)pars.triggers.size(); trig_index++) {
		if(mu_event->triggerIdCollection().nominal().isTrigger(pars.triggers[trig_index])) {
			good_trig = true;
			break;
		}
	}
	if(!good_trig) { return true; }
	event_cut_hist->Fill("Good Trigger", 1);


    // Check if run number is good
    event.run_num = mu_event->runId();
    vector<int> bad_runs_energy = pars.bad_runs;
    int num_bad_runs = (int) bad_runs_energy.size();
    for(int bad_run_index = 0; bad_run_index < num_bad_runs; bad_run_index++) {
    	if(event.run_num == bad_runs_energy[bad_run_index]) {
    		return true;
    	}
    }
    event_cut_hist->Fill("Good Run", 1);

    // Get x,y,z components of primary vertex
	event.vx = mu_event->primaryVertexPosition().x();
	event.vy = mu_event->primaryVertexPosition().y();
	event.vz = mu_event->primaryVertexPosition().z();

	// Check vertex is within pars.vz_max cm of detector center along beam pipe
	if(fabs(event.vz) > pars.vz_max) { return true; }
	event_cut_hist->Fill("Good Vz", 1);

	// Check that vertex is within x cm radially (x-y plane) of detector axis
	if(sqrt(pow(event.vx, 2) + pow(event.vy + pars.vy_offset, 2)) > pars.vr_max) {
		return true;
	}
	event_cut_hist->Fill("Good Vr", 1);

	// On old tapes, no-vertex gets reported as VtxPosition=(0,0,0)
	if(fabs(event.vx) < pars.vertex_min &&
			fabs(event.vy) < pars.vertex_min &&
			fabs(event.vz) < pars.vertex_min) {
		return true;
	}
	event_cut_hist->Fill("Vertex Non-Zero", 1);

	// Filter out events with disagreement between vpd and vertex reconstruction.
	if(pars.vpd_vz_max_diff > 0) {  // -1 error code
		if(muDst->btofHeader()) {
			float vpd_vz = muDst->btofHeader()->vpdVz();
			if(fabs(vpd_vz - event.vz) > pars.vpd_vz_max_diff) {
				return true;
			}
		} else {
			return true;
		}
	}
	event_cut_hist->Fill("Good VPD Vz", 1);


	// Add other event variables to event
	event.event_id = mu_event->eventId();
	event.refmult = mu_event->refMult();
//	muDst->setVertexIndex(0);
	muDst_maker->printArrays();
	StMuPrimaryVertex *pv = muDst->primaryVertex();
	StMuPrimaryVertex *pv2 = (StMuPrimaryVertex*)muDst->primaryVertices()->UncheckedAt(0);
	cout << muDst->primaryVertices()->GetEntriesFast() << " " << pv << " " << pv2 << " " << pv->nBTOFMatch() << endl;
	if(pv) {
		event.btof = pv->nBTOFMatch();
	} else { event.btof = 0; }

	return false;  // If all above checks are passed, event is good
}


bool TreeMaker::is_bad_event(StPicoEvent *pico_event) {
	if(!pico_event) { return true; }
	event_cut_hist->Fill("Is dstEvent", 1);

	// Check for good trigger
	vector<int> good_triggers = pars.triggers;
	bool good_trig = false;
	for(int trig_index = 0; trig_index < (int)pars.triggers.size(); trig_index++) {
		if(pico_event->isTrigger(pars.triggers[trig_index])) {
			good_trig = true;
			break;
		}
	}
	if(!good_trig) { return true; }
	event_cut_hist->Fill("Good Trigger", 1);


    // Check if run number is good
    event.run_num = pico_event->runId();
    vector<int> bad_runs_energy = pars.bad_runs;
    int num_bad_runs = (int) bad_runs_energy.size();
    for(int bad_run_index = 0; bad_run_index < num_bad_runs; bad_run_index++) {
    	if(event.run_num == bad_runs_energy[bad_run_index]) {
    		return true;
    	}
    }
    event_cut_hist->Fill("Good Run", 1);

    // Get x,y,z components of primary vertex
	event.vx = pico_event->primaryVertex().X();
	event.vy = pico_event->primaryVertex().Y();
	event.vz = pico_event->primaryVertex().Z();

	// Check vertex is within pars.vz_max cm of detector center along beam pipe
	if(fabs(event.vz) > pars.vz_max) { return true; }
	event_cut_hist->Fill("Good Vz", 1);

	// Check that vertex is within x cm radially (x-y plane) of detector axis
	if(sqrt(pow(event.vx, 2) + pow(event.vy + pars.vy_offset, 2)) > pars.vr_max) {
		return true;
	}
	event_cut_hist->Fill("Good Vr", 1);

	// On old tapes, no-vertex gets reported as VtxPosition=(0,0,0)
	if(fabs(event.vx) < pars.vertex_min &&
			fabs(event.vy) < pars.vertex_min &&
			fabs(event.vz) < pars.vertex_min) {
		return true;
	}
	event_cut_hist->Fill("Vertex Non-Zero", 1);

	// Filter out events with disagreement between vpd and vertex reconstruction.
	if(pars.vpd_vz_max_diff > 0) {  // -1 ignore code
		float vpd_vz = pico_event->vzVpd();
		if(fabs(vpd_vz - event.vz) > pars.vpd_vz_max_diff) {
			return true;
		}
	}
	event_cut_hist->Fill("Good VPD Vz", 1);


	// Add other event variables to event
	event.event_id = pico_event->eventId();
	event.refmult = pico_event->refMult();
	event.btof = pico_event->nBTOFMatch();
	event.refmult2 = pico_event->refMult2();
	event.refmult3 = pico_event->refMult3();


	return false;  // If all above checks are passed, event is good
}


void TreeMaker::track_loop(StMuEvent *mu_event) {
	int num_primary = muDst->primaryTracks()->GetEntries();
	StMuTrack *track, *track_glob;

	int index_2g, nHitsFit, btofMatch, tofmatched = 0, tofmatchedbeta = 0, dca_xy_count = 0;
	float dca, dca_prim, eta, pt, nsigmapr, nsigmapi, phi, dca_xy_avg = 0, dca_xy_err = 0.;
	double ratio; // Important that this is double, 13/25 = 0.52 = cut!!!
	double beta, p, m;
	short charge;

	for(int track_index = 0; track_index < num_primary; track_index++) {
		track_cut_hist->Fill("Tracks Read", 1);
		track = (StMuTrack*) muDst->primaryTracks(track_index);

		// Initial track cuts

		if(!track) continue;  // Check that track not NULL
		track_cut_hist->Fill("Is Track", 1);

		if(track->vertexIndex() != 0) continue;  // Check that vertex index is zero

		index_2g = track->index2Global();
		if(index_2g < 0) continue;  // Check that global index non negative

		track_glob = (StMuTrack*) muDst->globalTracks(index_2g);

		if(track->flag() < 0) continue;  // Check primary track flag, still unsure what it is

		if(track_glob->flag() < 0) continue;  // Check global track flag, still unsure what it is

		track_cut_hist->Fill("Primary Track/Flags", 1);

		charge = track->charge();
		if(fabs(charge) != 1) continue;  // Eliminates neutral/exotic particles
		track_cut_hist->Fill("Charge", 1);

		// Temp QA plots
//		flag_diff_hist->Fill(fabs(track->flag() - track_glob->flag()));
//		nHitsFit_diff_hist->Fill(track->nHitsFit() - track_glob->nHitsFit());
//		nHitsPoss_diff_hist->Fill(track->nHitsPoss() - track_glob->nHitsPoss());
//		dca_diff_hist->Fill(track->dca().mag() - track->dcaGlobal().mag());
//		dca_prim_glob_hist->Fill(track->dca().mag(), track->dcaGlobal().mag());


		// Get main track variables

		p = track->p().mag();
		pt = track->pt();
		eta = track->eta();
		phi = track->phi();  if(phi < 0) { phi += 2*M_PI; }
		dca = track->dcaGlobal().mag();
		dca_prim = track->dca().mag();
		nsigmapr = track->nSigmaProton();

		nHitsFit = track_glob->nHitsFit();

		btofMatch = track->btofPidTraits().matchFlag();
		beta = track->btofPidTraits().beta();
		m = (beta > 1.e-5) ? p*p*(1./beta/beta - 1.) : -999;


		// Event track counters

		if(btofMatch > 0 && fabs(eta) < 0.5 && dca_prim < 3.0 && nHitsFit > 10) {
			tofmatched++;
			if(beta > 0.1) tofmatchedbeta++;
		}

//		if(beta > 0.1 && fabs(eta) < 1. && dca_prim < 3. && nHitsFit > 10) { } //betamatch

		if(fabs(eta) > 0.5 && fabs(eta) < 1. && dca_prim <= 3. && nHitsFit >= 10 && p >= 1.e-10) event.refmult2++;
		if(fabs(eta) < 1. && nHitsFit >= 10 && dca_prim <= 3. && nsigmapr < -3. && m < 0.4 && p >= 1.e-10) event.refmult3++;

		// Cut on ratio of nHitsFit to nHitsPossible
		ratio = (double) nHitsFit / (double) track_glob->nHitsPoss();
		if(ratio < 0.52) continue;
		track_cut_hist->Fill("nHitsRatio Min", 1);
		if(ratio > 1.05) continue;
		track_cut_hist->Fill("nHitsRatio Max", 1);

		// Event Plane Q vector
		if(nHitsFit > 15 && dca < 2.0 && fabs(eta) < 1.0 && pt > 0.2 && pt < 2.) {
			event.qx += cos(2*phi); event.qy += sin(2*phi);
		}

		// Fill PID plots
		de_dx_pq_hist->Fill(charge*p, track->dEdx());
		if(beta > 1.e-5) {
			beta_pq_hist->Fill(charge*p, 1 / beta);
		}

		// Calculate dca_xy variables
		if(track->dcaD() < 4 && track->dcaD() >= -4) {
			dca_xy_avg += track->dcaD();
			dca_xy_err += pow(track->dcaD(), 2);  // Calculate second raw moment first
			dca_xy_count++;
		}

		if(fabs(eta) > 1.0) continue;
		track_cut_hist->Fill("eta", 1);

		if(nHitsFit <= 20) continue;
		track_cut_hist->Fill("nHitsFit", 1);
		if(track->nHitsDedx() <= 5) continue;
		track_cut_hist->Fill("nHitsDedx", 1);

		if(dca < 0 || dca > 1.2) continue;
		track_cut_hist->Fill("dca", 1);

		if(pt < 0.3) continue;
		track_cut_hist->Fill("pt_low", 1);
		if(pt > 2.2) continue;
		track_cut_hist->Fill("pt_high", 1);

		// Temp QA plots
//		nHitsFit_diff_post_hist->Fill(track->nHitsFit() - track_glob->nHitsFit());
//		nHitsPoss_diff_post_hist->Fill(track->nHitsPoss() - track_glob->nHitsPoss());
//		dca_diff_post_hist->Fill(track->dca().mag() - track->dcaGlobal().mag());
//		dca_prim_glob_post_hist->Fill(track->dca().mag(), track->dcaGlobal().mag());

		nsigmapi = track->nSigmaPion();

		if(energy == 27) {
			if(fabs(nsigmapr) <= 1.2) {
				track_cut_hist->Fill("nsigma_proton", 1);
				if( (m > 0.6 && m < 1.2) || m == -999) {
					track_cut_hist->Fill("m_proton", 1);
					protons.add_event(pt, phi, eta, dca, nsigmapr, beta, charge);
				}
			} if(fabs(nsigmapi) <= 1.0) {
				track_cut_hist->Fill("nsigma_pion", 1);
				if( (m > -0.15 && m < 0.15) || m == -999) {
					track_cut_hist->Fill("m_pion", 1);
					pions.add_event(pt, phi, eta, dca, nsigmapi, beta, charge);
				}
			}
		} else {
			if(fabs(nsigmapr) <= 2.2) {
				track_cut_hist->Fill("nsigma_proton", 1);
				if( (m > 0.6 && m < 1.2) || m == -999) {
					track_cut_hist->Fill("m_proton", 1);
					protons.add_event(pt, phi, eta, dca, nsigmapr, beta, charge);
				}
			} if(fabs(nsigmapi) <= 2.0) {
				track_cut_hist->Fill("nsigma_pion", 1);
				if( (m > -0.15 && m < 0.15) || m == -999) {
					track_cut_hist->Fill("m_pion", 1);
					pions.add_event(pt, phi, eta, dca, nsigmapi, beta, charge);
				}
			}
		}

	}

	// Calculate and set dca_xy variables in event
	if(dca_xy_count > 0) { event.dca_xy_avg = dca_xy_avg / dca_xy_count; event.dca_xy_err = pow((dca_xy_err / dca_xy_count - pow(event.dca_xy_avg, 2)) / dca_xy_count, 0.5); }
	else { event.dca_xy_avg = -899; event.dca_xy_err = -899; }

}


void TreeMaker::track_loop(StPicoEvent *pico_event) {
	int num_tracks = picoDst->numberOfTracks();
	StPicoTrack *track;

	int nHitsFit, dca_xy_count = 0;
	float dca, eta, pt, nsigmapr, nsigmapi, phi, dcas, dca_xy_avg = 0, dca_xy_err = 0.;
	double ratio; // Important that this is double, 13/25 = 0.52 = cut!!!
	double beta, p, m;
	short charge;

	for(int track_index = 0; track_index < num_tracks; track_index++) {
		track_cut_hist->Fill("Tracks Read", 1);
		track = (StPicoTrack*) picoDst->track(track_index);

		// Initial track cuts

		if(!track) continue;  // Check that track not NULL
		track_cut_hist->Fill("Is Track", 1);

		if(!track->isPrimary()) continue;  // Check track is primary track
		track_cut_hist->Fill("Primary Track/Flags", 1);

		charge = track->charge();
		if(fabs(charge) != 1) continue;  // Eliminates neutral/exotic particles
		track_cut_hist->Fill("Charge", 1);

		// Temp QA plots
//		flag_diff_hist->Fill(fabs(track->flag() - track_glob->flag()));
//		nHitsFit_diff_hist->Fill(track->nHitsFit() - track_glob->nHitsFit());
//		nHitsPoss_diff_hist->Fill(track->nHitsPoss() - track_glob->nHitsPoss());
//		dca_diff_hist->Fill(track->dca().mag() - track->dcaGlobal().mag());
//		dca_prim_glob_hist->Fill(track->dca().mag(), track->dcaGlobal().mag());


		// Get main track variables

		p = track->pMom().Mag();
		pt = track->pMom().Perp();
		eta = track->pMom().PseudoRapidity();
		phi = track->pMom().Phi();  if(phi < 0) { phi += 2*M_PI; }
		dca = track->gDCA(event.vx, event.vy, event.vz);
		dcas = 0; // track->gDCAs(pico_event->primaryVertex());  Still in dev only version, hopefully will go to pro before I need it
//		dca_prim = track->dca().mag();
		nsigmapr = track->nSigmaProton();

		nHitsFit = track->nHitsFit();

		int btof_pid_traits_index = track->bTofPidTraitsIndex();
		if(btof_pid_traits_index >= 0) {
			StPicoBTofPidTraits *btof_pid_traits = picoDst->btofPidTraits(btof_pid_traits_index);
			beta = btof_pid_traits->btofBeta();
			m = (beta > 1.e-5) ? p*p*(1./beta/beta - 1.) : -999;
		}

		// Event track counters

//		if(btofMatch > 0 && fabs(eta) < 0.5 && dca_prim < 3.0 && nHitsFit > 10) {
//			tofmatched++;
//			if(beta > 0.1) tofmatchedbeta++;
//		}

//		if(beta > 0.1 && fabs(eta) < 1. && dca_prim < 3. && nHitsFit > 10) { } //betamatch

//		if(fabs(eta) > 0.5 && fabs(eta) < 1. && dca_prim < 3. && nHitsFit > 10) refmult2++;
//		if(fabs(eta) < 1. && nHitsFit > 10 && dca_prim < 3. && nsigmapr < -3. && m < 0.4) refmult3++;

		// Cut on ratio of nHitsFit to nHitsPossible
		ratio = (double) nHitsFit / (double) track->nHitsMax();
		if(ratio < 0.52) continue;
		track_cut_hist->Fill("nHitsRatio Min", 1);
		if(ratio > 1.05) continue;
		track_cut_hist->Fill("nHitsRatio Max", 1);

		// Event Plane Q vector
		if(nHitsFit > 15 && dca < 2.0 && fabs(eta) < 1.0 && pt > 0.2 && pt < 2.) {
			event.qx += cos(2*phi); event.qy += sin(2*phi);
		}

		// Fill PID plots
		de_dx_pq_hist->Fill(charge*p, track->dEdx());
		if(beta > 1.e-5) {
			beta_pq_hist->Fill(charge*p, 1 / beta);
		}

		// Calculate dca_xy variables
		if(dcas < 4 && dcas >= -4) {
			dca_xy_avg += dcas;
			dca_xy_err += pow(dcas, 2);  // Calculate second raw moment first
			dca_xy_count++;
		}

		if(fabs(eta) > 1.0) continue;
		track_cut_hist->Fill("eta", 1);

		if(nHitsFit <= 20) continue;
		track_cut_hist->Fill("nHitsFit", 1);
		if(track->nHitsDedx() <= 5) continue;
		track_cut_hist->Fill("nHitsDedx", 1);

		if(dca < 0 || dca > 3.0) continue;
		track_cut_hist->Fill("dca", 1);

		if(pt < 0.3) continue;
		track_cut_hist->Fill("pt_low", 1);
		if(pt > 2.2) continue;
		track_cut_hist->Fill("pt_high", 1);

		// Temp QA plots
//		nHitsFit_diff_post_hist->Fill(track->nHitsFit() - track_glob->nHitsFit());
//		nHitsPoss_diff_post_hist->Fill(track->nHitsPoss() - track_glob->nHitsPoss());
//		dca_diff_post_hist->Fill(track->dca().mag() - track->dcaGlobal().mag());
//		dca_prim_glob_post_hist->Fill(track->dca().mag(), track->dcaGlobal().mag());

		nsigmapi = track->nSigmaPion();

		if(energy == 27) {
			if(fabs(nsigmapr) <= 1.2) {
				track_cut_hist->Fill("nsigma_proton", 1);
				if( (m > 0.6 && m < 1.2) || m == -999) {
					track_cut_hist->Fill("m_proton", 1);
					protons.add_event(pt, phi, eta, dca, nsigmapr, beta, charge);
				}
			} if(fabs(nsigmapi) <= 1.0) {
				track_cut_hist->Fill("nsigma_pion", 1);
				if( (m > -0.15 && m < 0.15) || m == -999) {
					track_cut_hist->Fill("m_pion", 1);
					pions.add_event(pt, phi, eta, dca, nsigmapi, beta, charge);
				}
			}
		} else {
			if(fabs(nsigmapr) <= 2.2) {
				track_cut_hist->Fill("nsigma_proton", 1);
				if( (m > 0.6 && m < 1.2) || m == -999) {
					track_cut_hist->Fill("m_proton", 1);
					protons.add_event(pt, phi, eta, dca, nsigmapr, beta, charge);
				}
			} if(fabs(nsigmapi) <= 2.0) {
				track_cut_hist->Fill("nsigma_pion", 1);
				if( (m > -0.15 && m < 0.15) || m == -999) {
					track_cut_hist->Fill("m_pion", 1);
					pions.add_event(pt, phi, eta, dca, nsigmapi, beta, charge);
				}
			}
		}

	}

	// Calculate and set dca_xy variables in event
	if(dca_xy_count > 0) { event.dca_xy_avg = dca_xy_avg / dca_xy_count; event.dca_xy_err = pow((dca_xy_err / dca_xy_count - pow(event.dca_xy_avg, 2)) / dca_xy_count, 0.5); }
	else { event.dca_xy_avg = -899; event.dca_xy_err = -899; }

}
