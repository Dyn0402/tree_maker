/*
 * TreeMaker.cpp
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */


#include "TreeMaker.h"

ClassImp(TreeMaker)


// Structors
TreeMaker::TreeMaker(StMuDstMaker *maker) : StMaker("TreeMaker") {
	mudst_maker = maker;

	out_file_name = "";
	out_file = NULL;
	tree = NULL;

	event_cut_hist = NULL;
	track_cut_hist = NULL;
	de_dx_pq_hist = NULL;
	beta_pq_hist = NULL;

	events_read = 0;
	events_processed = 0;
	energy = 0;
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

	event.clear();
	protons.clear();
	pions.clear();

	tree->Branch("run_num", &event.run_num, "run_num/I");
	tree->Branch("event_id", &event.event_id, "event_id/I");
	tree->Branch("refmult", &event.refmult, "refmult/I");
	tree->Branch("refmult2", &event.refmult2, "refmult2/I");
	tree->Branch("refmult3", &event.refmult3, "refmult3/I");
	tree->Branch("btof", &event.btof, "btof/I");
	tree->Branch("vx", &event.vx, "vx/F");
	tree->Branch("vy", &event.vy, "vy/F");
	tree->Branch("vz", &event.vz, "vz/F");
	tree->Branch("qx", &event.qx, "qx/F");
	tree->Branch("qy", &event.qy, "qy/F");
	tree->Branch("dca_xy_avg", &event.dca_xy_avg, "dca_xy_avg/F");
	tree->Branch("dca_xy_err", &event.dca_xy_err, "dca_xy_err/F");

	tree->Branch("proton.pt", &protons.pt, pars::branch_buffer, pars::branch_split);
	tree->Branch("proton.phi", &protons.phi, pars::branch_buffer, pars::branch_split);
	tree->Branch("proton.eta", &protons.eta, pars::branch_buffer, pars::branch_split);
	tree->Branch("proton.dca", &protons.dca, pars::branch_buffer, pars::branch_split);
	tree->Branch("proton.nsigma", &protons.nsigma, pars::branch_buffer, pars::branch_split);
	tree->Branch("proton.beta", &protons.beta, pars::branch_buffer, pars::branch_split);
	tree->Branch("proton.charge", &protons.charge, pars::branch_buffer, pars::branch_split);

	tree->Branch("pion.pt", &pions.pt, pars::branch_buffer, pars::branch_split);
	tree->Branch("pion.phi", &pions.phi, pars::branch_buffer, pars::branch_split);
	tree->Branch("pion.eta", &pions.eta, pars::branch_buffer, pars::branch_split);
	tree->Branch("pion.dca", &pions.dca, pars::branch_buffer, pars::branch_split);
	tree->Branch("pion.nsigma", &pions.nsigma, pars::branch_buffer, pars::branch_split);
	tree->Branch("pion.beta", &pions.beta, pars::branch_buffer, pars::branch_split);
	tree->Branch("pion.charge", &pions.charge, pars::branch_buffer, pars::branch_split);

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
	track_cut_hist->GetXaxis()->SetBinLabel(1, "Expected");
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

	de_dx_pq_hist = new TH2F("dedx_pq_pid", "Dedx PID", 1000, -3, 3, 1000, 0, 0.5e-4);
	beta_pq_hist = new TH2F("beta_pq_pid", "Beta PID", 1000, -3, 3, 1000, 0, 5);

	return kStOK;
}


Int_t TreeMaker::Make() {
	events_read++;
	event_cut_hist->Fill("Original", 1);

	StMuEvent* mu_event = mudst_maker->muDst()->event();
	if(is_bad_event(mu_event)) { return kStOk; }

}


Int_t TreeMaker::Finish() {
	//
}



// Doers
bool TreeMaker::is_bad_event(StMuEvent *mu_event) {
	if(!mu_event) { return true; }
	event_cut_hist->Fill("Is muEvent", 1);

	// Check for good trigger
	vector<int> good_triggers = pars::triggers[energy];
	bool good_trig = false;
	for(int trig_index = 0; trig_index < (int)pars::triggers[energy].size(); trig_index++) {
		if(mu_event->triggerIdCollection().nominal().isTrigger(pars::triggers[energy][trig_index])) {
			good_trig = true;
			break;
		}
	}
	if(!good_trig) { return true; }
	event_cut_hist->Fill("Good Trigger", 1);


    // Check if run number is good
    event.run_num = mu_event->runId();
    vector<int> bad_runs_energy = pars::bad_runs[energy];
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

	// Check vertex is within pars::vz_max[energy] cm of detector center along beam pipe
	if(fabs(event.vz) > pars::vz_max[energy]) { return true; }
	event_cut_hist->Fill("Good Vz", 1);

	// Check that vertex is within x cm radially (x-y plane) of detector axis
	if(sqrt(pow(event.vx, 2) + pow(event.vy + pars::vy_offset[energy], 2)) > pars::vr_max[energy]) {
		return true;
	}
	event_cut_hist->Fill("Good Vr", 1);

	// On old tapes, no-vertex gets reported as VtxPosition=(0,0,0)
	if(fabs(event.vx) < pars::vertex_min &&
			fabs(event.vy) < pars::vertex_min &&
			fabs(event.vz) < pars::vertex_min) {
		return true;
	}
	event_cut_hist->Fill("Vertex Non-Zero", 1);


	// Filter out events with disagreement between vpd and vertex reconstruction.
	if(energy >= 39) {
		if(mudst_maker->muDst()->btofHeader()) {
			VpdVzPos    =  mMuDstMaker->muDst()->btofHeader()->vpdVz();
			VertexZPos  =  muEvent-> primaryVertexPosition().z();
			if(fabs(VpdVzPos-VertexZPos) > 3) return kStOK; // for 39,62 GeV
		} else {
			return kStOK; }
	}


	return false;  // If all above checks are passed, event is good
}

void track_loop(StMuEvent *mu_event) {

}
