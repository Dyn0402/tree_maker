/*
 * FlattenerPhiEp.cpp
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */


#include "FlattenerPhiEp.h"

ClassImp(FlattenerPhiEp)


FlattenerPhiEp::FlattenerPhiEp(StMuDstMaker *maker) : StMaker("FlattenerPhiEp") {
	muDst_maker = maker;
	muDst = NULL;

	picoDst_maker = NULL;
	picoDst = NULL;

	out_file_name = "";
	out_file = NULL;

	events_read = 0;
	events_processed = 0;
	energy = 0;
	bes_phase = 1;
	ref_num = 3;

	refmultCorrUtil = new StRefMultCorr(("refmult" + to_string(ref_num)).data());

	pars.set_energy_bes(energy, bes_phase);
}

FlattenerPhiEp::FlattenerPhiEp(StMuDstMaker *maker, string name, int energy_in, int bes_phase) : StMaker("FlattenerPhiEp") {
	muDst_maker = maker;
	muDst = NULL;

	picoDst_maker = NULL;
	picoDst = NULL;

	out_file_name = name;
	out_file = NULL;

	events_read = 0;
	events_processed = 0;
	energy = energy_in;
	this->bes_phase = bes_phase;
	ref_num = 3;

	refmultCorrUtil = new StRefMultCorr(("refmult" + to_string(ref_num)).data());

	pars.set_energy_bes(energy, bes_phase);
}


// Structors
FlattenerPhiEp::FlattenerPhiEp(StPicoDstMaker *maker) : StMaker("FlattenerPhiEp") {
	picoDst_maker = maker;
	picoDst = NULL;

	muDst_maker = NULL;
	muDst = NULL;

	out_file_name = "";
	out_file = NULL;

	events_read = 0;
	events_processed = 0;
	energy = 0;
	bes_phase = 1;
	ref_num = 3;

	refmultCorrUtil = new StRefMultCorr(("refmult" + to_string(ref_num)).data());

	pars.set_energy_bes(energy, bes_phase);
}

FlattenerPhiEp::FlattenerPhiEp(StPicoDstMaker *maker, string name, int energy_in, int bes_phase) : StMaker("FlattenerPhiEp") {
	picoDst_maker = maker;
	picoDst = NULL;

	muDst_maker = NULL;
	muDst = NULL;

	out_file_name = name;
	out_file = NULL;

	events_read = 0;
	events_processed = 0;
	energy = energy_in;
	this->bes_phase = bes_phase;
	ref_num = 3;

	refmultCorrUtil = new StRefMultCorr(("refmult" + to_string(ref_num)).data());

	pars.set_energy_bes(energy, bes_phase);
}

FlattenerPhiEp::~FlattenerPhiEp() {
	// Nothing
}



// Setters
void FlattenerPhiEp::set_energy(int energy_in) {
	energy = energy_in;
}

void FlattenerPhiEp::set_out_file_name(string name) {
	out_file_name = name;
}



// St Doers
Int_t FlattenerPhiEp::Init() {
	out_file = new TFile(out_file_name.data(),"UPDATE") ;

	for (string phi_type : phi_types) {
		for (int cent_bin : cent_bins) {
			for (int eta_bin = 0; eta_bin < eta_bins; eta_bin++) {
				string name = "phi_dist_" + phi_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin);
				phi_dists[phi_type][cent_bin].push_back(new TH1D(name.data(), "Phi_Dist", 1000, 0, 2 * M_PI);
			}
		}
	}

	return kStOK;
}


Int_t FlattenerPhiEp::Make() {
	events_read++;

	event.clear();  // Clear event objects before processing new event

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

	events_processed++;

	return kStOk;
}


//Int_t FlattenerPhiEp::Finish() {
//	cout << endl;
//	cout << "Finishing and writing histograms to file... " << endl;
//	cout << endl;
//
//	out_file->Write();
//	out_file->Close();
//
//	cout <<"\n ======> Finished <======"<<endl;
//	cout<<" Acutal #Events Read = " << events_read << endl ;
//	cout<<" Acutal #Events Processed = " << events_processed << endl ;
//
//	return kStOk;
//}



// Doers
//bool FlattenerPhiEp::is_bad_event(StMuEvent *mu_event) {
//	if(!mu_event) { return true; }
//	event_cut_hist->Fill("Is dstEvent", 1);
//
//	// Check for good trigger
//	vector<int> good_triggers = pars.triggers;
//	bool good_trig = false;
//	for(int trig_index = 0; trig_index < (int)pars.triggers.size(); trig_index++) {
//		if(mu_event->triggerIdCollection().nominal().isTrigger(pars.triggers[trig_index])) {
//			good_trig = true;
//			break;
//		}
//	}
//	if(!good_trig) { return true; }
//	event_cut_hist->Fill("Good Trigger", 1);
//
//
//    // Check if run number is good
//    event.run_num = mu_event->runId();
//    vector<int> bad_runs_energy = pars.bad_runs;
//    int num_bad_runs = (int) bad_runs_energy.size();
//    for(int bad_run_index = 0; bad_run_index < num_bad_runs; bad_run_index++) {
//    	if(event.run_num == bad_runs_energy[bad_run_index]) {
//    		return true;
//    	}
//    }
//    if(energy == 14) {
//    	if(event.run_num <= pars.min_14GeV_run) { return true; }
//    }
//    event_cut_hist->Fill("Good Run", 1);
//
//    // Get x,y,z components of primary vertex
//	event.vx = mu_event->primaryVertexPosition().x();
//	event.vy = mu_event->primaryVertexPosition().y();
//	event.vz = mu_event->primaryVertexPosition().z();
//
//	// Check vertex is within pars.vz_max cm of detector center along beam pipe
//	if(fabs(event.vz) > pars.vz_max) { return true; }
//	event_cut_hist->Fill("Good Vz", 1);
//
//	// Check that vertex is within x cm radially (x-y plane) of detector axis
//	if(sqrt(pow(event.vx, 2) + pow(event.vy + pars.vy_offset, 2)) > pars.vr_max) {
//		return true;
//	}
//	event_cut_hist->Fill("Good Vr", 1);
//
//	// On old tapes, no-vertex gets reported as VtxPosition=(0,0,0)
//	if(fabs(event.vx) < pars.vertex_min &&
//			fabs(event.vy) < pars.vertex_min &&
//			fabs(event.vz) < pars.vertex_min) {
//		return true;
//	}
//	event_cut_hist->Fill("Vertex Non-Zero", 1);
//
//	// Filter out events with disagreement between vpd and vertex reconstruction.
//	if(pars.vpd_vz_max_diff > 0) {  // -1 error code
//		if(muDst->btofHeader()) {
//			float vpd_vz = muDst->btofHeader()->vpdVz();
//			if(fabs(vpd_vz - event.vz) > pars.vpd_vz_max_diff) {
//				return true;
//			}
//		} else {
//			return true;
//		}
//	}
//	event_cut_hist->Fill("Good VPD Vz", 1);
//
//
//	// Add other event variables to event
//	event.event_id = mu_event->eventId();
//	event.refmult = mu_event->refMult();
//	event.btof_multi = mu_event->btofTrayMultiplicity();
//
//	return false;  // If all above checks are passed, event is good
//}


//bool FlattenerPhiEp::is_bad_event(StPicoEvent *pico_event) {
//	if(!pico_event) { return true; }
//	event_cut_hist->Fill("Is dstEvent", 1);
//
//	// Check for good trigger
//	vector<int> good_triggers = pars.triggers;
//	bool good_trig = false;
//	for(int trig_index = 0; trig_index < (int)pars.triggers.size(); trig_index++) {
//		if(pico_event->isTrigger(pars.triggers[trig_index])) {
//			good_trig = true;
//			break;
//		}
//	}
//	if(!good_trig) { return true; }
//	event_cut_hist->Fill("Good Trigger", 1);
//
//
//    // Check if run number is good
//    event.run_num = pico_event->runId();
//    vector<int> bad_runs_energy = pars.bad_runs;
//    int num_bad_runs = (int) bad_runs_energy.size();
//    for(int bad_run_index = 0; bad_run_index < num_bad_runs; bad_run_index++) {
//    	if(event.run_num == bad_runs_energy[bad_run_index]) {
//    		return true;
//    	}
//    }
//    event_cut_hist->Fill("Good Run", 1);
//
//    // Get x,y,z components of primary vertex
//	event.vx = pico_event->primaryVertex().X();
//	event.vy = pico_event->primaryVertex().Y();
//	event.vz = pico_event->primaryVertex().Z();
//
//	// Check vertex is within pars.vz_max cm of detector center along beam pipe
//	if(fabs(event.vz) > pars.vz_max) { return true; }
//	event_cut_hist->Fill("Good Vz", 1);
//
//	// Check that vertex is within x cm radially (x-y plane) of detector axis
//	if(sqrt(pow(event.vx, 2) + pow(event.vy + pars.vy_offset, 2)) > pars.vr_max) {
//		return true;
//	}
//	event_cut_hist->Fill("Good Vr", 1);
//
//	// On old tapes, no-vertex gets reported as VtxPosition=(0,0,0)
//	if(fabs(event.vx) < pars.vertex_min &&
//			fabs(event.vy) < pars.vertex_min &&
//			fabs(event.vz) < pars.vertex_min) {
//		return true;
//	}
//	event_cut_hist->Fill("Vertex Non-Zero", 1);
//
//	// Filter out events with disagreement between vpd and vertex reconstruction.
//	if(pars.vpd_vz_max_diff > 0) {  // -1 ignore code
//		float vpd_vz = pico_event->vzVpd();
//		if(fabs(vpd_vz - event.vz) > pars.vpd_vz_max_diff) {
//			return true;
//		}
//	}
//	event_cut_hist->Fill("Good VPD Vz", 1);
//
//
//	// Add other event variables to event
//	event.event_id = pico_event->eventId();
//	event.refmult = pico_event->refMult();
//	event.btof_match = pico_event->nBTOFMatch();
//	event.btof_multi = pico_event->btofTrayMultiplicity();
//	event.refmult2 = pico_event->refMult2();
//	event.refmult3 = pico_event->refMult3();
//
//
//	return false;  // If all above checks are passed, event is good
//}


void FlattenerPhiEp::track_loop(StMuEvent *mu_event) {
	int num_primary = muDst->primaryTracks()->GetEntries();
	StMuTrack *track, *track_glob;

	int index_2g, nHitsFit;
	float dca, dca_prim, eta, rapidity, pt, nsigmapr, phi;
	float nsigmapr_eff;
	double ratio; // Important that this is double, 13/25 = 0.52 = cut!!!
	double beta, p, m;
	short charge;

	for (int track_index = 0; track_index < num_primary; track_index++) {  // Do refmult counting to get centrality
		track = (StMuTrack*)muDst->primaryTracks(track_index);

		// Initial track cuts
		if (!track) continue;  // Check that track not NULL

		if (track->vertexIndex() != 0) continue;  // Check that vertex index is zero

		index_2g = track->index2Global();
		if (index_2g < 0) continue;  // Check that global index non negative

		track_glob = (StMuTrack*)muDst->globalTracks(index_2g);

		if (track->flag() < 0) continue;  // Check primary track flag, still unsure what it is

		if (track_glob->flag() < 0) continue;  // Check global track flag, still unsure what it is

		charge = track->charge();
		if (fabs(charge) != 1) continue;  // Eliminates neutral/exotic particles

		// Get main track variables
		p = track->p().mag();
		pt = track->pt();
		eta = track->eta();
		phi = track->phi();  if (phi < 0) { phi += 2 * M_PI; }
		dca_prim = track->dca().mag();
		nsigmapr = track->nSigmaProton();
		nsigmapr_eff = nsigmapr;
		if (energy == 27) { nsigmapr_eff *= 2; }  // BES I 27GeV calibration issue, have to scale nsigmapr by 2

		nHitsFit = track_glob->nHitsFit();

		beta = track->btofPidTraits().beta();
		m = (beta > 1.e-5) ? p * p * (1. / beta / beta - 1.) : -999;

		if (fabs(eta) > 0.5 && fabs(eta) < 1. && dca_prim <= 3. && nHitsFit >= 10 && p >= 1.e-10) event.refmult2++;
		if (fabs(eta) < 1. && nHitsFit >= 10 && dca_prim <= 3. && nsigmapr_eff < -3. && m < 0.4 && p >= 1.e-10) event.refmult3++;
	}

	// Get centrality bin for event from ref_multn value
	refmultCorrUtil->init(event.get_run());
	refmultCorrUtil->initEvent((int)event.get_refn(), (double)event.get_vz());
	int cent9_corr = refmultCorrUtil->getCentralityBin9();

	int eta_bin;

	for (int track_index = 0; track_index < num_primary; track_index++) {  // Get phi distribution
		track = (StMuTrack*)muDst->primaryTracks(track_index);

		// Initial track cuts
		if (!track) continue;  // Check that track not NULL

		if (track->vertexIndex() != 0) continue;  // Check that vertex index is zero

		index_2g = track->index2Global();
		if (index_2g < 0) continue;  // Check that global index non negative

		track_glob = (StMuTrack*)muDst->globalTracks(index_2g);

		if (track->flag() < 0) continue;  // Check primary track flag, still unsure what it is

		if (track_glob->flag() < 0) continue;  // Check global track flag, still unsure what it is

		charge = track->charge();
		if (fabs(charge) != 1) continue;  // Eliminates neutral/exotic particles

		// Get main track variables
		p = track->p().mag();
		pt = track->pt();
		eta = track->eta();
		phi = track->phi();  if (phi < 0) { phi += 2 * M_PI; }
		dca = track->dcaGlobal().mag();
		dca_prim = track->dca().mag();
		nsigmapr = track->nSigmaProton();
		nsigmapr_eff = nsigmapr;
		if (energy == 27) { nsigmapr_eff *= 2; }  // BES I 27GeV calibration issue, have to scale nsigmapr by 2

		nHitsFit = track_glob->nHitsFit();

		beta = track->btofPidTraits().beta();
		m = (beta > 1.e-5) ? p * p * (1. / beta / beta - 1.) : -999;

		// Cut on ratio of nHitsFit to nHitsPossible
		ratio = (double)nHitsFit / (double)track_glob->nHitsPoss();
		if (ratio < 0.52) continue;
		if (ratio > 1.05) continue;

		// Fill phi distributions
		if (nHitsFit > 15 && dca < 2.0 && fabs(eta) < 1.0 && pt > 0.2 && pt < 2.) {
			eta_bin = get_eta_bin(eta);

			rapidity = log((sqrt(pow(pars.m_proton, 2) + pow(pt, 2) * pow(cosh(eta), 2)) + pt * sinh(eta)) / sqrt(pow(pars.m_proton, 2) + pow(pt, 2)));
			if (track->nHitsDedx() > 5 && dca < 1.0 && pt >= 0.3 && fabs(nsigmapr_eff) < 2.0 && ((m > 0.6 && m < 1.2) || m == -999) && fabs(rapidity) <= 0.5) {
				phi_dists["protons"][cent9_corr][eta_bin]->Fill(phi);
			}
			else {
				phi_dists["non-protons"][cent9_corr][eta_bin]->Fill(phi);
			}
		}
	}

}


void FlattenerPhiEp::track_loop(StPicoEvent *pico_event) {
	int num_tracks = picoDst->numberOfTracks();
	StPicoTrack *track;

	// Get centrality bin for event from ref_multn value
	refmultCorrUtil->init(event.get_run());
	refmultCorrUtil->initEvent((int)event.get_refn(), (double)event.get_vz());
	int cent9_corr = refmultCorrUtil->getCentralityBin9();

	int nHitsFit;
	float dca, eta, rapidity, pt, nsigmapr, phi;
	float nsigmapr_eff;
	double ratio; // Important that this is double, 13/25 = 0.52 = cut!!!
	double beta, p, m;
	short charge;

	int eta_bin;

	for(int track_index = 0; track_index < num_tracks; track_index++) {
		track = (StPicoTrack*) picoDst->track(track_index);

		// Initial track cuts
		if(!track) continue;  // Check that track not NULL

		if(!track->isPrimary()) continue;  // Check track is primary track

		charge = track->charge();
		if(fabs(charge) != 1) continue;  // Eliminates neutral/exotic particles

		// Get main track variables
		p = track->pMom().Mag();
		pt = track->pMom().Perp();
		eta = track->pMom().PseudoRapidity();
		phi = track->pMom().Phi();  if(phi < 0) { phi += 2*M_PI; }
		dca = track->gDCA(event.vx, event.vy, event.vz);
		nsigmapr = track->nSigmaProton();
		nsigmapr_eff = nsigmapr;
		if(energy == 27) { nsigmapr_eff *= 2; }  // BES I 27GeV calibration issue, have to scale nsigmapr by 2

		nHitsFit = track->nHitsFit();

		int btof_pid_traits_index = track->bTofPidTraitsIndex();
		if(btof_pid_traits_index >= 0) {
			StPicoBTofPidTraits *btof_pid_traits = picoDst->btofPidTraits(btof_pid_traits_index);
			beta = btof_pid_traits->btofBeta();
			m = (beta > 1.e-5) ? p*p*(1./beta/beta - 1.) : -999;
		}

		// Cut on ratio of nHitsFit to nHitsPossible
		ratio = (double) nHitsFit / (double) track->nHitsMax();
		if(ratio < 0.52) continue;
		if(ratio > 1.05) continue;

		// Fill phi distributions
		if(nHitsFit > 15 && dca < 2.0 && fabs(eta) < 1.0 && pt > 0.2 && pt < 2.) {
			eta_bin = get_eta_bin(eta);

			rapidity = log((sqrt(pow(pars.m_proton, 2) + pow(pt, 2) * pow(cosh(eta), 2)) + pt * sinh(eta)) / sqrt(pow(pars.m_proton, 2) + pow(pt, 2)));
			if (track->nHitsDedx() > 5 && dca < 1.0 && pt >= 0.3 && fabs(nsigmapr_eff) < 2.0 && ((m > 0.6 && m < 1.2) || m == -999) && fabs(rapidity) <= 0.5) {
				phi_dists["protons"][cent9_corr][eta_bin]->Fill(phi);
			}
			else {
				phi_dists["non-protons"][cent9_corr][eta_bin]->Fill(phi);
			}
		}
	}
}


// Calculate eta bin for given eta value
int FlattenerPhiEp::get_eta_bin(float eta) {
	if (eta == eta_max)  return eta_bins - 1;
	float eta_range = eta_max - eta_min;
	return int((eta - eta_min) / eta_range * eta_bins);
}