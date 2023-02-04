#include "Flattener.h"


// Structors
Flattener::Flattener() {
	phi_file_name = "";
	phi_file = NULL;
	ep_file_name = "";
	ep_file = NULL;

	cent_bins = vector<int>(10);
	iota(cent_bins.begin(), cent_bins.end(), -1);  // Populate cent bins with -1, 0, 1, 2, ..., 8
	phi_types = { "protons", "non-protons" };
}

Flattener::Flattener(string phi_name) {
	phi_file_name = phi_name;
	phi_file = NULL;
	ep_file_name = "";
	ep_file = NULL;

	cent_bins = vector<int>(10);
	iota(cent_bins.begin(), cent_bins.end(), -1);  // Populate cent bins with -1, 0, 1, 2, ..., 8
	phi_types = { "protons", "non-protons" };
}

Flattener::Flattener(string phi_name, string ep_name) {
	phi_file_name = phi_name;
	phi_file = NULL;
	ep_file_name = ep_name;
	ep_file = NULL;

	cent_bins = vector<int>(10);
	iota(cent_bins.begin(), cent_bins.end(), -1);  // Populate cent bins with -1, 0, 1, 2, ..., 8
	phi_types = { "protons", "non-protons" };
}

Flattener::~Flattener() {
	// Nothing
}


// Setters


// Getters


// Doers

void Flattener::test() {
	cout << "Flattener working" << endl;
}

// Initialize outfile and sin/cos_terms for getting phi coefficients
void Flattener::init_phi_flattener() {
	phi_file = new TFile(phi_file_name.data(), "UPDATE");
	init_phi_terms();
}


// Initialize input phi coefficient file and output event plane coefficient file
void Flattener::init_ep_flattener() {
	phi_file = new TFile(phi_file_name.data(), "READ");
	ep_file = new TFile(ep_file_name.data(), "UPDATE");
	init_phi_terms();
}

// Initialize input phi and event plane coefficient files
void Flattener::init_treemaker() {
	phi_file = new TFile(phi_file_name.data(), "READ");
	ep_file = new TFile(ep_file_name.data(), "READ");
	init_phi_terms();
}

// Initialize phi_sin/cos_terms
void Flattener::init_phi_terms() {
	for (string phi_type : phi_types) {
		for (int cent_bin : cent_bins) {
			for (int eta_bin = 0; eta_bin < eta_bins; eta_bin++) {
				phi_sin_terms[phi_type][cent_bin].push_back({});
				phi_cos_terms[phi_type][cent_bin].push_back({});
			}
		}
	}
}

// Fill TProfiles with harmonic terms of phi
void Flattener::calc_phi_terms(string particle_type, int cent_bin, int eta_bin, int run_key, float phi) {
	if (phi_sin_terms[particle_type][cent_bin][eta_bin].count(run_key) < 1) {
		phi_file->cd();
		string sin_name = "sine_terms_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
		phi_sin_terms[particle_type][cent_bin][eta_bin][run_key] = new TProfile(sin_name.data(), "Sine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
		string cos_name = "cosine_terms_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
		phi_cos_terms[particle_type][cent_bin][eta_bin][run_key] = new TProfile(cos_name.data(), "Cosine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
	}
	for (int n = n_harmonic_low; n <= n_harmonic_high; n++) {
		phi_sin_terms[particle_type][cent_bin][eta_bin][run_key]->Fill(n, sin(n * phi));
		phi_cos_terms[particle_type][cent_bin][eta_bin][run_key]->Fill(n, cos(n * phi));
	}
	return;
}


// Write and close phi output file
void Flattener::write_phi() {
	phi_file->Write();
	phi_file->Close();
}

// Write and close event plane output file, close input phi file
void Flattener::write_ep() {
	ep_file->Write();
	ep_file->Close();
	phi_file->Close();
}

// Close input event plane and phi files
void Flattener::close_phi_ep() {
	ep_file->Close();
	phi_file->Close();
}


// Calculate eta bin for given eta value
int Flattener::get_eta_bin(float eta) {
	return 1;
	if (eta == eta_max)  return eta_bins - 1;
	float eta_range = eta_max - eta_min;
	return int((eta - eta_min) / eta_range * eta_bins);
}


// Get map key for run_num. For now just truncate last digit of run
int Flattener::get_run_bin_key(int run_num) {
	return int(run_num / 10);
}