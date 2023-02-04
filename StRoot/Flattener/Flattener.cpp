#include "Flattener.h"

ClassImp(Flattener)


// Structors
Flattener::Flattener() {
	phi_file_name = "";
	phi_file = NULL;
	ep_file_name = "";
	ep_file = NULL;
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

// Fill TProfiles with harmonic terms of phi
void Flattener::calc_phi_terms(string particle_type, int cent_bin, int eta_bin, int run_key, float phi) {
	//if (sin_terms[particle_type][cent_bin][eta_bin].count(run_key) < 1) {
	//	phi_file->cd();
	//	string sin_name = "sine_terms_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
	//	sin_terms[particle_type][cent_bin][eta_bin][run_key] = new TProfile(sin_name.data(), "Sine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
	//	string cos_name = "cosine_terms_" + particle_type + "_cent_" + to_string(cent_bin) + "_eta_bin_" + to_string(eta_bin) + "_runkey_" + to_string(run_key);
	//	cos_terms[particle_type][cent_bin][eta_bin][run_key] = new TProfile(cos_name.data(), "Cosine Terms", n_harmonic_high - n_harmonic_low + 1, n_harmonic_low - 0.5, n_harmonic_high + 0.5);
	//}
	//for (int n = n_harmonic_low; n <= n_harmonic_high; n++) {
	//	sin_terms[particle_type][cent_bin][eta_bin][run_key]->Fill(n, sin(n * phi));
	//	cos_terms[particle_type][cent_bin][eta_bin][run_key]->Fill(n, cos(n * phi));
	//}
	return;
}


// Calculate eta bin for given eta value
int Flattener::get_eta_bin(float eta) {
	return 1;
	//if (eta == eta_max)  return eta_bins - 1;
	//float eta_range = eta_max - eta_min;
	//return int((eta - eta_min) / eta_range * eta_bins);
}


// Get map key for run_num. For now just truncate last digit of run
int Flattener::get_run_bin_key(int run_num) {
	return int(run_num / 10);
}