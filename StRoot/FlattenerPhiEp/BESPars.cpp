/*
 * BES_Pars.cpp
 *
 *  Created on: Jul 20, 2020
 *      Author: Dylan Neff
 */



#include "BESPars.h"



// Structors
BESPars::BESPars() {
	energy = 0;
	bes_phase = 1;
	set_pars();
}

BESPars::BESPars(int energy) {
	this->energy = energy;
	bes_phase = 1;
	set_pars();
}

BESPars::BESPars(int energy, int bes_phase) {
	this->energy = energy;
	this->bes_phase = bes_phase;
	set_pars();
}

BESPars::~BESPars() {
	// Nothing
}


// Getters


// Setters
void BESPars::set_energy(int energy) {
	this->energy = energy;
	set_pars();
}

void BESPars::set_bes_phase(int bes_phase) {
	this->bes_phase = bes_phase;
	set_pars();
}

void BESPars::set_energy_bes(int energy, int bes_phase) {
	this->energy = energy;
	this->bes_phase = bes_phase;
	set_pars();
}


// Doers
void BESPars::set_pars() {
	vz_max = vz_max_map[energy][bes_phase];
	vr_max = vr_max_map[energy][bes_phase];
	vy_offset = vy_offset_map[energy][bes_phase];
	if(vpd_vz_max_diff_map.count(energy) > 0 && vpd_vz_max_diff_map[energy].count(bes_phase) > 0) {
		vpd_vz_max_diff = vpd_vz_max_diff_map[energy][bes_phase];
	} else {
		vpd_vz_max_diff = -1;  // -1 ignore code, had to put here since getting compile error with - sign in initialization.
	}
	triggers = triggers_map[energy][bes_phase];
	bad_runs = bad_runs_map[energy][bes_phase];
}

