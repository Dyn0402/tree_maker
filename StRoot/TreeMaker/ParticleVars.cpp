/*
 * ParticleVars.cpp
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */


#include "ParticleVars.h"


// Structors
ParticleVars::ParticleVars() {
	clear();
}

ParticleVars::~ParticleVars() {
	// Nothing
}



// Doers
void ParticleVars::add_event(float pt, float phi, float eta, float dca, float dca_z, float nsigma, float beta, short charge, short nhits_fit) {
	this->pt.push_back(pt);
	this->phi.push_back(phi);
	this->eta.push_back(eta);
	this->dca.push_back(dca);
	this->dca_z.push_back(dca_z);
	this->nsigma.push_back(nsigma);
	this->beta.push_back(beta);
	this->charge.push_back(charge);
	this->nhits_fit.push_back(nhits_fit);
}

void ParticleVars::clear() {
	pt.clear();
	phi.clear();
	eta.clear();
	dca.clear();
	dca_z.clear();
	nsigma.clear();
	beta.clear();
	charge.clear();
	nhits_fit.clear();
}
