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
void ParticleVars::add_event(float pt, float phi, float eta, float dca, float nsigma, float beta, short charge) {
	this->pt.push_back(pt);
	this->phi.push_back(phi);
	this->eta.push_back(eta);
	this->dca.push_back(dca);
	this->nsigma.push_back(nsigma);
	this->beta.push_back(beta);
	this->charge.push_back(charge);
}

void ParticleVars::clear() {
	pt.clear();
	phi.clear();
	eta.clear();
	dca.clear();
	nsigma.clear();
	beta.clear();
	charge.clear();
}
