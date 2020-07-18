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
void ParticleVars::clear() {
	pt.clear();
	phi.clear();
	eta.clear();
	dca.clear();
	nsigma.clear();
	beta.clear();
	charge.clear();
}
