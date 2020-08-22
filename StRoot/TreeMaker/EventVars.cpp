/*
 * EventVars.cpp
 *
 *  Created on: Jul 16, 2020
 *      Author: Dylan Neff
 */


#include "EventVars.h"



// Structors
EventVars::EventVars() {
	clear();
}

EventVars::~EventVars() {
	// Nothing
}



// Doers
void EventVars::clear() {
	run_num = 0;
	event_id = 0;
	refmult = 0;
	refmult2 = 0;
	refmult3 = 0;
	btof_multi = 0;
	btof_match = 0;

	vx = 0.;
	vy = 0.;
	vz = 0.;
	qx = 0.;
	qy = 0.;
	dca_xy_avg = 0.;
	dca_xy_err = 0.;
}
