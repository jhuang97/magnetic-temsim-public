#pragma once

#ifndef PARAMGROUP_HPP   // only include this file if its not already
#define PARAMGROUP_HPP   // remember that this has been included

using namespace std;
#include <string>   // STD string class

struct opts_4dsave {
	bool save_wvfn;
	bool save_cbed;
	bool save_probe_int_distr; // projected probe intensity distributions
	string wvfn_dir;
	string cbed_dir;
	vectori pos_ix;
	vectori pos_iy;
	int cbed_nxout;
	int cbed_nyout;
};

#endif