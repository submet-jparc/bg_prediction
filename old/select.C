////////////////////////////////////////////////////////////////////////////////
//
//   How to use:
// - In beam-on  case, `root select.C\(1\)`
// - In beam-off case, `root seelct.C\(0\)`
//
//   Jul 2025, Hoyong Jeong (hoyong5419@korea.ac.kr)
//
////////////////////////////////////////////////////////////////////////////////



//------------------------------------------------------------------------------
// Required headers
//------------------------------------------------------------------------------
#include <iostream>
#include <vector>
#include <algorithm>

#include "TString.h"
#include "TChain.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TMath.h"



void select(const bool isbeamon = false, const unsigned int region_mask = 131071, const unsigned short int imod = 0)
{
	//------------------------------------------------------------------------------
	// Init
	//------------------------------------------------------------------------------
	TString sTemp;
	unsigned long int nEntries      = 0;
	unsigned long int nEntriesClean = 0;
	gInterpreter -> GenerateDictionary("vector<vector<unsigned short>>", "vector");
	const unsigned short int dtmax = 14;
	double dt = 0.;
	double p  = 0.;
	double dp = 0.;
	std::vector<unsigned short int> event_npulses_sel(160, 0);
	unsigned short int event_npulses_sel_total = 0;
	double hmax = 0.;


	//----------------------------------------------------------
	// Region of interest
	//----------------------------------------------------------
	double trange = 0.;
	short int last_region = -1;
	std::vector<bool> region(17);
	//--------------------------------------
	// Bits extraction
	//--------------------------------------
	for ( unsigned short int i = 0; i < 17; i++ )
	{
		region[i] = (region_mask >> i) & 1;
		if ( region[i] )
		{
			last_region = i;
			if      ( i     ==  0 ) trange += 132.305; // Pre-bunch
			else if ( i     == 16 ) trange += 292.895; // After-bunch
			else if ( i % 2 ==  0 ) trange += 336.4  ; // Inter-bunch
			else                    trange += 140.   ; // In-bunch
		}
	}
	if ( last_region == -1 ) return;
	double roiend = 0.;
	if      ( last_region     == 16 ) roiend = 4000.;
	else if ( last_region % 2 ==  0 ) roiend = 232.305        + 238.2* last_region   ;
	else                              roiend = 232.305 + 140. + 238.2*(last_region-1);



	//------------------------------------------------------------------------------
	// Data files
	//------------------------------------------------------------------------------
	TChain* tree = new TChain("tree");


	//----------------------------------------------------------
	// Load files: beam-off case
	//----------------------------------------------------------
	if ( !isbeamon )
	{
		for ( unsigned short int i = 11; i <= 22; i++ )
//		for ( unsigned short int i = 11; i <= 11; i++ )
		{
			sTemp . Form("/data3/submet_exp/e00006/tree/2.1.0_speinfo/r%05hu_spill.root", i);
			tree -> Add(sTemp . Data());
		}
	}


	//----------------------------------------------------------
	// Load files: beam-on case
	//----------------------------------------------------------
	if ( isbeamon )
	{
		for ( unsigned short int i = 22; i <= 56; i++ )
		{
			if ( i == 34 ) continue;
			sTemp . Form("/data3/submet_exp/e00000/tree/2.1.0_bsdspe_info/r%05hu_spill.root", i);
			tree -> Add(sTemp . Data());
		}
	}


	std::cout << "[Intro]" << std::endl;
	std::cout << "Region of interest: ";
	for ( bool bit : region )
	{
		std::cout << (bit ? "1 " : "0 ");
	}
	std::cout << std::endl;
	std::cout << "Length of roi  = " << trange << std::endl;
	std::cout << "The end of roi = " << roiend << std::endl;
	if ( !isbeamon ) std::cout << "Beam off case" << std::endl;
	else             std::cout << "Beam on  case" << std::endl;
	nEntries = tree -> GetEntries();
	std::cout << "Entries = " << nEntries << std::endl;


	//----------------------------------------------------------
	// Branches
	//----------------------------------------------------------
	// Original
	unsigned long int event_id = 0;
	std::vector<unsigned short int>* event_imod          = 0;
	std::vector<unsigned short int>* event_ix            = 0;
	std::vector<unsigned short int>* event_iy            = 0;
	std::vector<unsigned short int>* event_il            = 0;
	std::vector<double>*             event_volt_mean     = 0;
	std::vector<double>*             event_volt_rms      = 0;
	std::vector<unsigned short int>* event_npulses       = 0;
	unsigned short int               event_npulses_total = 0;
	std::vector<double>*             event_volt_mean_abe = 0;
	std::vector<double>*             event_volt_rms_abe  = 0;
	std::vector<unsigned short int>* event_iter_abe      = 0;
	unsigned short int event_spillnum_short = 0;
	unsigned int event_trg_sec  = 0;
	unsigned int event_trg_msec = 0;
	unsigned short int event_trg_type = 0;
	std::vector<unsigned short int>* event_stopaddr = 0;
	std::vector<unsigned short int>* pulse_imod = 0;
	std::vector<unsigned short int>* pulse_ix   = 0;
	std::vector<unsigned short int>* pulse_iy   = 0;
	std::vector<unsigned short int>* pulse_il   = 0;
	std::vector<double>* pulse_time_zc0    = 0;
	std::vector<double>* pulse_time_zc1    = 0;
	std::vector<double>* pulse_time_tc0    = 0;
	std::vector<double>* pulse_time_tc1    = 0;
	std::vector<double>* pulse_time_is     = 0;
	std::vector<double>* pulse_time_peak   = 0;
	std::vector<double>* pulse_volt_height = 0;
	std::vector<double>* pulse_area        = 0;
	std::vector<double>* pulse_time_n      = 0;
	std::vector<double>* pulse_time_mean   = 0;
	std::vector<double>* pulse_time_rms    = 0;
	std::vector<double>* pulse_time_skew   = 0;
	std::vector<double>* pulse_time_kurt   = 0;
	std::vector<double>* pulse_time_fwhm   = 0;
	std::vector<double>* pulse_time_fwtc   = 0;
	std::vector<std::vector<unsigned short int>>* pulse_vec_time = 0;
	std::vector<std::vector<unsigned short int>>* pulse_vec_volt = 0;

	tree -> SetBranchAddress("event_id"            , &event_id            );
	tree -> SetBranchAddress("event_spillnum_short", &event_spillnum_short);
	tree -> SetBranchAddress("event_trg_sec"       , &event_trg_sec       );
	tree -> SetBranchAddress("event_trg_msec"      , &event_trg_msec      );
	tree -> SetBranchAddress("event_trg_type"      , &event_trg_type      );
	tree -> SetBranchAddress("event_stopaddr"      , &event_stopaddr      );
	tree -> SetBranchAddress("event_imod"          , &event_imod          );
	tree -> SetBranchAddress("event_ix"            , &event_ix            );
	tree -> SetBranchAddress("event_iy"            , &event_iy            );
	tree -> SetBranchAddress("event_il"            , &event_il            );
	tree -> SetBranchAddress("event_volt_mean"     , &event_volt_mean     );
	tree -> SetBranchAddress("event_volt_rms"      , &event_volt_rms      );
	tree -> SetBranchAddress("event_npulses"       , &event_npulses       );
	tree -> SetBranchAddress("event_npulses_total" , &event_npulses_total );
	tree -> SetBranchAddress("event_volt_mean_abe" , &event_volt_mean_abe );
	tree -> SetBranchAddress("event_volt_rms_abe"  , &event_volt_rms_abe  );
	tree -> SetBranchAddress("event_iter_abe"      , &event_iter_abe      );
	tree -> SetBranchAddress("pulse_imod"          , &pulse_imod          );
	tree -> SetBranchAddress("pulse_ix"            , &pulse_ix            );
	tree -> SetBranchAddress("pulse_iy"            , &pulse_iy            );
	tree -> SetBranchAddress("pulse_il"            , &pulse_il            );
	tree -> SetBranchAddress("pulse_time_zc0"      , &pulse_time_zc0      );
	tree -> SetBranchAddress("pulse_time_zc1"      , &pulse_time_zc1      );
	tree -> SetBranchAddress("pulse_time_tc0"      , &pulse_time_tc0      );
	tree -> SetBranchAddress("pulse_time_tc1"      , &pulse_time_tc1      );
	tree -> SetBranchAddress("pulse_time_is"       , &pulse_time_is       );
	tree -> SetBranchAddress("pulse_time_peak"     , &pulse_time_peak     );
	tree -> SetBranchAddress("pulse_volt_height"   , &pulse_volt_height   );
	tree -> SetBranchAddress("pulse_area"          , &pulse_area          );
	tree -> SetBranchAddress("pulse_time_n"        , &pulse_time_n        );
	tree -> SetBranchAddress("pulse_time_mean"     , &pulse_time_mean     );
	tree -> SetBranchAddress("pulse_time_rms"      , &pulse_time_rms      );
	tree -> SetBranchAddress("pulse_time_skew"     , &pulse_time_skew     );
	tree -> SetBranchAddress("pulse_time_kurt"     , &pulse_time_kurt     );
	tree -> SetBranchAddress("pulse_time_fwhm"     , &pulse_time_fwhm     );
	tree -> SetBranchAddress("pulse_time_fwtc"     , &pulse_time_fwtc     );
	tree -> SetBranchAddress("pulse_vec_time"      , &pulse_vec_time      );
	tree -> SetBranchAddress("pulse_vec_volt"      , &pulse_vec_volt      );

	// Additional
	std::vector<double>* pulse_spe_area  = 0;
	std::vector<double>* pulse_spe_width = 0;
	tree -> SetBranchAddress("pulse_spe_area" , &pulse_spe_area );
	tree -> SetBranchAddress("pulse_spe_width", &pulse_spe_width);



	//------------------------------------------------------------------------------
	// Define counters and containers
	//------------------------------------------------------------------------------
	TH1D* ha = new TH1D("ha", "Area", 1000, 0, 1000);



	//------------------------------------------------------------------------------
	// Looping over entries
	//------------------------------------------------------------------------------
	std::vector<unsigned short int> ll0;
	std::vector<unsigned short int> ll1;
	for ( unsigned long int i = 0; i < nEntries; i++ )
	{
		tree -> GetEntry(i);
		if ( (i + 1) % 100000 == 0 ) std::cout << i + 1 << "-th events being processed." << std::endl;


		//----------------------------------------------------------
		// Event cleaning
		//----------------------------------------------------------
		//--------------------------------------
		// Too noisy pedestal
		//--------------------------------------
		if ( *std::max_element(event_volt_rms_abe -> begin(), event_volt_rms_abe -> end()) > 2. ) continue;

		//--------------------------------------
		// Large hit in the time window or before?
		//--------------------------------------
		hmax = 0.;
		for ( unsigned short int j = 0; j < pulse_imod -> size(); j++ )
		{
			if ( pulse_time_is -> at(j) < roiend )
			{
				if ( pulse_volt_height -> at(j) > hmax ) hmax = pulse_volt_height -> at(j);
			}
		}
		if ( hmax > 500. ) continue;


		//----------------------------------------------------------
		// Layer 0 pulses / layer 1 pulses, extra cleaning
		//----------------------------------------------------------
		ll0 . clear();
		ll1 . clear();
		std::fill(event_npulses_sel . begin(), event_npulses_sel . end(), 0);
		event_npulses_sel_total = 0;
		for ( unsigned short int j = 0; j < pulse_imod -> size(); j++ )
		{
			//--------------------------------------
			// ROI cut
			//--------------------------------------
			if ( !((region_mask >> GetTimeRegion(pulse_time_is -> at(j))) & 1) ) continue;

			//--------------------------------------
			// SPE pulse selection
			//--------------------------------------
			// SPE area +/- 2 sigma
//			if ( pulse_area -> at(j) * ADCToV * SampleToNs < pulse_spe_area -> at(j) - 2. * pulse_spe_width -> at(j) ) continue;
//			if ( pulse_area -> at(j) * ADCToV * SampleToNs > pulse_spe_area -> at(j) + 2. * pulse_spe_width -> at(j) ) continue;
			if ( pulse_area -> at(j) > 1000. ) continue;

			// 1 < width < 10
			if ( pulse_time_fwhm -> at(j) < 2 ) continue;
			if ( pulse_time_fwhm -> at(j) > 9 ) continue;

			//--------------------------------------
			// Are you in layer0 or 1?
			//--------------------------------------
			if ( pulse_il -> at(j) == 0 )
			{
				ll0 . push_back(j);
			}
			else
			{
				ll1 . push_back(j);
			}

			//--------------------------------------
			// # of selected pulses
			//--------------------------------------
			event_npulses_sel[pulse_imod -> at(j)]++;
			event_npulses_sel_total++;
		}


		//----------------------------------------------------------
		// Event cleaning 2
		//----------------------------------------------------------
		// Large number of after-pulses
		if ( *std::max_element(event_npulses_sel . begin(), event_npulses_sel . end()) > 2 ) continue;

		// Flash bang
		if ( event_npulses_sel_total > 20 ) continue;


		//----------------------------------------------------------
		// Count clean events; how many events survived?
		//----------------------------------------------------------
		nEntriesClean++;


		//----------------------------------------------------------
		// Counting
		//----------------------------------------------------------
		for ( unsigned short int j = 0; j < ll0 . size(); j++ )
		{
			if ( pulse_imod -> at(ll0[j]) == imod )
			{
				ha -> Fill(pulse_area -> at(ll0[j]));
//				cout << "area = "    << pulse_spe_area  -> at(ll0[j]) * NsToSample * VToADC
//				     << ", width = " << pulse_spe_width -> at(ll0[j]) * NsToSample * VToADC << endl;
			}
		}
		for ( unsigned short int j = 0; j < ll1 . size(); j++ )
		{
			if ( pulse_imod -> at(ll1[j]) == imod ) ha -> Fill(pulse_area -> at(ll1[j]));
		}
	}



	//------------------------------------------------------------------------------
	// Draw plot
	//------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas("c1", "Select", 1500, 1000);
	c1 -> Divide(3, 2);
	c1 -> cd(1);
	ha -> Draw();



	//------------------------------------------------------------------------------
	// Finalize
	//------------------------------------------------------------------------------
	return;
}
