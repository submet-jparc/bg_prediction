////////////////////////////////////////////////////////////////////////////////
//
//   bg_pred.C
//
//   Script to predict coincidence counts and compare with the actual.
//
//   How to use:
// - In beam-on  case, `root bg_pred.C\(1\)`
// - In beam-off case, `root bg_pred.C\(0\)`
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



void bg_pred(const bool isbeamon = false, const unsigned int region_mask = 131071)
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
	const double bunch0from = 238.;
	const double bunch0to   = 380.;
	const double bunch1from = 716.;
	const double bunch1to   = 856.;
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
	// ┌───┬───┐
	// │ B │ C │ x: whether pointing or not
	// ├───┼───┤
	// │ A │ D │ y: dT(max,min) for us it will be just dT
	// └───┴───┘ 
	//------------------------------------------------------------------------------
	unsigned long int A = 0;
	unsigned long int B = 0;
	unsigned long int C = 0;
	unsigned long int D = 0;
	TH1D* h1dt = new TH1D("h1dt", "dt Distribution", 256, -4095., 4095.);
	TH1D* h1t0 = new TH1D("h1t0", "t0 Distribution", 256,     0., 4095.);
	TH1D* h1t1 = new TH1D("h1t1", "t1 Distribution", 256,     0., 4095.);
	TH2D* h2t    = new TH2D("h2t"   , "t1 vs t0"  , 256, 0. , 4095. , 256, 0. , 4095. );
	TH2D* h2abcd = new TH2D("h2abcd", "ABCD Plane",   2, -.5,     .5,   2, -.5,     .5);
	TH2D* h2map  = new TH2D("h2map" , "Coin Found",  10, -.5,    9.5,   8, -.5,    7.5);



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
//		if ( *std::max_element(event_volt_rms_abe -> begin(), event_volt_rms_abe -> end()) > 2. ) continue;
		if ( *std::max_element(event_volt_rms_abe -> begin(), event_volt_rms_abe -> end()) > 1.3 ) continue;

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
			if ( pulse_area -> at(j) * ADCToV * SampleToNs < pulse_spe_area -> at(j) - 2. * pulse_spe_width -> at(j) ) continue;
			if ( pulse_area -> at(j) * ADCToV * SampleToNs > pulse_spe_area -> at(j) + 2. * pulse_spe_width -> at(j) ) continue;

			// 1 < width < 10
			if ( pulse_time_fwhm -> at(j) < 2 ) continue;
			if ( pulse_time_fwhm -> at(j) > 9 ) continue;

			//--------------------------------------
			// Are you in layer0 or 1?
			//--------------------------------------
			if ( pulse_il -> at(j) == 0 )
			{
				ll0 . push_back(j);
				h1t0 -> Fill(pulse_time_is -> at(j));
			}
			else
			{
				ll1 . push_back(j);
				h1t1 -> Fill(pulse_time_is -> at(j));
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
			for ( unsigned short int k = 0; k < ll1 . size(); k++ )
			{
				//--------------------------------------
				// Calculate delta t
				//--------------------------------------
				dt = pulse_time_is -> at(ll1[k]) - pulse_time_is -> at(ll0[j]);

				//--------------------------------------
				// Fill hists
				//--------------------------------------
				h1dt -> Fill(dt);
				h2t -> Fill(pulse_time_is -> at(ll0[j]), pulse_time_is -> at(ll1[k]));

				//--------------------------------------
				// Aligned?
				//--------------------------------------
				if ( pulse_ix -> at(ll0[j]) == pulse_ix -> at(ll1[k])
				  && pulse_iy -> at(ll0[j]) == pulse_iy -> at(ll1[k]) )
				{
					//--------------------------------------
					// In time?
					//--------------------------------------
					if ( TMath::Abs(dt) < dtmax )
					{
						if ( !isbeamon )
						{
							A++;
							h2map -> Fill(pulse_ix -> at(ll0[j]), pulse_iy -> at(ll0[j]));
						}
					}
					else
					{
						B++;
					}
				}
				else
				{
					//--------------------------------------
					// In time?
					//--------------------------------------
					if ( TMath::Abs(dt) < dtmax )
					{
						D++;
					}
					else
					{
						C++;
					}
				}
			}
		}
	}



	//------------------------------------------------------------------------------
	// Fill ABCD hist
	//------------------------------------------------------------------------------
	h2abcd -> SetBinContent(1, 1, A);
	h2abcd -> SetBinContent(1, 2, B);
	h2abcd -> SetBinContent(2, 1, D);
	h2abcd -> SetBinContent(2, 2, C);


	//------------------------------------------------------------------------------
	// Prediction result
	//------------------------------------------------------------------------------
	p = 1. * B * D / C;
	dp = TMath::Sqrt(p * (p+D+B) / C);



	//------------------------------------------------------------------------------
	// Print info
	//------------------------------------------------------------------------------
	std::cout << "[Info]"                                       << std::endl;
	std::cout << "Clean entries      = " << nEntriesClean << " (" << 100. * nEntriesClean / nEntries << "%)" << std::endl;
	std::cout << "dtmax              = " << dtmax << " samples" << std::endl;
	std::cout << "ABCD prediction    = " << p << " +/- " << dp  << std::endl;
	std::cout << "Expected nBG per evt per sample = " << p / nEntriesClean / trange / SampleToNs << " +/- " << dp / nEntriesClean / trange / SampleToNs << std::endl;
	if ( !isbeamon ) std::cout << "Yields in region A = " << A << std::endl;



	//------------------------------------------------------------------------------
	// Draw plot
	//------------------------------------------------------------------------------
	TCanvas* c1 = new TCanvas("c1", "Prediction", 1500, 1000);
	c1 -> Divide(3, 2);
	c1 -> cd(1) -> SetLogz();
	h2t -> Draw("colz");
	c1 -> cd(2);
	h1dt -> Draw();
	c1 -> cd(3) -> SetLogz();
	h2abcd -> Draw("coltext");
	c1 -> cd(4);
	h1t0 -> Draw();
	c1 -> cd(5);
	h1t1 -> Draw();
	c1 -> cd(6);
	if ( !isbeamon ) h2map -> Draw("colztext");



	//------------------------------------------------------------------------------
	// Finalize
	//------------------------------------------------------------------------------
	return;
}
