#include<iostream>
#include<fstream>
#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<TSpectrum.h>
#include<TH1F.h>
#include<TCanvas.h>
#include<TPad.h>
#include<TFile.h>
#include<time.h>

struct dithermap{
	double prob[10];
	int assign[10];
};

#ifndef __TDC_CH_VALUES_H__
#define __TDC_CH_VALUES_H__

class tdc_ch_values{

   public:

	//FUNCTIONS TO RETURN THE CREATED HISTOGRAMS

	TH1F* get_cdt() {return this->cdt;};
	TH1F* get_DNL() {return this->DNL;};
	TH1F* get_INL() {return this->INL;};
//	TH1F* get_DNL_c() {return this->DNL_c};
//	TH1F* get_INL_c() {return this->INL_c};
//	TH1F* get_corrected() {return this->corrected};
//	TH1F* get_corrected_dith() {return this->corrected_dith};
	TH1F* get_realtimedist() {return this->realtimedist;};
	TH1F* get_mps() {return this->mps;};

	//FUNCTIONS TO FILL AND CREATE THE HISTOGRAMS

	tdc_ch_values(int min,int max); //CONSTRUCTOR: ALLOCATE ALL THE NEEDED HISTOGRAMS

	void fill(unsigned int tdc_value);	
	double get_ps_value(unsigned int tdc_value);
	unsigned int get_dith_value(unsigned int tdc_value);
	int update_histos();

   private:
	void create_map(TH1F* realtime); //CREATES THE MAP FOR PSEUDORANDOM BIN DITHERING
	unsigned int bin_dice(int bin); //RETURNS THE ASSIGNED VALUE OF THE GIVEN BIN
	
	//START AND STOP OF THE BINS CONTAINING HITS
	int start_bin;
	int stop_bin;

	struct dithermap *interpol;
	TH1F *cdt;
	TH1F *DNL;
	TH1F *INL;
//	TH1F *DNL_c;	//THESE SHOULD BE CREATED OUTSIDE OF THE CLASS TO SEE WHETHER THE CORRECTIONS ARE WORKING
//	TH1F *INL_c;
//	TH1F *corrected;
//	TH1F *corrected_dith;

	TH1F* realtimedist;
	TH1F* mps;

};
#endif
