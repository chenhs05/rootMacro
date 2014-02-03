#include "tdc_ch_values.h"



tdc_ch_values::tdc_ch_values(int min,int max) //Constructor: allocate all the needed histograms
{

	cdt = new TH1F("cdt","TDC Channel CDT",(max-min)+1,min,max+1);
	DNL = new TH1F("dnl","DNL Evaluation without remapping",(max-min)+1,min,max+1);
	INL = new TH1F("inl","INL Evaluation",(max-min)+1,min,max+1);

	realtimedist= new TH1F("realdelays","real bin delays calculated from CDT",(max-min)+1,min,max+1);
	mps = new TH1F("probable_binsize","probability distribution of binsizes",400,0,200e-12);

	//allocate the map for bin dithering
	interpol = new dithermap[1024];
	//printf("------	%d\n",interpol[1]);
	//printf("------	%d\n",interpol[0].assign[0]);
	//printf("------	%d\n",interpol[0]);
}

void tdc_ch_values::create_map(TH1F* realtime)
{
	FILE* fout;
	fout =fopen("create_map.log","w");


	//memset(interpol,0,sizeof(interpol));	//can't be excused by root
	int i,j;

	for(i=0; i<1024;i++)
	{
//		memset(interpol[i].prob,0,sizeof(interpol[i].prob));
//		memset(interpol[i].assign,0,sizeof(interpol[i].assign));
		for(j=0;j<10;j++)
		{
			interpol[i].prob[0]=0.0;
			interpol[i].assign[0]=0;
		}
	//printf("	i=%d	coming here!!\n",i);
	}

	int nbins;
	double tlow,thigh,tsize;

	int binlow, binhigh;
	double blr, bhr;

	double binsize;

	//binsize=(1.0/(320e6)/512.0);	//Binsize of the bins we map to
	binsize=(1.0/(40e6)/1024.0);	//Binsize of the bins of CAEN TDC we map to
//	binsize=1.0;
	tlow=0;	//lower and higher timeborder of the bin
	thigh=0;

	for (i=0;i<realtime->GetNbinsX();i++)
	{
		tlow=thigh;
		thigh += realtime->GetBinContent(i+1);	//Get the Content of bin+1 because STUPID ROOT USES BIN 0 AS UNDERFLOW!!!! 
		tsize = realtime->GetBinContent(i+1); 
		binlow = (int) (tlow/binsize);
		binhigh = (int) (thigh/binsize);


		blr = binsize - (tlow-binlow*binsize); //much better calculation of the remainders
		bhr = thigh - binhigh*binsize;


		fprintf(fout,"creating dithermap: Real bin %u : size: %E tlow: %E thigh: %E binlow: %u binhigh: %u \n", i,tsize,tlow,thigh, binlow, binhigh);


		if(binlow == binhigh)	//THE BIN IS COMPLETLY WITHIN A BIN IT GETS MAPPED TO SO NO NEED TO CALCULATE THE PROBABILITIES
		{
			interpol[i].prob[0] = 1.0;
			interpol[i].assign[0] = binlow;
		}
		else
		{
			nbins = binhigh-binlow-1;	//how many bins except the start and end bin are occupied?
			interpol[i].prob[0] = blr/tsize;	//assign the probabilities of the low and high bin	
			interpol[i].assign[0] = binlow;
			if(nbins != 0)
			{
				for(j=0;j<nbins;j++)		//assign the probabilites for the rest of the occupied bins
				{
					interpol[i].prob[1+j] = blr/tsize+(j+1)*(tsize-bhr-blr)/(tsize*nbins); //Take care with the probabilities
					interpol[i].assign[1+j] = binlow+1+j;
				}
			}

			interpol[i].prob[nbins+1] = interpol[i].prob[nbins]+bhr/tsize;
			interpol[i].assign[nbins+1] = binhigh;


		}

	}
	fclose(fout);
}

double tdc_ch_values::get_ps_value(unsigned int tdc_value)
{
	double ps1=0;
	int i;
	for(i=start_bin+1;i<=(0x1FF&tdc_value)+1;i++)
	{
		ps1 += realtimedist->GetBinContent(i);
	}
	ps1 = (double)(tdc_value >> 6)/320e6 - ps1;	//SPECIFIC FOR FPGA TDC!!!
	return ps1;
}

//USE RANDOM NUMBER GENERATOR TO GET THE DITHERED BIN NUMBER
unsigned int tdc_ch_values::bin_dice(int bin)
{
	double rnumber;
	rnumber = (double)(rand() % 10000)/10000.0;
	float last_prob;
	last_prob = 0.0;
	for(int i = 0; i<10;i++)
	{
		if(last_prob <= rnumber && interpol[bin].prob[i] > rnumber) 
		{
			return interpol[bin].assign[i];
		}
		last_prob = interpol[bin].prob[i];
	}

	return 0;
}


unsigned int tdc_ch_values::get_dith_value(unsigned int tdc_value)
{
	unsigned int temp;
	temp = 0x3FF & tdc_value;
	temp = bin_dice(temp);
	temp = tdc_value + temp - (0x3FF&tdc_value);
	return temp;

}

int tdc_ch_values::update_histos()
{

	int i,bin_number;
	double sum;
	start_bin = 512;
	stop_bin = 0;

	for(i=1; i <= cdt->GetNbinsX();i++)
	{
		bin_number = i - 1;		//~₢!"☠! ROOT -- Maybe I should stick to the bin number of ROOT but...
		if(cdt->GetBinContent(i)>2 && bin_number < start_bin) start_bin = bin_number;
		if(cdt->GetBinContent(i)>2 && bin_number > stop_bin) stop_bin = bin_number;
	}


	float average_size = 0.0;
	unsigned int event_number;
	event_number = cdt->GetEntries();
	average_size = (double)event_number/(double)(stop_bin-start_bin+1);

	//CALCULATE THE DNL AND THE REALTIME DISTRIBUTION OF THE BINS
	for(i=1; i <= cdt->GetNbinsX();i++)
	{
		DNL->SetBinContent(i, ((cdt->GetBinContent(i) - average_size)/average_size));
		realtimedist->SetBinContent(i,(cdt->GetBinContent(i)*((1.0/40e6)/(double)event_number)));
	}

	//CALCULATE THE BIN WIDTH DISTRIBUTION
	for (i = 1; i <= realtimedist->GetNbinsX();i++)
	{
		if(realtimedist->GetBinContent(i) != 0) mps->Fill(realtimedist->GetBinContent(i));
	}


	//CREATE THE BIN DITHERING MAP
	create_map(realtimedist);


	//CALCULATE THE INL -- ONLY NEEDED FOR VECTOR CORRECTION 	
	sum=0;
	for(i=start_bin+1;i<=stop_bin+1;i++)	//First calculate the INL without remapping. Remapping should be the latest stage!!
	{
		sum = DNL->GetBinContent(i) + sum;
		INL->SetBinContent(i, sum);
	}

	return 0;

}


void tdc_ch_values::fill(unsigned int tdc_value)
{
	unsigned int time;
	time = 0x3FF&tdc_value;
	cdt->Fill(time);
}
