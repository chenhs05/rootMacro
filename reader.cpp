#include "TTree.h"
#include "TBranch.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TFile.h"
#include "stdio.h"
#include "TSystem.h"
#include "EventDict.h"
#include "EventType.h"
#include "TFile.h"
#include "TKDE.h"
#include "TMath.h"
#include "Math/RootFinder.h"
#include "Math/WrappedTF1.h"

#include "TSpectrum.h"
#include "TPad.h"
#include "TLine.h"
#include "TStyle.h"
#include "TArrow.h"
#include "TText.h"
#include "TPaveText.h"


int CH_1;  // first channel of the coincidence
int CH_2; // second channel of the coincidence

double tLimit=70; //time difference window to look for coincidence, it also define the range of the time difference plot, the range is [-tLimit,tLimit], unit is 50.2ps
double EHistHigh = 600; // the max for the energy histogram, the histogram range would be [0,EHistHigh], unit is 1.6ns

bool ManuECut=0; //ManuECut=0: energy cut range according to peak position and sigma; ManuECut=1: energy cut range is definde manually

double  NsigmaE=1.5; // the energy cut for the energy peak, the range is -+ NsigmaE*sigma around the 511kev photon peak, unit is coarse counter width, 1.6ns
double  NsigmaT=1.65; // fitting range for the coincidence time difference, the range is -+NsigmaT*sigma around the peak, unit is find counter width, 50.2ps

//define the the energy cut range manually, ManuECut need to be set to 1 
double  CH_1_E_LOW=130; 
double  CH_1_E_HIGH=170; 
double  CH_2_E_LOW=130;
double  CH_2_E_HIGH= 170; 

TF1* func; //function for the kernel density estimator to get the FWHM
FILE* fitResultOut;

Double_t sub_func(Double_t *x, Double_t *par)
{
	return (func->Eval(x[0]) - par[0]/2.0);
}
//get FWHM from kernel density estimator
double get_FWHM()
{
	double function_max;
	function_max=func->GetMaximum();

	double values[2];

	TF1* fwhm_func= new TF1("fwhmfunction",sub_func,-70,70,1);
	fwhm_func->SetParameter(0,function_max);

	printf("Sanity check: value at x_min: %f, function maximum/2.0: %f\n",fwhm_func->Eval(-10.0), function_max/2.0);
	printf("Value of unmodified function: %f \n",func->Eval(-10.0));

	ROOT::Math::WrappedTF1 wf1(*fwhm_func);
	ROOT::Math::RootFinder brf;

	brf.SetFunction(wf1,-10.0,0);
	brf.Solve();
	
	printf("Found Root: %f, status: %d, Iterations: \n",brf.Root(), (int)brf.Status());
	values[0]=brf.Root();
	brf.SetFunction(wf1,0,10.0);
	brf.Solve();
	values[1]=brf.Root();

	printf("Found Root: %f, status: %d, Iterations: \n",brf.Root(),brf.Status());
	printf("FWHM: %f TDC Bins = %f ps FWHM\n",values[1]-values[0],(values[1]-values[0])*50.2 );

	double val_fwhm = (values[1]-values[0])*50.2; 

	//////////////output fit result to txt
	fitResultOut=fopen("fitResult.txt","a");
	fprintf(fitResultOut,"%f	",val_fwhm);
	fclose(fitResultOut);
	/////////////////////////

	return val_fwhm;
}

//get the peak position and sigma of the 511kev photo peak
double getEPeak(double outPara[2],TH1F *totHist)
{
	TF1 *gausFit = new TF1("gausFit","gaus");
	gausFit->SetLineWidth(1);
	double para[16];

	TSpectrum *s = new TSpectrum(2,1.0);

	int nfound = s->Search(totHist,2,"",0.05);
	float *xpeak = s->GetPositionX();
	//if(xpeak[0]<65)	xpeak[0]=xpeak[1]; //if there is a noise peak with low energy , turn this on to take the right peak

	gausFit->SetRange(xpeak[0]-5,xpeak[0]+5);
	totHist->Fit("gausFit","RQ");
	gausFit->GetParameters(&para[0]);

	outPara[0] = para[1];
	outPara[1] = para[2];

	return 0;
}

//main funtion
int main(int argc, char *argv[]){
	if (argc<5){
		printf("use: %s file treename ch#1 ch#2\n",argv[0]);
		return -1;
	}
	else {
		fitResultOut=fopen("fitResult.txt","a");
		fprintf(fitResultOut,"%s	",argv[1]);
		fclose(fitResultOut);
	}


	TFile f(argv[1]);
	TTree *tree=(TTree*)f.Get(argv[2]);
	CH_1 = strtol(argv[3],NULL,0);
	CH_2 = strtol(argv[4],NULL,0);

	if (tree==NULL)
		return -1;

	printf("---------------------------------------------------------------------------\n");

	TObjArray *branches=tree->GetListOfBranches();

	for(int i=0; i<branches->GetEntries(); i++){
		TBranch * thisBranch=(TBranch*)(*branches)[i];
		printf("branch \"%s\" has %lld entries\n",thisBranch->GetName(),thisBranch->GetEntries());
	}

	gStyle->SetOptFit(0);
	gStyle->SetOptStat(111);
	TCanvas *c1 = new TCanvas("c1","Canvas",800,800);

	TH1F *tot_hist_ch1 = new TH1F("tot1","ToT Spectrum 1",(int)EHistHigh,0,EHistHigh);
	TH1F *tot_hist_ch2 = new TH1F("tot2","ToT Spectrum 2",(int)EHistHigh,0,EHistHigh);
	TH1F *timediff_hist = new TH1F("timediff","Timedifference",(int)(2*tLimit),-tLimit,tLimit);

	//printf("%lld\n",tree->GetEntries());
	stic_data_t* event=NULL;
	stic_data_t last_event_1;
	stic_data_t last_event_2;

	tree->SetBranchAddress("br",&event);

	//Coincidence Search:
	//	printf("%d has CCM=%d, CCS=%d, ch=%u\n",i,event->T_CCM,event->T_CCS,event->channel);
	for (int i=0; i<tree->GetEntries(); i++){
	//for (int i=0; i< (int) (tree->GetEntries()/2); i++){

		tree->GetEntry(i);
		//Remove entries with invalid CCS or CCM measurements (time=-1 or 0);
		if(event->time == ((0x1 << 32) - 1) || event->time == 0)
		{
			//printf("Invalid event found \n");
			continue;
			printf("This should never occur\n");
		}
		
		if( event->E_CCM < event->E_CCS || event->E_CCM > event->E_CCS +1 ) continue;
		if( event->T_CCM < event->T_CCS || event->T_CCM > event->T_CCS +1 ) continue;

		if(event->errors.size() != 0) continue; //remove all the events that have error
		
		//Fill the energy spect of CH_1 
		if(event->channel==CH_1)
		{
			if(event->energy !=0 && event->T_fine != 0) 
			{
				last_event_1.time = event->time;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_1.energy = event->energy;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_1.T_fine = event->T_fine;
				last_event_1.T_CCM = event->T_CCM;
				last_event_1.T_CCS = event->T_CCS;
				last_event_1.frame_number = event->frame_number;

				tot_hist_ch1->Fill(event->energy);

			}
		}

		//Fill the energy spect of CH_2 
		if(event->channel==CH_2)
		{
			if(event->energy !=0 && event->T_fine != 0) 
			{
				last_event_2.time = event->time;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_2.energy = event->energy;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_2.T_fine = event->T_fine;
				last_event_2.T_CCM = event->T_CCM;
				last_event_2.T_CCS = event->T_CCS;
				last_event_2.frame_number = event->frame_number;

				tot_hist_ch2->Fill(event->energy);
			}
		}
	}

	//find the 511 keV peak and plot
	double peakEpara[2],tot1Range[2],tot2Range[2];
	double maxY;
	TLine *l1, *l2;
	c1->Divide(2,2);
	c1->cd(1);
	tot_hist_ch1->Draw();
	tot_hist_ch1->GetXaxis()->SetTitle("TOT /1.6ns");
	tot_hist_ch1->GetYaxis()->SetTitle("Counts");
	tot_hist_ch1->GetXaxis()->SetLabelSize(0.05);
	tot_hist_ch1->GetYaxis()->SetLabelSize(0.05);
	tot_hist_ch1->GetXaxis()->SetTitleOffset(1.2);
	tot_hist_ch1->GetYaxis()->SetTitleOffset(1.5);
	tot_hist_ch1->SetTitle("TOT1");
	c1->cd(1)->Update();
	getEPeak(peakEpara,tot_hist_ch1);
	//if(peakEpara[0]>250)	NsigmaE = 1.0; //smaller cut range for higher HV since worse E resolution
	tot1Range[0] = peakEpara[0] - NsigmaE*peakEpara[1];
	tot1Range[1] = peakEpara[0] + NsigmaE*peakEpara[1];
	if(ManuECut == 1 )  //manually set the energy cut range
	{
		tot1Range[0] = CH_1_E_LOW;
		tot1Range[1] = CH_1_E_HIGH;
	}

	//printf("peak range: %f - %f\n",tot1Range[0],tot1Range[1]);

	maxY = gPad->GetUymax(); 
	l1 = new TLine(tot1Range[0],0,tot1Range[0],maxY);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(1);
	l1->Draw();
	l2 = new TLine(tot1Range[1],0,tot1Range[1],maxY);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(1);
	l2->Draw();
	c1->cd(1)->Update();

	c1->cd(2);
	tot_hist_ch2->Draw();
	tot_hist_ch2->GetXaxis()->SetTitle("TOT /1.6ns");
	tot_hist_ch2->GetYaxis()->SetTitle("Counts");
	tot_hist_ch2->GetXaxis()->SetLabelSize(0.05);
	tot_hist_ch2->GetYaxis()->SetLabelSize(0.05);
	tot_hist_ch2->GetXaxis()->SetTitleOffset(1.2);
	tot_hist_ch2->GetYaxis()->SetTitleOffset(1.5);
	tot_hist_ch2->SetTitle("TOT2");
	c1->cd(2)->Update();
	getEPeak(peakEpara,tot_hist_ch2);
	//if(peakEpara[0]>250)	NsigmaE = 1.0; //smaller cut range for higher HV
	tot2Range[0] = peakEpara[0] - NsigmaE*peakEpara[1];
	tot2Range[1] = peakEpara[0] + NsigmaE*peakEpara[1];
	if(ManuECut == 1 )  //manually set the energy cut range
	{
		tot2Range[0] = CH_2_E_LOW;
		tot2Range[1] = CH_2_E_HIGH;
	}

	double nEvnetCH[2],nEventPE[2];
	nEvnetCH[0] = tot_hist_ch1->GetEntries();
	nEvnetCH[1] = tot_hist_ch2->GetEntries();
	nEventPE[0] = tot_hist_ch1->Integral(tot1Range[0],tot1Range[1]);
	nEventPE[1] = tot_hist_ch2->Integral(tot2Range[0],tot2Range[1]);
	printf("Event number in CH %d: %.0f; Events selected for CH %d: %.0f;\nEvent number in CH %d: %.0f; Events selected for CH %d: %.0f;\n",CH_1,nEvnetCH[0],CH_1,nEventPE[0],CH_2,nEvnetCH[1],CH_2,nEventPE[1]);

	maxY = gPad->GetUymax(); 
	l1 = new TLine(tot2Range[0],0,tot2Range[0],maxY);
	l1->SetLineColor(kRed);
	l1->SetLineWidth(1);
	l1->Draw();
	l2 = new TLine(tot2Range[1],0,tot2Range[1],maxY);
	l2->SetLineColor(kRed);
	l2->SetLineWidth(1);
	l2->Draw();
	c1->cd(2)->Update();

	double* data = new double[tree->GetEntries()/2]; //Array containing the coincidences
	unsigned int n;				//number of found coincidences

	double* tot = new double[tree->GetEntries()]; //Array of energy of the coincidences
	int nTotCount = (int) tree->GetEntries()/2;

	n=0;
	for (int i=0; i<tree->GetEntries(); i++){
	//for (int i=0; i<nTotCount; i++){

		tree->GetEntry(i);

		//Remove entries with invalid CCS or CCM measurements (time=-1 or 0);
		if(event->time == ((0x1 << 32) - 1) || event->time == 0)
		{
			continue;
			printf("This should never occur\n");
		}

		if( event->E_CCM < event->E_CCS || event->E_CCM > event->E_CCS +1 ) continue;
		if( event->T_CCM < event->T_CCS || event->T_CCM > event->T_CCS +1 ) continue;
		//if( (event->T_fine >2 && event->T_fine < 6) || (event->T_fine >18 && event->T_fine < 22) ) continue;

		if(event->errors.size() != 0) continue; //remove all the events that have error


		//looking for the 511keV conincidence events
		if(event->channel==CH_1)
		{
			if(event->energy !=0  && event->T_fine != 0 && event->energy < tot1Range[1] && event->energy > tot1Range[0]) 
			{
				last_event_1.time = event->time;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_1.energy = event->energy;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_1.T_fine = event->T_fine;
				last_event_1.T_badhit= event->T_badhit;
				last_event_1.T_CCM = event->T_CCM;
				last_event_1.T_CCS = event->T_CCS;
				last_event_1.E_badhit= event->E_badhit;
				last_event_1.E_CCM = event->E_CCM;
				last_event_1.E_CCS = event->E_CCS;
				last_event_1.channel= event->channel;
				last_event_1.frame_number = event->frame_number;

				if( abs((long int)((int)(event->time) - (int)(last_event_2.time))) < tLimit && event->frame_number == last_event_2.frame_number)
				{
					data[n]=(int)(last_event_1.time) - (int)(last_event_2.time);
					n++;
					timediff_hist->Fill((int)(last_event_1.time) - (int)(last_event_2.time));
					//print out the data of the coincidence event, for debugging
//					if( ((int)(last_event_1.time) - (int)(last_event_2.time) < -10 ) &&((int)(last_event_1.time) - (int)(last_event_2.time) > -20 ) )
//					{
//						printf("event_1\n");
//			 			printf("E_CCM: %u\n",last_event_1.E_CCM);
//			 			printf("E_CCS: %u\n",last_event_1.E_CCS);
//			 			printf("E_BADHIT: %u\n",last_event_1.E_badhit);
//			 			printf("T_FINE: %u\n",last_event_1.T_fine);
//			 			printf("T_CCM: %u\n",last_event_1.T_CCM);
//			 			printf("T_CCS: %u\n",last_event_1.T_CCS);
//			 			printf("T_BADHIT: %u\n",last_event_1.T_badhit);
//			 			printf("CHANNEL: %u\n",last_event_1.channel);
//			 			printf("FRAME: %u\n",last_event_1.frame_number);
//			 			printf("ENERGY: %u\n",last_event_1.energy);
//			 			printf("TIME: %u\n\n",last_event_1.time);
//
//						printf("event_2\n");
//			 			printf("E_CCM: %u\n",last_event_2.E_CCM);
//			 			printf("E_CCS: %u\n",last_event_2.E_CCS);
//			 			printf("E_BADHIT: %u\n",last_event_2.E_badhit);
//			 			printf("T_FINE: %u\n",last_event_2.T_fine);
//			 			printf("T_CCM: %u\n",last_event_2.T_CCM);
//			 			printf("T_CCS: %u\n",last_event_2.T_CCS);
//			 			printf("T_BADHIT: %u\n",last_event_2.T_badhit);
//			 			printf("CHANNEL: %u\n",last_event_2.channel);
//			 			printf("FRAME: %u\n",last_event_2.frame_number);
//			 			printf("ENERGY: %u\n",last_event_2.energy);
//			 			printf("TIME: %u\n\n",last_event_2.time);
//					}
				}
			}
		}

		if(event->channel==CH_2)
		{
			if(event->energy !=0 && event->T_fine != 0 && event->energy < tot2Range[1] && event->energy > tot2Range[0]) 
			{
				last_event_2.time = event->time;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_2.energy = event->energy;	//SAVE THE EVENT AS THE LAST 511 keV event 
				last_event_2.T_fine = event->T_fine;
				last_event_2.T_badhit= event->T_badhit;
				last_event_2.T_CCM = event->T_CCM;
				last_event_2.T_CCS = event->T_CCS;
				last_event_2.E_badhit= event->E_badhit;
				last_event_2.E_CCM = event->E_CCM;
				last_event_2.E_CCS = event->E_CCS;
				last_event_2.channel= event->channel;
				last_event_2.frame_number = event->frame_number;

				if( abs((long int)((int)(event->time) - (int)(last_event_1.time))) < tLimit && event->frame_number == last_event_1.frame_number)
				{
					data[n]=(int)(last_event_1.time) - (int)(last_event_2.time);
					n++;
					timediff_hist->Fill((int)(last_event_1.time) - (int)(last_event_2.time));
					//print out the data of the coincidence event, for debugging
//					if( ((int)(last_event_1.time) - (int)(last_event_2.time) < -10 ) &&((int)(last_event_1.time) - (int)(last_event_2.time) > -20 ) )
//					{
//						printf("event_1\n");
//			 			printf("E_CCM: %u\n",last_event_1.E_CCM);
//			 			printf("E_CCS: %u\n",last_event_1.E_CCS);
//			 			printf("E_BADHIT: %u\n",last_event_1.E_badhit);
//			 			printf("T_FINE: %u\n",last_event_1.T_fine);
//			 			printf("T_CCM: %u\n",last_event_1.T_CCM);
//			 			printf("T_CCS: %u\n",last_event_1.T_CCS);
//			 			printf("T_BADHIT: %u\n",last_event_1.T_badhit);
//			 			printf("CHANNEL: %u\n",last_event_1.channel);
//			 			printf("FRAME: %u\n",last_event_1.frame_number);
//			 			printf("ENERGY: %u\n",last_event_1.energy);
//			 			printf("TIME: %u\n\n",last_event_1.time);
//
//						printf("event_2\n");
//			 			printf("E_CCM: %u\n",last_event_2.E_CCM);
//			 			printf("E_CCS: %u\n",last_event_2.E_CCS);
//			 			printf("E_BADHIT: %u\n",last_event_2.E_badhit);
//			 			printf("T_FINE: %u\n",last_event_2.T_fine);
//			 			printf("T_CCM: %u\n",last_event_2.T_CCM);
//			 			printf("T_CCS: %u\n",last_event_2.T_CCS);
//			 			printf("T_BADHIT: %u\n",last_event_2.T_badhit);
//			 			printf("CHANNEL: %u\n",last_event_2.channel);
//			 			printf("FRAME: %u\n",last_event_2.frame_number);
//			 			printf("ENERGY: %u\n",last_event_2.energy);
//			 			printf("TIME: %u\n\n",last_event_2.time);
//					}
				}
			}
		}
	}
	printf("Coincidence event number: %d. \n\n",n);

	TKDE* kde = new TKDE(n,&data[0],-30,30,"",1.0);
	printf("KDE Sigma: %f\n",kde->GetSigma());

	//TF1 *func;
	func = kde->GetFunction(600,-30,30);

	double HM = func->GetMaximum()/2.0;

	//Now try to find the FWHM
	double fwhm_kde = get_FWHM();



	double mean,sigma,err_gaus,para[16];

	TF1 *f1 = new TF1("f1","gaus");
	//TF1 *f1 = new TF1("f1","gaus(0)+gaus(3)");
	TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",-tLimit,tLimit);
	f1->SetLineColor(2);
	f1->SetLineWidth(2);
	ftemp->SetLineWidth(2);
	ftemp->SetLineStyle(1);
	ftemp->SetLineColor(2);
	c1->cd(3);
	timediff_hist->Draw("pe1");
	timediff_hist->GetXaxis()->SetTitle("time Diff /50.2ps");
	timediff_hist->GetYaxis()->SetTitle("Counts");
	timediff_hist->GetXaxis()->SetLabelSize(0.05);
	timediff_hist->GetYaxis()->SetLabelSize(0.05);
	timediff_hist->GetXaxis()->SetTitleOffset(1.2);
	timediff_hist->GetYaxis()->SetTitleOffset(1.5);
	timediff_hist->SetTitle("CTR");

  	mean=timediff_hist->GetMean();
  	sigma=timediff_hist->GetRMS();
	f1->SetRange(mean-NsigmaT*sigma, mean+NsigmaT*sigma);
	f1->SetParameter(0,700 );
	f1->SetParameter(1, 2) ;
	f1->SetParameter(2,1 );
	f1->SetParameter(3,250 );
	f1->SetParameter(4,15 );
	f1->SetParameter(5,2);
	f1->SetParLimits(0,4 ,1000 );
	f1->SetParLimits(1, -3,9) ;
	f1->SetParLimits(2,0.1 ,10);
	f1->SetParLimits(3,1 ,500 );
	f1->SetParLimits(4,10 ,40 );
	f1->SetParLimits(5,0.1 , 10);
	//timediff_hist->Fit("f1","NQR");
	//f1->GetParameters(&para[0]);
	//f1->SetRange(para[1]-NsigmaT*para[2], para[1]+NsigmaT*para[2]);
	//timediff_hist->Fit("f1","NQR");
	//f1->GetParameters(&para[0]);
	//f1->SetRange(para[1]-NsigmaT*para[2], para[1]+NsigmaT*para[2]);
	timediff_hist->Draw();
	timediff_hist->Fit("f1","");
	f1->GetParameters(&para[0]);
	err_gaus=f1->GetParError(2)*2.35*50.2;
	double fwhm_gaus = 2.35*para[2]*50.2;

	printf("\nFWHM from gaussian fit: %f ps, error: %f ps\n",fwhm_gaus,err_gaus);
	
	ftemp->SetParameters(para[0],para[1],para[2]);
	ftemp->Draw("same");

	double xMax_gaus = ftemp->GetXmax();
	double xMin_gaus = ftemp->GetXmin();
	double maxX_gaus = ftemp->GetMaximumX();
	double max_gaus  = ftemp->GetMaximum();
	double xhm1_gaus = ftemp->GetX(0.5*max_gaus,xMin_gaus,maxX_gaus);
	double xhm2_gaus = ftemp->GetX(0.5*max_gaus,maxX_gaus,xMax_gaus);
	printf("FWHM from gaussian fit (manually search): %f ps\n",50.2*(xhm2_gaus-xhm1_gaus));
	printf("---------------------------------------------------------------------------\n\n");
	TArrow *ar1 = new TArrow(xhm1_gaus-10,0.5*max_gaus,xhm1_gaus-1,0.5*max_gaus,0.01,"|>");
	ar1->SetLineWidth(1);
	ar1->Draw();
	TArrow *ar2 = new TArrow(xhm2_gaus+10,0.5*max_gaus,xhm2_gaus+1,0.5*max_gaus,0.01,"|>");
	ar2->SetLineWidth(1);
	ar2->Draw();

	TPaveText *pt = new TPaveText(xhm2_gaus+3.5,0.5*max_gaus,xMax_gaus,0.65*max_gaus,"bl");
	pt->AddText("FWHM = ");
	pt->AddText(Form("%.2f ps",fwhm_gaus));
	pt->SetBorderSize(0);
	pt->SetFillStyle(0);
	//pt->SetFillColor(0);
	pt->SetTextSize(0.05);
	pt->SetTextAlign(11);
	pt->SetTextFont(12);
	pt->Draw();

	c1->cd(3)->Update();
	//
	//
	c1->cd(4);

	TFile* fout=new TFile("coincidence_result.root","RECREATE"); // create the file to save the plot
//	timediff_hist->Write();
//	tot_hist_ch1->Write();
//	tot_hist_ch2->Write();
	//func->Write("func");
	c1->Write();
	fout->Close();

	fitResultOut=fopen("fitResult.txt","a");
	fprintf(fitResultOut,"%f	%f",fwhm_gaus,err_gaus);
	fprintf(fitResultOut,"\n");
	fclose(fitResultOut);

	
	return 0;
}
