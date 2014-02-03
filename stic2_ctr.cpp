gROOT->Reset();
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "math.h"
#include "string.h"
#include "include/myCali.h"
#include "include/EventType.h"

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/201310_CTR/";
char *fileName = "hvScan_73.50";

const unsigned short ch1 = 5;
const unsigned short ch2 = 10;

const long int tWin = 10000;
const float tLimit[2]={-1.0,2.0};
const int tBins = (tLimit[1]-tLimit[0])*20;
const float nSigma = 1.0;
const float nTotSigma = 1.0;
const float resTDC = 0;

const bool modeSelect = 0;
const unsigned int nSelect = 1000;

void ReadFile()
{
	stic_data_t *current_event, *last_event;
	current_event = new stic_data_t;
	last_event = new stic_data_t;

	unsigned short frameNum;

	TFile *fIn = new TFile(Form("%s%s.root",directory,fileName));
	TTree *tree = (TTree*)fIn->Get("dump");
	tree->SetBranchAddress("packet_number",&(current_event->packet_number));
	tree->SetBranchAddress("frame_number",&(current_event->frame_number));
	tree->SetBranchAddress("channel",&(current_event->channel));
	tree->SetBranchAddress("T_CCM",&(current_event->T_CCM));
	tree->SetBranchAddress("T_CCS",&(current_event->T_CCS));
	tree->SetBranchAddress("T_badhit",&(current_event->T_badhit));
	tree->SetBranchAddress("T_fine",&(current_event->T_fine));
	tree->SetBranchAddress("E_CCM",&(current_event->E_CCM));
	tree->SetBranchAddress("E_CCS",&(current_event->E_CCS));
	tree->SetBranchAddress("E_badhit",&(current_event->E_badhit));
	tree->SetBranchAddress("time",&(current_event->time));
	tree->SetBranchAddress("energy",&(current_event->energy));
	tree->SetBranchAddress("errors",&(current_event->errors));
	int nEvents = tree->GetEntries();

	unsigned int tot1,tot2;
	double tDiff;
	TFile *fOut = new TFile(Form("%s%s_coin.root",directory,fileName),"RECREATE","Root file for coincidance data");
	TTree *h1 = new TTree("h1","Data for CTR");
	h1->Branch("tDiff",&tDiff,"tDiff/D");
	h1->Branch("tot1",&tot1,"tot1/i");
	h1->Branch("tot2",&tot2,"tot2/i");
	for(int i=0;i<nEvents;i++)
	{
		tree->GetEntry(i);
		if(i==0) *last_event=*current_event;
		if(current_event->errors.size()==0 )
		{
			cout<<current_event->frame_number<<'\t'<<last_event->frame_number<<endl;
			if(last_event->errors.size()==0 && current_event->frame_number == last_event->frame_number && abs( (long int)((int)(current_event->time) - (int)(last_event[fmonitors->at(i_monitor).fChannel2].time)) ) < tWin)
			{
				if(current_event->channel == ch1 && last_event->channel == ch2)
				{
					tot1=current_event->energy;
					tot2=last_event->energy;
					tDiff = 0.05*((int)current_event->time - (int)last_event->time);
					h1->Fill();
				}
				else if(current_event->channel == ch2 && last_event->channel == ch1)
				{
					tot1=last_event->energy;
					tot2=current_event->energy;
					tDiff = 0.05*(int)last_event->time - (int)current_event->time;
					h1->Fill();
				}
			printf("Coming here !!\n");
			}
			*last_event=*current_event;
		}
	}
	h1->Write();
	fOut->Close();
}

/************  Read .root file and plot the spectrum *****************************************
*********************************************************************************************/
void PlotSpectrum()
{
	double mean,sigma,para[16],rang[2],err;
	double tot1Range[2],tot2Range[2];
	gStyle->SetOptStat(1111);
	gStyle->SetOptFit(111);
	TFile *f;

	f = new TFile(Form("%s%s_coin.root",directory,fileName));
	TTree *T = (TTree*)gDirectory->Get("h1");
	T->Print();	

	//gROOT->ProcessLine("gROOT->SetBatch()");
	TCanvas *c1 = new TCanvas("CTR measurement","CTR measurement results",40,20,1000,1000);
	c1->Divide(2,2);
	c1->cd(1);
	T->Draw(Form("tot1>>TOT1(%d,%f,%f)",400,0.0,400.0));
	TSpectrum *s = new TSpectrum(2);
	int nfound = s->Search(TOT1,2,"",0.05);
	float *xpeak = s->GetPositionX();
	TOT1->Fit("gaus","NQ");
	gaus->SetRange(xpeak[0]-50,xpeak[0]+50);
	TOT1->Fit("gaus","RQ");
	gaus->GetParameters(&para[0]);
	tot1Range[0]=para[1]-nTotSigma*para[2];
	tot1Range[1]=para[1]+nTotSigma*para[2];
	TCut totcut1 = Form("tot1>%.1f && tot1<%.1f ",tot1Range[0],tot1Range[1]);
	double maxY = gPad->GetUymax(); 
	TLine *l1 = new TLine(tot1Range[0],0,tot1Range[0],maxY);
	l1->SetLineColor(kRed);
	l1->Draw();
	TLine *l2 = new TLine(tot1Range[1],0,tot1Range[1],maxY);
	l2->SetLineColor(kRed);
	l2->Draw();
	TOT1->GetXaxis()->SetTitle("TOT (CC)");
	TOT1->GetYaxis()->SetTitle("Counts");
	TOT1->GetXaxis()->SetLabelSize(0.05);
	TOT1->GetYaxis()->SetLabelSize(0.05);
	TOT1->SetTitle("TOT1");
	c1->cd(1)->Update();

	c1->cd(2);
	T->Draw(Form("tot2>>TOT2(%d,%f,%f)",400,0.0,400.0));
	nfound = s->Search(TOT2,2,"",0.05);
	xpeak = s->GetPositionX();
	TOT2->Fit("gaus","NQ");
	gaus->SetRange(xpeak[0]-50,xpeak[0]+50);
	TOT2->Fit("gaus","RQ");
	gaus->GetParameters(&para[0]);
	tot2Range[0]=para[1]-nTotSigma*para[2];
	tot2Range[1]=para[1]+nTotSigma*para[2];
	TCut totcut2 = Form("tot2>%.1f && tot2<%.1f ",tot2Range[0],tot2Range[1]);
	maxY = gPad->GetUymax(); 
	TLine *l1 = new TLine(tot2Range[0],0,tot2Range[0],maxY);
	l1->SetLineColor(kRed);
	l1->Draw();
	TLine *l2 = new TLine(tot2Range[1],0,tot2Range[1],maxY);
	l2->SetLineColor(kRed);
	l2->Draw();
	TOT2->GetXaxis()->SetTitle("TOT (CC)");
	TOT2->GetYaxis()->SetTitle("Counts");
	TOT2->GetXaxis()->SetLabelSize(0.05);
	TOT2->GetYaxis()->SetLabelSize(0.05);
	TOT2->SetTitle("TOT2");
	c1->cd(2)->Update();

	c1->cd(3);
	if(modeSelect ==1)	{
		TCut totcut = totcut1 && totcut2 && Form("Entry$ < %d", nSelect);
	}
	else	{
		TCut totcut = totcut1 && totcut2;	
	}
	T->Draw(Form("tDiff>>CTR1(%d,%f,%f)",tBins,tLimit[0],tLimit[1]),totcut);
	cout<<"BinWidth"<<CTR1->GetBinWidth(1)<<endl;
	CTR1->Fit("gaus","NQ");
  	mean=CTR1->GetMean();
  	sigma=CTR1->GetRMS();
	gaus->SetRange(mean-nSigma*sigma,mean+nSigma*sigma);
	gaus->GetRange(rang[0],rang[1]);
	//cout<<"range:	"<<rang[0]<<"	"<<rang[1]<<endl;
	CTR1->Fit("gaus","NQR");
	gaus->GetParameters(&para[0]);
	gaus->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
	gaus->GetRange(rang[0],rang[1]);
	cout<<"range:	"<<rang[0]<<"	"<<rang[1]<<endl;
	CTR1->Fit("gaus","NQR");
	gaus->GetParameters(&para[0]);
	gaus->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
	gaus->GetRange(rang[0],rang[1]);
	cout<<"range:	"<<rang[0]<<"	"<<rang[1]<<endl;
	CTR1->Fit("gaus","R");
	CTR1->GetXaxis()->SetTitle("t Diff(ns)");
	CTR1->GetYaxis()->SetTitle("Counts");
	CTR1->GetXaxis()->SetLabelSize(0.05);
	CTR1->GetYaxis()->SetLabelSize(0.05);
	CTR1->SetTitle("CTR");

	gaus->GetParameters(&para[0]);
	err=gaus->GetParError(2);
	
	TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",tLimit[0],tLimit[1]);
	ftemp->SetParameters(para[0],para[1],para[2]);
	ftemp->SetLineWidth(3);
	ftemp->SetLineStyle(2);
	ftemp->SetLineColor(2);
	ftemp->Draw("same");
	c1->cd(3)->Update();


//	///////////save histogram to .root file
//	TFile *fHist = new TFile(Form("%s%s-hist.root",directory,fileName),"recreate");
//	c1->Write();
//	fHist->ls();
//	fHist->Close();
//	//////////
//
//	//////////////output fit result to txt
//	double ctrRes,ctrErr;
//	ctrRes=sqrt(para[2]*para[2]-resTDC*resTDC)*2.35*1000.0;
//	ctrErr=err*2.35*1000.0;
//
//	FILE *ofile;
//	ofile = fopen(Form("%sfitResult.txt",directory),"a");
//	fprintf(ofile,"########	%s	#######\n",fileName);
//	fprintf(ofile,"nSigma	NLCorr	tLimit[0]	tLimit[1]	tBins	tot1Range[0]	tot1Range[1]	tot2Range[0]	tot2Range[1]\n");
//	fprintf(ofile,"%.2f	%.2f	%.2f	%d	%.2f	%.2f	%.2f	%.2f\n",nSigma,tLimit[0],tLimit[1],tBins,tot1Range[0],tot1Range[1],tot2Range[0],tot2Range[1]);
//	fprintf(ofile,"Sigma[ps]	Err_of_Sigma[ps]\n");
//	fprintf(ofile,"%f	%f\n\n",ctrRes,ctrErr);
//	fclose(ofile);
//	/////////////////////////
	
	//wait to close the file
//	wait();
	f->Close();	
	c1->Close();
	//gROOT->ProcessLine("gROOT->SetBatch(kFALSE)");
	
}


/************ wait funtion !!  *************************************************************** 
*********************************************************************************************/
void wait(){

	TTimer *timer = new TTimer("gSystem->ProcessEvents();", 50, kFALSE);
	char *input;
	Bool_t done = kFALSE;
	do{
		timer->TurnOn();
		timer->Reset();
		input=Getline("Type <return> to continue : ");
		timer->TurnOff();
		if(input){			done =kTRUE;		}	
	}while(!done);
}

/************ Main funtion !!  *************************************************************** 
*********************************************************************************************/
void stic2_ctr()
{
	char buf[1024];
	sprintf(buf,"%s%s.root",directory,fileName);	
	TFile *f1 = new TFile(buf);
	if(f1->IsZombie())	{
		cout<<"Error openning file!!"<<"\n";
		return;
	}
	else {
		cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n\n";
		f1->Close();
		ReadFile(); 
		PlotSpectrum();
		//DoCali();
	}
}
