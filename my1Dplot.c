gROOT->Reset();
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
//#include"myWait.c"

char *fileName = "/home/huangshan/Documents/measurement/SiPM_Measurements/201310_ctr/20131127_hvScan/caliResult.txt"; 
void my1Dplot()	{
	int n,buf;
	double ctr[1024],err[1024],bias1[1024],bias2[1024],errX[1024];
	ifstream infile;
	infile.open(fileName,ios::in);

	if(!infile) {
		cout<<"\n"<<fileName<<" doesn't exist!!"<<"\n"<<"\n";
		return;
	}
	else {
		n=0;
		while(!infile.eof())	{
			infile>>bias1[n]>>ctr[n]>>err[n];
			//if(bias1[n]<60)	continue;
			if(ctr[n]>800 || err[n]>150) continue;
			n++;
		}
	}
	infile.close();

	TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,600);

//	TGraph *g = new TGraph(n-1, bias1, ctr);
//	g->Draw("AL");
	for(int i=0;i<1024;i++)	errX[i]=0;
	TGraphErrors *g = new TGraphErrors(n-1,bias1,ctr,errX,errX);
	g->SetTitle("CTR_scan");
	g->GetXaxis()->SetTitle("HV /V");
	//g->GetXaxis()->SetTitle("DAC");
	g->GetYaxis()->SetTitle("CTR FWHM /ps");
	g->GetXaxis()->CenterTitle();
	g->GetYaxis()->CenterTitle();
	g->Draw("APL");
	//g->Fit("pol2");
	g->SetMarkerSize(0.5);
	g->SetMarkerStyle(4);
	//pol2->SetLineWidth(0.5);

	c1->Update();

	wait();
	c1->Close();
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
