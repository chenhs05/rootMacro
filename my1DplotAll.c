gROOT->Reset();
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
//#include"myWait.c"

//const char* file[14] = {"1120_1","1120_2","1120_3","1118_2","1119_1","1126_1","1122_2","1122_3","1121_3","1121_4","1121_5","1122_1"};
//const char* capa[14] = {"0.1nF","0.1nF","0.1nF","1nF","1nF","10nF","100nF","100nF","cml=0","cml=15","cml=45","cml=60"};
const char* file[14] = {"1120_1","1120_2","1120_3","1118_2","1119_1","1126_1","1122_2","1122_3"};
const char* capa[14] = {"0.1nF","0.1nF","0.1nF","1nF","1nF","10nF","100nF","100nF"};

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/201310_ctr/capa/"; 
char fileName[1024]; 
void my1DplotAll()	{
	int n,buf,nFile;
	double ctr[1024],err[1024],bias1[1024],bias2[1024],errX[1024];
	TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,600);
	TGraph *g;
	TLegend *leg = new TLegend(0.8,0.4,0.9,0.9);
	int nLoop=0;

// set color patten	

	int MyPalette[100];
	double red[]  = { 0.0, 1.0, 0.0};
	double green[]= { 1.0, 0.0, 0.0};
	double blue[] = { 0.0, 0.0, 1.0};
	double stop[] = { 0.0, .50, 1.0};
	int FI = TColor::CreateGradientColorTable(3, stop, red, green, blue, 100);
//	for (int i=0;i<100;i++) MyPalette[i] = FI+i;
//	gStyle->SetPalette(100, MyPalette);

//

	for(nFile=0;nFile<8;nFile++)
	{
		sprintf(fileName,"%s%s.txt",directory,file[nFile]);
		printf("\n%s\n",fileName);
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
	
		g = new TGraph(n-1, bias1, ctr);
	//	g->Draw("AL");
	//	for(int i=0;i<1024;i++)	errX[i]=0;
	//	TGraphErrors *g = new TGraphErrors(n-1,bias1,ctr,errX,err);
		if(nLoop == 0)	{
			g->Draw("APL");
			g->SetTitle("CTR_scan");
			g->GetXaxis()->SetTitle("HV /V");
			//g->GetXaxis()->SetTitle("DAC");
			g->GetYaxis()->SetTitle("CTR FWHM /ps");
			g->GetXaxis()->CenterTitle();
			g->GetYaxis()->CenterTitle();
			g->GetXaxis()->SetLimits(72.2,74.6);
			g->GetYaxis()->SetRangeUser(200.0,350.0);
		}
		else g->Draw("PL");
		//g->Fit("pol2");
		g->SetMarkerSize(0.5);
		g->SetMarkerStyle(4);
		g->SetMarkerColor(FI+10*nLoop);
		g->SetLineColor(FI+10*nLoop);
		//pol2->SetLineWidth(0.5);
		leg->AddEntry(g,capa[nFile],"l");

		c1->Update();

		nLoop++;
	}
		leg->Draw();
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
