#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>

#include "math.h"
#include "string.h"

#include "TTree.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TStyle.h"


char* directory="/home/huangshan/Dropbox/measurement/STiC3_measurements/20140517_thr_scan/scan_result/";

void ReadPlot(char* inFile, char* outFile, char* outTxtFile)
{
	int i;
	int nPoint=0;
	double inbuf[7];
	char fileName[1024];
	TGraph* gr[3];
	for(i=0;i<3;i++) gr[i]=new TGraph();
	sprintf(fileName,"%s%s.txt",directory,inFile);
	ifstream indata;
	indata.open(fileName,ios::in);
	if(!indata) {
		printf("%s doesn't existed!!\n",inFile);
		return;
	}
	else {
		while(!indata.eof())  {
			indata>>inbuf[0]>>inbuf[1]>>inbuf[2]>>inbuf[3]>>inbuf[4]>>inbuf[5]>>inbuf[6];
			gr[0]->SetPoint(nPoint,inbuf[0]*1e3*33.0,inbuf[3]/inbuf[1]);
			gr[1]->SetPoint(nPoint,inbuf[0]*1e3*33.0,inbuf[2]*1e12);
			if(inbuf[3]>0) gr[2]->SetPoint(nPoint,inbuf[0]*1e3*33.0,inbuf[6]/inbuf[3]);
			nPoint++;
		}
	}

//	TF1 *erf = new TF1("erf","([0]*TMath::Erf((x-[1])/[2]))+[3]");
//	erf->SetRange(0,1000);
//	erf->SetParameter(0,0.5);
//	erf->SetParameter(1,100);
//	erf->SetParameter(2,1.0);
//	erf->SetParameter(3,0.5);
//	erf->SetParLimits(0,0.47,0.53);
//	erf->SetParLimits(1,0,500);
//	erf->SetParLimits(2,0.01,100);
//	erf->SetParLimits(3,0.47,0.53);
//
	TF1 *erf = new TF1("erf","(0.5*TMath::Erf((x-[0])/[1]))+0.5");
	erf->SetRange(0,1000);
	erf->SetParameter(0,200);
	erf->SetParameter(1,2.0);
	erf->SetParLimits(0,0,500);
	erf->SetParLimits(1,0.1,100);
	double par[4];

	//gROOT->Reset();
	gStyle->SetMarkerStyle(24);
	gStyle->SetMarkerSize(1);
	gStyle->SetMarkerColor(kBlue);

	TCanvas *c1 = new TCanvas(inFile,"c1",800,600);
	c1->Divide(2,1);
	c1->cd(1);
	gr[0]->GetXaxis()->SetTitle("Input Charge [fC]");
	gr[0]->GetYaxis()->SetTitle("Efficiency");
	gr[0]->GetYaxis()->SetRangeUser(-0.05,1.1);
	gr[0]->SetTitle(Form("%s_efficiency",inFile));
	//gr[0]->Draw("AP");
	erf->SetLineColor(kRed);
	for(i=0;i<5;i++) gr[0]->Fit("erf","NQMR+");
	gr[0]->Fit("erf","MR");
	gr[0]->Draw("AP");
	erf->GetParameters(par);
	

	c1->cd(2);
	gr[1]->GetXaxis()->SetTitle("Input Charge [fC]");
	gr[1]->GetYaxis()->SetTitle("Time Jitter [ps]");
	gr[1]->SetTitle(Form("%s_time_jitter",inFile));
	gr[1]->Draw("AP");

//	c1->cd(3);
//	gr[2]->SetMarkerStyle(24);
//	gr[2]->SetMarkerSize(1);
//	gr[2]->SetMarkerColor(kRed);
//	gr[2]->GetXaxis()->SetTitle("Input Charge [fC]");
//	gr[2]->GetYaxis()->SetTitle("nTot/nT");
//	gr[2]->Draw("AP");

	TFile *fout = new TFile(outFile,"UPDATE");
	gr[0]->Write(Form("%s_efficiency",inFile));
	gr[1]->Write(Form("%s_time_jitter",inFile));
//	gr[2]->Write();
	c1->Write();
	fout->Close();

	FILE * foutTxt = fopen(outTxtFile,"a");
	fprintf(foutTxt,"%f\t%f\n",par[0],par[1]);
	fclose(foutTxt);
}

// main function
void thr_scan()
{
	int i;
	char inFile[1024];
	char outFile[1024];
	char outTxtFile[1024];

	sprintf(outFile,"%sresult.root",directory);
	sprintf(outTxtFile,"%sresult.txt",directory);
	
	TFile *fout = new TFile(outFile,"RECREATE");
	fout->Close();
	FILE * pFile;
	pFile = fopen(outTxtFile,"w"); 	
	fclose(pFile);

	for(i=0;i<54;i+=2)
	{
		pFile = fopen(outTxtFile,"a"); 	
		fprintf(pFile,"%d\t",i);
		fclose(pFile);

		sprintf(inFile,"thr_%d_scanResult",i);
		ReadPlot(inFile,outFile,outTxtFile);
	}
}
