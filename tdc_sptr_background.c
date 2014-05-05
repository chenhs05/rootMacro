//root macro for sptr analysis with HPTDC
//create by H.Chen
//
//
//#include <TROOT.h>
//gROOT->Reset();

#ifndef __CINT__
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
#endif // __CINT__

#include "tdc_ch_values.cpp"
#include "include/myCali.h"

#include "TTree.h"
#include "TGraph.h"
#include "TFile.h"
#include "TLine.h"

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/201401_sptr/20140212_thScan_5/";
char fileName[1024];
const int chStic = 4;
const int edgeStic[2] = {1,0};

const float binTDC = 25.0; ///24.414; //ps
const int totHyst = 0; //
const float maxDelay = 4096;	//*25ps

const int nMax = 1e7; 
int nEvents;

const float minTot = -50;
const float maxTot = 500;
const float binningAllTot =0.5;
int nBinsAllTot = int((maxTot-minTot)/binningAllTot);

double minPeakTot[2] = {0,0};
double maxPeakTot[2] = {400,0};
float binningPeakTot = 0.5;
int nBinsPeakTot = int((maxPeakTot[0]-minPeakTot[0])/binningPeakTot);

int NLCorr=0;
/*********** create the NL correction mapping for HPTDC ************************************
 ********************************************************************************************/
void CreatTdcCorrection(char *file,tdc_ch_values *chan,int channel,int edge)
{
	int chIn,tsIn,edgeIn;
	printf("####	Create TDC NL correction for file:\n%s\n",file);
	ifstream is;
	is.open(file,ios::in);
	while(!is.eof())
	{
		is>>chIn>>tsIn>>edgeIn;
		if(chIn == channel && edgeIn == edge)	chan->fill(tsIn);
	}
	chan->update_histos();
	is.close();
}

/************  Read .dat file and convert to .root file***************************************
*************  correction with class tdc_ch_value     ****************************************
*********************************************************************************************/
void ReadFile()
{
	long int chIn,tsIn,edgeIn;
	long int tsBuf,tsSticBuf,totSticBuf;
	double tot;
	char buf1[1024];
	char file1[1024];

	TH1D *histTOT = new TH1D("TOT","TOT",nBinsAllTot,minTot,maxTot);

	sprintf(file1,"%s%s.dat",directory,fileName);

	tdc_ch_values *chanStic;	
	if(NLCorr == 1)	{
		printf("\n####	Start to read file with TDC NL correction!	####\n");
		chanStic = new tdc_ch_values(1,1024);
		CreatTdcCorrection(file1,chanStic,chStic,edgeStic[0]); 
		chanStic->update_histos();
		printf("----> finish create TDC NL correction!	<---- \n");
	}

	cout<<"####	Read File	"<<file1<<"	####"<<"\n";
	ifstream is;
	is.open(file1,ios::in);

	is.clear();
	is.seekg(0,ios::beg);
	nEvents = 0;

	tsBuf=0;
	tsSticBuf  = 0;
	int i = 0;
	while(!is.eof())
	{
		is>>chIn>>tsIn>>edgeIn;
		if(tsIn<tsBuf-1e5)	{
			tsSticBuf  = 0;
		}	
		tsBuf=tsIn;
		if(i%100000 == 0)	cout<<"!!!!    Reading Data	"<<i<<"!!!!	nEvents = "<<nEvents<<"\n";
		i++;

		if(chIn == chStic)	{
			if(edgeIn == edgeStic[0])	{
				if(NLCorr == 1) tsSticBuf = chanStic->get_dith_value(tsIn);
				else tsSticBuf = tsIn;
			}
			else	{
				if(tsSticBuf != 0)	{
					if(NLCorr == 1 ) totSticBuf = chanStic->get_dith_value(tsIn)-tsSticBuf;
					else totSticBuf = tsIn-tsSticBuf;
					tot = totSticBuf * binTDC * 0.001; //ns 
					histTOT->Fill(tot);
					tsSticBuf =0;
				}
			}
		}
		else	{
			printf("channel = %d, not STiC channel.\n",chIn);
		}

	}
	is.close();
	
	cout<<"####	Complete reading file :"<<file1<<" !	####"<<"\n";
	cout<<"####	Start to plot histogram!	####"<<""<<"\n";

	TCanvas *c1 = new TCanvas("status1","status1",10,10,800,600);
	histTOT->Draw();
	histTOT->GetXaxis()->SetTitle("TOT [ns]");
	histTOT->GetXaxis()->CenterTitle();
	histTOT->GetYaxis()->SetTitle("Counts");
	histTOT->GetYaxis()->CenterTitle();
	c1->Update();
	double maxY;
	maxY = gPad->GetUymax(); 
//	cout<<"max of Y axis = "<<maxY<<"\n";
	TLine *l1 = new TLine(minPeakTot[0],0,minPeakTot[0],maxY);
	l1->SetLineColor(kRed);
	l1->Draw();
	TLine *l2 = new TLine(maxPeakTot[0],0,maxPeakTot[0],maxY);
	l2->SetLineColor(kRed);
	l2->Draw();
	c1->Update();

	sprintf(buf1,"%s%s_RunStatus.root",directory,fileName);
	TFile *hf = new TFile(buf1,"recreate");
	c1->Write();

	hf->Close();
	c1->Close();

	cout<<"####	Complete ploting spectrums!	####"<<"\n"<<"\n";
	cout<<"####	Function ReadFile complete!	####"<<"\n";
}


void tdc_sptr_background()
{
	char buf[1024];
	int i,j,k;
	for(i=0;i<5;i++)	{
		for(j=0;j<5;j++)	{
			for(k=0;k<5;k++)	{
				sprintf(fileName,"sptr_%d_%d_%d",14+i*3,2+j*1,1+k);
				sprintf(buf,"%s%s.dat",directory,fileName);	
			
				ifstream isTest;
				isTest.open(buf,ios::in);
				if(!isTest) {
					cout<<"\n"<<buf<<" doesn't exist!!"<<"\n"<<"\n";
					continue;
				}
				else {
					isTest.close();
					cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
					ReadFile();
				}
			}
		}
	}
}
