gROOT->Reset();
#include <Riostream.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include "math.h"
#include "string.h"
//#include "TROOT.h"
//#include "Riostream.h"
//#include "TFile.h"
//#include "TChain.h"
//#include "TH1.h"
//#include "TH2.h"
//#include "TStyle.h"
//#include "TCanvas.h"
//#include "TProfile.h"
//#include "TTree.h"
//#include "TNtuple.h"  
//#include "TRandom.h"
//#include "TF1.h"
//#include "TMath.h"
//#include "TGraph.h"
//#include "TGraphErrors.h"
//#include "TPostScript.h"

const int nMax=1e9;
float var[nMax],indx[nMax];
const float nSigma[4] = {1.0,1.5,2.0,2.5};
const float tLimits[2]={-1.0,1.0};
int nNumBins[4]={20,40,80,160};
//int nBins= (tLimits[1]-tLimits[0])*80;

const int nFillMax=1.0e3;

char *directory="/home/huangshan/Documents/measurement/SiPM_Measurements/20130603_HPTDC_SCOPE_compar/scope_data/";
char fileName[1024];
char buf1[1024];

float varUp,varLow,varMax,varMin,var1,varbuf,indx1,indxMin,indxMax,indxStart;
int nEntrs;
//const int nBins=1000;

//float hv[16] = {72.7, 72.5, 72.7, 72.9, 72.9, 72.9, 73.1, 73.1, 73.3, 73.3, 73.5, 73.5, 72.7, 72.5, 72.3, 72.3};
//float hv[17] = {72.3, 72.1, 72.6, 72.5, 72.7, 72.9, 73.1, 73.3, 73.5, 71.9, 71.7, 72.5, 72.7, 72.5, 72.7, 72.5, 72.7};
float hv[3] = {72.3, 72.6, 72.9};

void read_scope_spect()
{
	int i,j,k;
	//for(i=0;i<17;i++)	{	
	cin>>i;
		ifstream indata;
	
		sprintf(fileName,"F620130603-%d-%.1f-00000",i+1,hv[i]);
	//	sprintf(fileName,"F6-21-73.9v-00000");
		sprintf(buf1,"%s%s.txt",directory,fileName);
		indata.open(buf1,ios::in);
		if(!indata) {
	    	cout<<endl<<buf1<<" doesn't exist!!"<<endl<<endl;
	//		exit(-1);
			continue;
	  	}
	  	else {
			int nCount;		
			varMax=0;
			varMin=nMax;
	
			indxMin=nMax;
			indxMax=0;
			nCount=0;
			while(!indata.eof())  {
				indata>>indx1>>var1;
				//cout<<indx1<<"	"<<var1<<endl;
	
				if(indxMin>indx1)	indxMin=indx1;
				if(indxMax<indx1)	indxMax=indx1;
				var[nCount]=var1*1e9;
				indx[nCount]=indx1;
				if(varMax<var[nCount])  varMax=var[nCount];
				if(varMin>var[nCount])  varMin=var[nCount];
				nCount++;
				
			}
			//h1->Write();
			indata.close();
			//file1->Close();
		
				varbuf=varMax-varMin;
				varMax=varMax+0.025*varbuf;
				varMin=varMin-0.025*varbuf;
	
			//TH1D *spec = new TH1D("SiPM","Spectrum",nBins,varMin,varMax);
			//TH1D *spec = new TH1D("SiPM","Spectrum",nBins,tLimits[0],tLimits[1]);
	//		TGraph *spec = new TGraph(nCount,indx,var);
	//		for(j=0;j<nCount;j++)	{
	//			spec->Fill(var[j]);
	//		}
			
			gStyle->SetOptFit(1111);
			gStyle->SetOptStat(1111);
			
			double mean,sigma,para[16],err;
			
			//////////////output fit result to txt
			ofstream outResult;
			outResult.open(Form("%sfitResult.txt",directory),ios::app);
			cout<<fileName<<endl;
			outResult<<"########    file: "<<fileName<<".txt    #######"<<"\n";
			/////////////////////////////////

			int nBins,nFill;
			TH1D *spec;
			TCanvas *c1 = new TCanvas(fileName,fileName,10,10,700,600);
			for(k=0;k<3;k++)	{
				nBins= (tLimits[1]-tLimits[0])*1.0*nNumBins[k];
				spec = new TH1D("SiPM","Spectrum",nBins,tLimits[0],tLimits[1]);
				nFill=0;
				for(j=0;j<nCount;j++)	{
					if(var[j]>=tLimits[0] && var[j]<=tLimits[1])	{
						spec->Fill(var[j]);
						nFill++;
					}
					if(nFill==nFillMax)	break;
				}
				c1->Divide(2,2);
			for(j=0;j<4;j++)	{
				c1->cd(j+1);
				spec->Draw();
			//	spec->Draw("AB");
			//	spec->SetFillColor(38);
			//	spec->GetXaxis()->SetLimits(tLimits[0],tLimits[1]);
		
				spec->GetXaxis()->SetTitle("t [ns]");
		//		spec->GetXaxis()->SetTitle("Maximum of SiPM signals [mV]");
				spec->GetXaxis()->CenterTitle();
				spec->GetYaxis()->SetTitle("Counts");
				spec->GetYaxis()->CenterTitle();
				spec->SetTitle("CTR_2_STiC");
				
				spec->Fit("gaus","NQ");
				gaus->GetParameters(&para[0]);
				gaus->SetRange(para[1]-nSigma[j]*para[2], para[1]+nSigma[j]*para[2]);
				spec->Fit("gaus","NQR");
				gaus->GetParameters(&para[0]);
				gaus->SetRange(para[1]-nSigma[j]*para[2], para[1]+nSigma[j]*para[2]);
				spec->Fit("gaus","NQR");
				gaus->GetParameters(&para[0]);
				gaus->SetRange(para[1]-nSigma[j]*para[2], para[1]+nSigma[j]*para[2]);
				spec->Fit("gaus","QR");
				//spec->Fit("gaus");
				gaus->GetParameters(&para[0]);
				err=gaus->GetParError(2);

				cout<<nNumBins[k]<<"	"<<nSigma[j]<<"	FWHM:	"<<para[2]*2.35<<"\t"<<err*2.35<<"\n";
			
				TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",tLimits[0],tLimits[1]);
				//TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",varMin,varMax);
				ftemp->SetParameters(para[0],para[1],para[2]);
				ftemp->SetLineWidth(3);
				ftemp->SetLineStyle(2);
				ftemp->SetLineColor(2);
				ftemp->Draw("same");
				c1->cd(j+1)->Update();
				
				outResult<<nSigma[j]<<'	'<<tLimits[0]<<'	'<<tLimits[1]<<"	"<<nBins<<"\n";
				outResult<<"Sigma[ps]	Err_of_Sigma[ps]	CTR_FWHM	Err_of_FWHM"<<"\n";
			  	outResult<<para[2]<<"\t"<<err<<"\t"<<para[2]*2.35<<"\t"<<err*2.35<<"\n\n";
			}
			wait();
			c1->Clear();
			}
			c1->Close();
	
	//		///////////save histogram to .root file
	//		TFile *fHist = new TFile(Form("%s%s-hist.root",directory,fileName),"recreate");
	//		c1->Write();
	//		fHist->ls();
	//		fHist->Close();
	//		//////////
		
			////////close output file
			outResult<<"\n\n";
			outResult.close();
			/////////////////////////
		
	//		wait();
	//		c1->Close();
		}
	//}
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

