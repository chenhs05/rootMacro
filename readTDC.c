gROOT->Reset();
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>

#include<stdio.h>
#include<unistd.h>
#include<stdlib.h>
#include<TSpectrum.h>
#include<TPad.h>
#include<time.h>

#include "math.h"
#include "string.h"
#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TTree.h"
#include "TNtuple.h"  
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TPostScript.h"
#include "TFrame.h"
//#include "tdc_ch_values.h"
#include "tdc_ch_values.cpp"

char *directory = "/home/huangshan/Documents/measurement/SiPM_Measurements/20121203TDC_Scope";
char *fileName = "20121203";

const float binTDC = 24.414; //ps
const int totHyst = 409; //10ns
const float maxDelay = 2048;	//50ns
const int chLaser = 1;
const int chStic = 4;
const int edgeLaser[2] = {0,1};
const int edgeStic[2] = {1,0};

const int nMax = 1e7; 
float t[nMax],tot[nMax],totLaser[nMax];
int nEvents;
int nSelect = 13000;
bool modeSelect = 1;

const float minT = 23.50;
const float maxT = 29.50;
const float binningAllT = 0.05;
int nBinsAllT = int((maxT-minT)/binningAllT);
const float minTot = 20;
const float maxTot = 250;
const float binningAllTot =0.65;
int nBinsAllTot = int((maxTot-minTot)/binningAllTot);

const float shift_t = -27;
const int nCor = 5;
const float nSigma = 1.5; 
const float halfPeakWidth = 1.500;

const float binningPeakT = 0.025;
int nBinsPeakT = int(2*halfPeakWidth/binningPeakT);
const float minPeakTot[2] = {52,0};
const float maxPeakTot[2] = {65,0};
const float binningPeakTot = 0.65;
int nBinsPeakTot = int((maxPeakTot[0]-minPeakTot[0])/binningPeakTot);
//const int nBinsPeakTot = 20;
const float fitMinTot[2] = {minPeakTot[0],0};
const float fitMaxTot[2] = {maxPeakTot[0],0};

int nCharge = 1;

//TF1 *fitf = new TF1("fitf","[0]+[1]/sqrt(x)+[2]/x+[3]*x");  //½ÃÕýº¯Êý
//TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6");
TF1 *fitf = new TF1("fitf","[0]+[1]*x+[2]*x**2+[3]*x**3+[4]*x**4+[5]*x**5+[6]*x**6+[7]*x**7+[8]*x**8+[9]*x**9+[10]*x**10");
/************************ Do the T-Tot correction ****************************************************************************
//½ÃÕýº¯Êý by Huangshan Chen, ¸ÄÐ´/·â×°×Ô Íõ¾°²¨µÄ wjb_test.c
//ttm£ºÊ±¼äÊý¾ÝµÄÊ×µØÖ·£¬Ò»Î¬Êý×éµÄ»°ÌîÊý×éÃû£»
//tot£ºµçºÉÊý¾ÝµÄÊ×µØÖ·£¬¶þÎ¬Êý¾ÝµÄ»°Ìîtot[0];
//nCount: ÓÐÐ§eventµÄ¸öÊý£» 
//nCharge£ºÐèÒª½ÃÕýµÄµçºÉ¸öÊý£»
//nCor£º½ÃÕý´ÎÊý£»
//psName£ºÊä³öpsÎÄ¼þµÄÎÄ¼þÃû£¬Õâ¸öº¯Êý»áÔÚdirectory Ä¿Â¼ÏÂÉú³ÉpsName_Calibration.psÎÄ¼þ£»
//tLimit£ºfloat¸ñÊ½£¬Ê±¼ä·Ö²¼Í¼µÄ·¶Î§ÊÇ [-tLimit,tLimit]£»
//tBins: Ê±¼ä·Ö²¼µÄbinÊý£»
//totMin: ·ù¶È»òÕßtime over threshold·Ö²¼µÄ×îÐ¡Öµ£¬ÓÃÓÚ¶¨T-QÍ¼µÄxÖáµÄ·¶Î§£¬¶ÔÓ¦²»Í¬µÄµçºÉµÄ·ù¶È·¶Î§£¬ÌîÊý×éµÄÊ×µØÖ·£»
//totMin: ·ù¶È»òÕßtime over threshold·Ö²¼µÄ×î´óÖµ£¬ÓÃÓÚ¶¨T-QÍ¼µÄxÖáµÄ·¶Î§£¬ÌîÊý×éµÄÊ×µØÖ·£»
//totBins: T-Q¹ØÏµÍ¼xÖáµÄbinÊý£»
//fitMin£º ½ÃÕýµÄ·¶Î§µÄ×îÐ¡Öµ£¬ÌîÊý×éµÄÊ×µØÖ·£»
//fitMax£º ½ÃÕýµÄ·¶Î§µÄ×î´óÖµ£¬ÌîÊý×éµÄÊ×µØÖ·£»
//fitf£º ½ÃÕýº¯Êý£»
//Nsigma£º ¶ÔÊ±¼ä·Ö²¼½øÐÐGaussion fitµÄ·¶Î§ÊÇÈ¡Æ½¾ùÖµ×óÓÒNsigma¸ösigma
//nTotMax£ºttm»òÕßtot[0]µÄ³¤¶È£¬ÓÃÓÚÈ¡µ½ÕýÈ·µÄtotÖµ£¬ÒòÎªtot[1][0]ºÍtot[0][0]µÄµØÖ·Ïà¸ônTotMax
//IsSplit£ºÊÇ·ñ½«RPCµÄ·ù¶È·Ö³Éavalanche ºÍ streamer ·Ö±ð½ÃÕýÔÙ·ÅÔÚÒ»ÆðfitµÃµ½Ê±¼ä·Ö±æÂÊ£¬0--·ñ£¬1--ÊÇ£¬Ä¬ÈÏÎª0£¬Èç¹û²»ÐèÒªÎÞÐèÌî
//totSplit£ºavalancheºÍstreamer·ù¶ÈµÄ·Ö½çµã£¬Ä¬ÈÏÎª0£¬Èç¹û²»ÐèÒª¾ÍÎÞÐèÌî
//IdxRpcCharge£ºRPCµÄµçºÉÔÚtot¶þÎ¬Êý×éÉÏµÄÐòºÅ£¬¼´RPCµçºÉ´Ótot[IdxRpcCharge][0]¿ªÊ¼£¬Ä¬ÈÏÎª4£¬Èç¹û²»ÐèÒª¾ÍÎÞÐèÌî
******************************************************************************************************************************/
void DoCorrection(float *ttm,float *tot,int nCount,int nCharge,int nCor,char *psName,float tLimit,int tBins,float *totMin,float *totMax,int totBins,float *fitMin,float *fitMax,TF1 *fitf,float Nsigma,int nTotMax,int IsSplit=0,float totSplit=0,int IdxRpcCharge=4)
{
	cout<<"######## "<<psName<<" Calibration ########"<<"\n";
	cout<<"Good events:    "<<nCount<<"\n";
	
	double mean,sigma,para[16],para2[16],err;
	TH1D *T;
	TH2D *Corr;
	TH1D *htemp;
	int nPad;
	double cor;
	char buf[1024];
  
	TF1 *f1 = new TF1("f1","gaus");
	TF1 *ftemp = new TF1("ftemp","[0]*exp(-0.5*((x-[1])/[2])**2)",-tLimit[0],tLimit[0]);
	f1->SetLineColor(2);
	ftemp->SetLineWidth(3);
	ftemp->SetLineStyle(2);
	ftemp->SetLineColor(2);
	
	sprintf(buf,"%s\/%s_Calibration.ps",directory,psName);
	TPostScript *ps = new TPostScript(buf,112);
	ps->Range(24,16);
  
	TCanvas *c1 = new TCanvas("c1","c1",0,0,1000,800);
	c1->SetFillColor(kWhite);
	c1->GetFrame()->SetFillColor(21);
 	c1->GetFrame()->SetBorderSize(6);
	c1->GetFrame()->SetBorderMode(-1);
  
	//gStyle->SetOptFit(110);
	//gStyle->SetOptStat(110);
	gStyle->SetOptFit(0);
	gStyle->SetOptStat(0);
	
	//////////////output fit result
  ofstream outResult;
  sprintf(buf,"%s\/%s_FitResult.txt",directory,fileName);
  outResult.open(buf,ios::app);
  outResult<<"########    "<<fileName<<"    #######"<<"\n";
  outResult<<"nSigma	nCor	nBinsPeakT	nBinsPeakTot	minPeakTot[0]	maxPeakTot[0]	fitMinTot[0]	fitMaxTot[0]	halfPeakWidth"<<"\n";
  outResult<<nSigma<<'	'<<nCor<<'	'<<nBinsPeakT<<'	'<<nBinsPeakTot<<'	'<<minPeakTot[0]<<'	'<<maxPeakTot[0]<<'	'<<fitMinTot[0]<<'	'<<fitMaxTot[0]<<'	'<<halfPeakWidth<<"\n";
  outResult<<"nCorr	Sigma	Err of Sigma"<<"\n";
  /////////////////////////
  
  for(int i=0;i<nCor;i++)
  {
  	cout<<"########### The No."<<i+1<<" correction. ##########"<<"\n";
  	ps->NewPage();
  	c1->Clear();
  	c1->Divide(1,1);
  	c1->cd(1)->Update();
  	//c1->cd(1)->SetLogy();
  	
  	sprintf(buf,"%s_%d_Original",psName,i);
//		sprintf(buf," ");
  	T= new TH1D(buf,buf,tBins,-tLimit,tLimit); 
	  T->SetLineWidth(2);
	  
  	for(int m=0;m<nCount;m++)
  	{
  		T->Fill(ttm[m]);
  	}
  	mean=T->GetMean();
  	sigma=T->GetRMS();
	  f1->SetRange(mean-Nsigma*sigma, mean+Nsigma*sigma);
	  T->Fit("f1","NQR");
	  f1->GetParameters(&para[0]);
	  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
	  T->Fit("f1","NQR");
	  f1->GetParameters(&para[0]);
	  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
	  T->Draw();
	  T->Fit("f1","QR");
	  f1->GetParameters(&para[0]);
		err=f1->GetParError(2);
	  
	  ftemp->SetParameters(para[0],para[1],para[2]);
	  ftemp->Draw("same");
	  
//	  T->GetXaxis()->SetTitle("((PMT1+PMT2)-(PMT3+PMT4))/4 [¡Á35ps]");
		T->GetXaxis()->SetTitle("t [ns]");
		T->GetXaxis()->CenterTitle();
		T->GetYaxis()->SetTitle("Counts");
		T->GetYaxis()->CenterTitle();
		/////////////////////
		if(i==0)
			{
  			outResult<<"No.0"<<'\t'<<para[2]*1000.0<<'\t'<<err*1000.0<<"\n";
  		}
  	////////////////////
		
		c1->cd(1)->Update();
		
  	for(int j=0;j<nCharge;j++)
  	{
  		if(j%2==0) 
  		{
  			ps->NewPage();
  			c1->Clear();
  			c1->Divide(2,2);
  			nPad=1;
  		}
  		else 
  		{
  			nPad=2;
  		}
  		sprintf(buf,"%s_Correction(%d)_%d",psName,i,j);
  		Corr = new TH2D(buf,buf,totBins,totMin[j],totMax[j],tBins,-tLimit,tLimit);
  		
  		for(m=0;m<nCount;m++)
  		{
  			Corr->Fill(*(tot+j*nTotMax+m),ttm[m]);
  		}
  		c1->cd(nPad);
  		Corr->Draw("colz");
  			
  		Corr->GetXaxis()->SetTitle("TOT [ns]");
			Corr->GetXaxis()->CenterTitle();
			Corr->GetYaxis()->SetTitle("t [ns]");
			Corr->GetYaxis()->CenterTitle();
			c1->cd(nPad)->Update();
  		
  		Corr->ProfileX();
			sprintf(buf,"%s_Correction(%d)_%d_pfx",psName,i,j);
			htemp=(TProfile *)gDirectory->Get(buf);
			c1->cd(nPad)->Update();
//			htemp->Draw("QRsame");
//			htemp->SetLineColor(kBlack);
//			htemp->SetMarkerColor(kBlack);
//			htemp->SetMarkerStyle(21);
//			htemp->SetMarkerSize(0.6);
//  		c1->cd(nPad)->Update();
  		
  		if(IsSplit==1&&j==IdxRpcCharge)
  		{
  			TF1 *fitf2= new TF1("fitf2","fitf");
	  		htemp->Fit("fitf","q","QRsame",fitMin[j],totSplit);
	  		fitf->GetParameters(&para[0]);
	  		fitf2->SetRange(fitMin[j],totSplit);
	  		fitf2->SetParameters(&para[0]);
				fitf2->SetLineWidth(2);
	  		fitf2->Draw("same");
  			c1->cd(nPad)->Update();
	  		htemp->Fit("fitf","q","QRsame",totSplit,fitMax[j]);
	  		fitf->GetParameters(&para2[0]);
  			c1->cd(nPad)->Update();	  		
	  		for(m=0;m<nCount;m++)
	  		{
	  			if((*(tot+j*nTotMax+m))<totSplit)
	  			{
//		  			cor=para[0]+para[1]/sqrt((*(tot+j*nTotMax+m)))+para[2]/(*(tot+j*nTotMax+m))+para[3]*(*(tot+j*nTotMax+m));
						cor=0;
		  			ttm[m]=ttm[m]-cor;
	  			}
	  			if((*(tot+j*nTotMax+m))>totSplit)
	  			{
//		  			cor=para2[0]+para2[1]/sqrt((*(tot+j*nTotMax+m)))+para2[2]/(*(tot+j*nTotMax+m))+para2[3]*(*(tot+j*nTotMax+m));
						cor=0;
		  			ttm[m]=ttm[m]-cor;
	  			}
	  		} 				
  		}
  		else
  		{
	  	//	htemp->Fit("fitf","q","QRsame",fitMin[j],fitMax[j]);
	  		htemp->Fit("fitf","Nq","QRsame",fitMin[j],fitMax[j]);
	  		c1->cd(nPad)->Update();
	  		fitf->GetParameters(&para[0]);
	  		for(m=0;m<nCount;m++)
	  		{
//	  			cor=para[0]+para[1]/sqrt((*(tot+j*nTotMax+m)))+para[2]/(*(tot+j*nTotMax+m))+para[3]*(*(tot+j*nTotMax+m));
	  			cor=fitf->Eval((*(tot+j*nTotMax+m)));
	  			ttm[m]=ttm[m]-cor;
	  		}
  		}
  		
//  		sprintf(buf,"%s-T(%d)_%d",psName,i,j);
  		sprintf(buf," ");
  		T= new TH1D(buf,buf,tBins,-tLimit,tLimit); 
  		for(m=0;m<nCount;m++)
  		{
  			T->Fill(ttm[m]);
  		}
  		c1->cd(nPad+2);
//  		c1->cd(nPad+2)->SetLogy();
			c1->cd(nPad+2)->SetBorderSize(6);
	  	mean=T->GetMean();
	  	sigma=T->GetRMS();
		  f1->SetRange(mean-Nsigma*sigma, mean+Nsigma*sigma);
  		T->Fit("f1","NQR");
		  f1->GetParameters(&para[0]);
		  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
		  T->Fit("f1","NQR");
		  f1->GetParameters(&para[0]);
		  f1->SetRange(para[1]-Nsigma*para[2], para[1]+Nsigma*para[2]);
		  T->Draw();
		  T->Fit("f1","QR");
		  f1->GetParameters(&para[0]);
		  err=f1->GetParError(2);
		  
		  ftemp->SetParameters(para[0],para[1],para[2]);
		  ftemp->Draw("same");
		  
		  T->GetXaxis()->SetTitle("t [ns]");
//		  T->GetXaxis()->SetTitle("t [\\times 35ps]");
//		  T->GetXaxis()->SetTitle("((PMT1+PMT2)-(PMT3+PMT4))/4 [\\times 35ps]");
			T->GetXaxis()->CenterTitle();
			T->GetYaxis()->SetTitle("Counts");
			T->GetYaxis()->CenterTitle();
  		c1->cd(nPad+2)->Update();
  		
  		/////////////////////
  		outResult<<"No."<<i+1<<'\t'<<para[2]*1000.0<<'\t'<<err*1000.0<<"\n";
  		////////////////////
  	}
  }
  
  outResult.close();
	ps->Close();
	delete c1;
}

/************  Read .dat file and convert to .root file***************************************
*********************************************************************************************/
void ReadFile(char *file1)
{
	cout<<"############## nMax = "<<nMax<<"\n";
	int chIn,tsIn,edgeIn,dt;
	float tsBuf,tsLaserBuf,totLaserBuf,tsSticBuf,tSticBuf,tsSticBuf2,totSticBuf;
	char buf1[1024];
	float tStic,totStic,totTrigger;
	cout<<"####	Read File	"<<file1<<"	####"<<"\n";
	ifstream is;
	is.open(file1,ios::in);
	nEvents = 0;

	tsBuf=0;
	tsLaserBuf = 0;
	tsSticBuf  = 0;
	tSticBuf   = 0;
	tsSticBuf2 = 0;
	totSticBuf = 0;

	int i = 0;
	while(is>>chIn>>tsIn>>edgeIn)
	{
		if(tsIn<tsBuf)	{
			tsLaserBuf = 0;
			tsSticBuf  = 0;
			tSticBuf   = 0;
			tsSticBuf2 = 0;
			totSticBuf = 0;
		}	
		tsBuf=tsIn;
		if(i%100000 == 0)	cout<<"!!!!    Reading Data	"<<i<<"!!!!	nEvents = "<<nEvents<<"\n";
		i++;
		switch(chIn)	{
			case chLaser:	{
				if(edgeIn == edgeLaser[0])	tsLaserBuf = tsIn;
				break;
			}
			case chStic:		{
				if(edgeIn == edgeStic[0])	{
					if(totSticBuf == 0)	{
						tsSticBuf = tsIn;
						tSticBuf = tsSticBuf - tsLaserBuf;
					}
					else	{	
						dt = tsIn-tsSticBuf2;
						if(dt>totHyst)	{
							if(tSticBuf<maxDelay)	{
								tot[nEvents] = totSticBuf*binTDC/1000.0;	//ns
								t[nEvents] = tSticBuf*binTDC/1000.0;	//ns
								nEvents++;
							}
							totSticBuf = 0;
							tsSticBuf = tsIn;
							tSticBuf = tsSticBuf - tsLaserBuf;
						}
					}
				}
				else if(tsSticBuf != 0)	{
					totSticBuf = tsIn-tsSticBuf;
					tsSticBuf2 = tsIn;
				}	
				break;		
			}
			default:	{
		//		cout<<chIn<<"	!!!!!	switch channel == default	!!!!!!!!!"<<"\n";
				break;
			}
		}	
	}
	is.close();
	
	cout<<"####	Complete reading file :"<<file1<<" !	####"<<"\n";
	cout<<"####	Start to create root file!	####"<<""<<"\n";

	//create root file
	sprintf(buf1,"%s\/%s.root",directory,fileName);
	TFile *rootFile1 = new TFile(buf1,"RECREATE","Root file from TDC data");
	TTree *h1 = new TTree("h1","Data from TDC");
	h1->Branch("t",&tStic,"t/F");
	h1->Branch("tot",&totStic,"tot/F");
	for(i=0;i<nEvents;i++)  {
		tStic = t[i];
		totStic = tot[i];
		h1->Fill();
	}
	h1->Print();
	h1->Write();

	//close the root file
	rootFile1->Close();
	cout<<"####	Function ReadFile complete	####"<<"\n";
}

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
	//	else continue;
	}
	chan->update_histos();
	//printf("!!!!!!!!coming here!!\n");
	is.close();
}
/************  Read .dat file and convert to .root file***************************************
*************  correction with class tdc_ch_value     ****************************************
*********************************************************************************************/
void ReadFile(char *file1, tdc_ch_values *chanLaser, tdc_ch_values *chanStic)
{
	printf("\n####	Start to read file with TDC NL correction!	####\n");
	int chIn,tsIn,edgeIn,dt;
	float tsBuf,tsLaserBuf,totLaserBuf,tsSticBuf,tSticBuf,tsSticBuf2,totSticBuf;
	char buf1[1024];
	float tStic,totStic,totTrigger;

	chanLaser = new tdc_ch_values(1,1024);
	chanStic = new tdc_ch_values(1,1024);
	CreatTdcCorrection(file1,chanLaser,chLaser,edgeLaser[0]); 
	CreatTdcCorrection(file1,chanStic,chStic,edgeStic[0]); 
	chanLaser->update_histos();
	chanStic->update_histos();
	printf("----> finish create TDC NL correction!	<---- \n");

	cout<<"####	Read File	"<<file1<<"	####"<<"\n";
	ifstream is;
	is.open(file1,ios::in);

	is.clear();
	is.seekg(0,ios::beg);
	nEvents = 0;

	tsBuf=0;
	tsLaserBuf = 0;
	tsSticBuf  = 0;
	tSticBuf   = 0;
	tsSticBuf2 = 0;
	totSticBuf = 0;
	int i = 0;
	while(is>>chIn>>tsIn>>edgeIn)
	{
		if(tsIn<tsBuf)	{
			tsLaserBuf = 0;
			tsSticBuf  = 0;
			tSticBuf   = 0;
			tsSticBuf2 = 0;
			totSticBuf = 0;
		}	
		tsBuf=tsIn;
		if(i%100000 == 0)	cout<<"!!!!    Reading Data	"<<i<<"!!!!	nEvents = "<<nEvents<<"\n";
		i++;
		switch(chIn)	{
			case chLaser:	{
				if(edgeIn == edgeLaser[0])	tsLaserBuf = chanLaser->get_dith_value(tsIn);
				break;
			}
			case chStic:		{
				if(edgeIn == edgeStic[0])	{
					if(totSticBuf == 0)	{
						tsSticBuf = chanStic->get_dith_value(tsIn);
						tSticBuf = tsSticBuf - tsLaserBuf;
					}
					else	{	
						dt = chanStic->get_dith_value(tsIn)-tsSticBuf2;
						if(dt>totHyst)	{
							if(tSticBuf<maxDelay)	{
								tot[nEvents] = totSticBuf*binTDC/1000.0;	//ns
								t[nEvents] = tSticBuf*binTDC/1000.0;	//ns
								nEvents++;
							}
							totSticBuf = 0;
							tsSticBuf = chanStic->get_dith_value(tsIn);
							tSticBuf = tsSticBuf - tsLaserBuf;
						}
					}
				}
				else if(tsSticBuf != 0)	{
					totSticBuf = chanStic->get_dith_value(tsIn)-tsSticBuf;
					tsSticBuf2 = chanStic->get_dith_value(tsIn);
				}	
				break;		
			}
			default:	{
		//		cout<<chIn<<"	!!!!!	switch channel == default	!!!!!!!!!"<<"\n";
				break;
			}
		}	
	}
	is.close();
	
	cout<<"####	Complete reading file :"<<file1<<" !	####"<<"\n";
	cout<<"####	Start to create root file!	####"<<""<<"\n";

	//create root file
	sprintf(buf1,"%s\/%s.root",directory,fileName);
	TFile *rootFile1 = new TFile(buf1,"RECREATE","Root file from TDC data");
	TTree *h1 = new TTree("h1","Data from TDC");
	h1->Branch("t",&tStic,"t/F");
	h1->Branch("tot",&totStic,"tot/F");
	for(i=0;i<nEvents;i++)  {
		tStic = t[i];
		totStic = tot[i];
		h1->Fill();
	}
	h1->Print();
	h1->Write();

	//close the root file
	rootFile1->Close();
	cout<<"####	Function ReadFile complete	####"<<"\n";
}

/************  Read .root file and put the data into variable*********************************
*********************************************************************************************/
void ReadRootFile(TFile *rootFile,char *treeName)
{
	float tBuf,totBuf;
 	TTree *h1 = rootFile->Get(treeName);
	h1->Print();
   	h1->SetBranchAddress("t",&tBuf);
   	h1->SetBranchAddress("tot",&totBuf);
	if(modeSelect == 1)	nEvents = nSelect;
	else	nEvents = (int)h1->GetEntries();
	for(int i=0;i<nEvents;i++)
	{
		h1->GetEntry(i);
		t[i]	= tBuf;
		tot[i]	= totBuf;
	}
}

/************ Plot the status of the data!  ************************************************** 
*********************************************************************************************/
void PlotStatus()
{
	cout<<"\n"<<"####	Plot the spectrums	####"<<"\n";
	int i,j,k;
	char buf2[1024];
	TObjArray histlist(0);
	TH1D *histT = new TH1D("t","t",nBinsAllT,minT,maxT);
	TH1D *histTOT = new TH1D("TOT","TOT",nBinsAllTot,minTot,maxTot);
	TH2D *tVsTot = new TH2D("tVsTot","t vs. TOT",nBinsPeakTot,minPeakTot[0],maxPeakTot[0],nBinsPeakT,-shift_t-halfPeakWidth,-shift_t+halfPeakWidth);
	TGraph *sigmaSlice; 
	histlist.Add(histT);
	histlist.Add(histTOT);
	histlist.Add(sigmaSlice);
	histlist.Add(tVsTot);
	
//	cout<<"\n"<<"nEvents = "<<nEvents<<"\n"<<"\n";
	for(i=0;i<nEvents;i++)
	{
		histT->Fill(t[i]);
		histTOT->Fill(tot[i]);
//		cout<<t[i]<<'\t'<<tot[i]<<"\n";
	}
	
	TCanvas *cStatus1 = new TCanvas("status1","status1",10,10,800,600);
//  cStatus->Clear();
//	cStatus->Divide(2,1);
//	cStatus->cd(1);
	histT->Draw();
	histT->GetXaxis()->SetTitle("delay [ns]");
	histT->GetXaxis()->CenterTitle();
	histT->GetYaxis()->SetTitle("Counts");
	histT->GetYaxis()->CenterTitle();
	cStatus1->Update();
	
	TCanvas *cStatus2 = new TCanvas("status2","status2",10,10,800,600);
//	cStatus->cd(2);
	histTOT->Draw();
	histTOT->GetXaxis()->SetTitle("TOT [ns]");
	histTOT->GetXaxis()->CenterTitle();
	histTOT->GetYaxis()->SetTitle("Counts");
	histTOT->GetYaxis()->CenterTitle();
	cStatus2->Update();
	double maxY;
	maxY = gPad->GetUymax(); 
//	cout<<"max of Y axis = "<<maxY<<"\n";
	TLine *l1 = new TLine(minPeakTot[0],0,minPeakTot[0],maxY);
	l1->SetLineColor(kRed);
	l1->Draw();
	TLine *l2 = new TLine(maxPeakTot[0],0,maxPeakTot[0],maxY);
	l2->SetLineColor(kRed);
	l2->Draw();
	cStatus2->Update();
	

	TF1 *f1 = new TF1("f1","gaus");
	const int maxNSlice = 1024;
	TH1D *histSlice[maxNSlice];
	double para[16],mean,sigma,err;
	double binSliceLow[maxNSlice],binSliceUp[maxNSlice],meanBinSlice[maxNSlice],meanT[maxNSlice],sigmaT[maxNSlice],sigmaSigmaT[maxNSlice],nEventsSlice[maxNSlice];
	for(i=0;i<nBinsPeakTot;i++)	{
		sprintf(buf2,"tSlice_%d",i);
		histSlice[i] = new TH1D(buf2,buf2,int((maxT-minT)/binningPeakT),minT,maxT);
		binSliceLow[i]=minPeakTot[0]+i*binningPeakTot;
		binSliceUp[i]=minPeakTot[0]+(i+1)*binningPeakTot;
		meanT[i]=0;
		sigmaT[i]=0;
		sigmaSigmaT[i]=0;
		nEventsSlice[i]=0;
	}
	for(i=0;i<nEvents;i++)	{
//		if(i%100000==0)	cout<<"!!!!!!!!!!!!!	"<<i<<"	!!!!!!!!!!!!!!!!!!!"<<"\n";
		for(int j=0;j<nBinsPeakTot;j++)	{
			if(tot[i]>=binSliceLow[j] && tot[i]<binSliceUp[j])	{
				histSlice[j]->Fill(t[i]);
				nEventsSlice[j]++;
			}
		}
	}
	ofstream outProfile;
	sprintf(buf2,"%s\/%s_SliceProfile.txt",directory,fileName);
	outProfile.open(buf2,ios::trunc);
	outProfile<<"Slice#	TOTLow	TOTUP	nEvents	meanT	sigmaT	sigmaSigmaT"<<"\n";
	for(i=0;i<nBinsPeakTot;i++)	{
		mean=histSlice[i]->GetMean();
		sigma=histSlice[i]->GetRMS();
		f1->SetRange(mean-nSigma*sigma, mean+nSigma*sigma);
		histSlice[i]->Fit("f1","NQR");
		f1->GetParameters(&para[0]);
		f1->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
		histSlice[i]->Fit("f1","NQR");
		f1->GetParameters(&para[0]);
		f1->SetRange(para[1]-nSigma*para[2], para[1]+nSigma*para[2]);
		histSlice[i]->Fit("f1","NQR");
		f1->GetParameters(&para[0]);
		err=f1->GetParError(2);
		meanT[i]=para[1];
		sigmaT[i]=para[2]*1000.0;
		sigmaSigmaT[i]=err*1000.0;
		meanBinSlice[i]=0.5*(binSliceLow[i]+binSliceUp[i]);
		outProfile<<i<<"	"<<binSliceLow[i]<<"	"<<binSliceUp[i]<<"	"<<nEventsSlice[i]<<"	"<<para[1]<<"	"<<para[2]*1000.0<<"	"<<err*1000.0<<"\n";
	}
	outProfile.close();
	TCanvas *cStatus3 = new TCanvas("status3","status3",10,10,800,600);
	sigmaSlice = new TGraph(nBinsPeakTot,meanBinSlice,sigmaT);
	sigmaSlice->SetTitle("sigmaT vs tot");
	sigmaSlice->Draw("A*");
	cStatus3->Update();

	for(i=0;i<nEvents;i++)	{
		if(tot[i]>minPeakTot[0] && tot[i]<maxPeakTot[0])	{	
			tVsTot->Fill(tot[i],t[i]);
		}
	}
	TCanvas *cStatus4 = new TCanvas("status4","status4",10,10,800,600);
	tVsTot->Draw("colz");
  	tVsTot->GetXaxis()->SetTitle("TOT [ns]");
	tVsTot->GetXaxis()->CenterTitle();
	tVsTot->GetYaxis()->SetTitle("t [ns]");
	tVsTot->GetYaxis()->CenterTitle();

	TProfile *tTotProfile = tVsTot->ProfileX();
	tTotProfile->Draw("QRsame");
	tTotProfile->SetMarkerStyle(21);
	tTotProfile->SetMarkerSize(0.6);
	cStatus4->Update();

	sprintf(buf2,"%s\/%s_RunStatus.root",directory,fileName);
	TFile *hf = new TFile(buf2,"recreate");
	histlist.Write();
	cout<<"####	Complete ploting spectrums !	####"<<"\n"<<"\n";
}

/************  Selecte the data and do the  calibration. *************************************
*********************************************************************************************/
void TQCorrection()
{
	cout<<"####	Start T-Q correction!	####"<<"\n";
	int i;
	int nPeakCounts=0;
	float t2[nMax],tot2[nMax];
	char psName[1024];
	sprintf(psName,"%s",fileName);
	
	for(i=0;i<nEvents;i++)
	{
//		cout<<tot[i]<<"\n";
		if(tot[i]>minPeakTot[0] && tot[i]<maxPeakTot[0])
			{
				t2[nPeakCounts]=t[i]+shift_t;
				tot2[nPeakCounts]=tot[i];
				nPeakCounts++;
			}
	}

	cout<<"	Number of good events = "<<nPeakCounts<<"\n";
	sprintf(psName,"TDC");
	DoCorrection(t2,tot2,nPeakCounts,nCharge,nCor,psName,halfPeakWidth,nBinsPeakT,minPeakTot,maxPeakTot,nBinsPeakTot,fitMinTot,fitMaxTot,fitf,nSigma,nPeakCounts);
	cout<<"####	Complete T-Q correction!	####"<<"\n"<<"\n";
}

/************ Main funtion !!  *************************************************************** 
*********************************************************************************************/
void readTDC()
{
	char buf[1024];
	char modeIn;
	cout<<"*******************************************************************"<<"\n";
	cout<<"*******************************************************************"<<"\n";
	cout<<"The mode of this macro:"<<"\n";
	cout<<"1. Read .dat file and convert to .root file."<<"\n";
	cout<<"2. Read .root file and plot the t and TOT spectrum."<<"\n";
	cout<<"3. Read .root file and do TOT-t calibration."<<"\n";
	cout<<"4. Read .dat file, convert to .root file and do TOT-t calibration."<<"\n";
	cout<<"*******************************************************************"<<"\n";
	cout<<"Please select the mode: ";
	cin>>modeIn;
	cout<<"\n";
	switch(modeIn)	{
		case '1':	{
			sprintf(buf,"%s\/%s.dat",directory,fileName);	
		
			ifstream isTest;
			isTest.open(buf,ios::in);
			if(!isTest) {
				cout<<"\n"<<buf<<" doesn't exist!!"<<"\n"<<"\n";
				return;
			}
			else {
				isTest.close();
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				tdc_ch_values *chanLaser,*chanStic;
				ReadFile(buf,chanLaser,chanStic);
				PlotStatus();
			}
			break;
		}
		case '2':	{
			sprintf(buf,"%s\/%s.root",directory,fileName);	
		
			TFile *f1 = new TFile(buf);
			if(f1->IsZombie())	{
				cout<<"Error openning file!!"<<"\n";
				return;
			}
			else {
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				ReadRootFile(f1,"h1");
				PlotStatus();
			}
			break;
		}
		case '3':	{
			sprintf(buf,"%s\/%s.root",directory,fileName);	
		
			TFile *f1 = new TFile(buf);
			if(f1->IsZombie())	{
				cout<<"Error openning file!!"<<"\n";
				return;
			}
			else {
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				ReadRootFile(f1,"h1");
			//	PlotStatus();
				TQCorrection();
			//	f1->Close();
			}
			break;
		}
		case '4':	{
			sprintf(buf,"%s\/%s.dat",directory,fileName);	
			ifstream isTest;
			isTest.open(buf,ios::in);
			if(!isTest) {
				cout<<"\n"<<buf<<" doesn't exist!!"<<"\n"<<"\n";
				return;
			}
			else {
				isTest.close();
				cout<<"File:"<<"\n"<<"########  "<<buf<<"  ########"<<"\n"<<"\n";
				ReadFile(buf);
				PlotStatus();
				TQCorrection();
			}
			break;
		}
		default:	{
			cout<<"Wrong Mode!!!!"<<"\n";
			break;
		}
	}
}
